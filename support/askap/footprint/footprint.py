#!/usr/bin/env python
# Copyright (c) 2015-2019 CSIRO
# Australia Telescope National Facility (ATNF)
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# PO Box 76, Epping NSW 1710, Australia
# atnf-enquiries@csiro.au
#
# This file is part of the ASKAP software distribution.
#
# The ASKAP software distribution is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the License
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
#
from __future__ import division
from builtins import map
from builtins import zip
from builtins import object
import numpy as np
import math
import re
import glob
import os
#from askap.parset import ParameterSet
#from askap import logging
#logger = logging.getLogger(__name__)
try:
    from askap import logging
except ImportError:
    import logging


#__all__ = ["Footprint", "FootprintFactory", "FootprintParsetDict", "Skypos", "from_parset", "ras_rad", 
#           "decs_rad", "dp_to_lm", "lm_to_dp"]
__all__ = ["Footprint",  "Skypos", "from_parset", "ras_rad", 
           "decs_rad", "dp_to_lm", "lm_to_dp"]

RAD2DEG = 180.0/math.pi
DEG2RAD = math.pi/180.0


class Skypos:
    """Defines a class that works with spherical geometry, specifically points
    in a unit sphere, such as the sky.

    This is general spherical geometry, with little to tie it to astronomy. The
    exceptions are the naming of longitude and latitude as RA,Dec
    """
    def_precra = 3
    def_precde = 2

    def __init__(self, ra, dec, precra=def_precra, precdec=def_precde):
        """
        Initialise a Skypos object defining a point on a unit sphere with longitude ra and latitude dec
        :param ra: right ascension (radians or hh:mm:ss.ss)
        :type ra: float or str
        :param dec: declination (radians or dd:mm:ss.ss)
        :type dec: float or str
        :param precra:
        :param precdec:
        """
        if isinstance(ra, str):
            self.ra = ras_rad(ra)
            self.dec = decs_rad(dec)
        else:
            self.ra = ra
            self.dec = dec
        self.precra = precra
        self.precdec = precdec
        self.rn = 12+self.precra-Skypos.def_precra
        self.dn = 12+self.precdec-Skypos.def_precde
        self.ras = None
        self.decs = None
        ps = math.pi * 0.5 - self.dec
        sps = math.sin(ps)
        cps = math.cos(ps)
        sra = math.sin(self.ra)
        cra = math.cos(self.ra)
        self._dvecx = [cps * cra, cps * sra, -sps]
        self._dvecy = [-sra, cra, 0.0]
        self._dvecz = [sps * cra, sps * sra, cps]
        self._vec = [cra*sps, sra*sps, cps]

    def get_vec(self):
        return self._vec
        
    def d_pa(self, other):
        ra = other.ra
        dec = other.dec
        xyz = rd_xyz(ra, dec)
        d = None
        x = self._dvecx[0] * xyz[0] + self._dvecx[1] * xyz[1] \
            + self._dvecx[2] * xyz[2]
        y = self._dvecy[0] * xyz[0] + self._dvecy[1] * xyz[1] \
            + self._dvecy[2] * xyz[2]
        z = self._dvecz[0] * xyz[0] + self._dvecz[1] * xyz[1] \
            + self._dvecz[2] * xyz[2]
        z = max(-1.0, min(z, 1.0))
        try:
            d = math.pi * 0.5 - math.asin(z)
        except ValueError:
            print("Can't perform asin(z) on z = {}, math domain error".format(z))
        pa = math.atan2(y, -x)
        return d, pa

    def offset(self, dpa):
        a2 = math.pi * 0.5 - dpa[0]
        a1 = -(math.pi + dpa[1])
        xyz = rd_xyz(a1, a2)
        x = self._dvecx[0] * xyz[0] + self._dvecy[0] * xyz[1] + \
            self._dvecz[0] * xyz[2]
        y = self._dvecx[1] * xyz[0] + self._dvecy[1] * xyz[1] + \
            self._dvecz[1] * xyz[2]
        z = self._dvecx[2] * xyz[0] + self._dvecy[2] * xyz[1] + \
            self._dvecz[2] * xyz[2]
        b2 = math.asin(z)
        b1 = (2.0 * math.pi + math.atan2(y, x)) % (2.0 * math.pi)
        return Skypos(b1, b2)

    # noinspection PyUnresolvedReferences
    def separation(self, other):
        """
        Return great circle angular separation between this Skypos and another
        :param other: position on the sky to determine separation to
        :type other: `:class:Skypos`
        :return: separation angle (radians)
        """

        if self.dec == other.dec:
            if self.ra == other.ra:
                return 0.

        # vincenty formula (https://en.wikipedia.org/wiki/Great-circle_distance)
        dra = np.abs(self.ra-other.ra)
        sep = np.arctan2(np.sqrt((np.cos(other.dec)*np.sin(dra))**2
                         + (np.cos(self.dec)*np.sin(other.dec)
                         - np.sin(self.dec)*np.cos(other.dec)*np.cos(dra))**2),
                         np.sin(self.dec)*np.sin(other.dec) + np.cos(self.dec)*np.cos(other.dec)*np.cos(dra))

        return sep

    def rotate_x(self, a):
        """return a skypos determined by rotating self about the X-axis by 
        angle a."""
        x, y, z = _rotate_v_x(self._vec, a)
        b2 = math.asin(z)
        b1 = (2 * math.pi + math.atan2(y, x)) % (2.0 * math.pi)
        return Skypos(b1, b2)

    def rotate_y(self, a):
        """return a skypos determined by rotating self about the X-axis by 
        angle a."""
        x, y, z = _rotatev_y(self._vec, a)
        b2 = math.asin(z)
        b1 = (2 * math.pi + math.atan2(y, x)) % (2.0 * math.pi)
        return Skypos(b1, b2)

    def rotate_z(self, a):
        """return a skypos determined by rotating self about the X-axis by 
        angle a."""
        x, y, z = _rotate_v_z(self._vec, a)
        b2 = math.asin(z)
        b1 = (2 * math.pi + math.atan2(y, x)) % (2.0 * math.pi)
        return Skypos(b1, b2)

    def shift(self, delta_lon, delta_lat):
        """
        Shift this direction (Skypos) in longitude and latitude.
        The longitude shift will be in radian units perpendicular to the direction to pole, along a great circle.
 
        :param float delta_lon: longitude (RA) offset in radians
        :param float delta_lat: latitude (DEC) offset in radians
        """
        lat = self.dec
        lon = self.ra
        # vector along X axis (first point of Aries)
        x0 = Skypos('0h0m0s', '0:0:0', 3, 3)
        shifted_direction = x0.rotate_z(delta_lon).rotate_y(lat + delta_lat).rotate_z(lon)
        return shifted_direction

    def get_ras(self):
        if self.ras is None:
            self.ras = ras(self.ra)
            self.decs = decs(self.dec)
        return self.ras[:self.rn]

    def get_decs(self):
        if self.ras is None:
            self.ras = ras(self.ra)
            self.decs = decs(self.dec)
        return self.decs[:self.dn]

    def __str__(self):
        return '{} {}'.format(self.get_ras(), self.get_decs())


def ras(ra):
    s = ra * (4.0 * 60.0 * RAD2DEG)
    hh = int(s / 3600.0)
    mm = int(s / 60.0) - hh * 60
    ss = s - 60 * (mm + 60 * hh)
    if "{:9.6f}".format(ss) == '60.000000':
        ss = 0.0
        mm += 1
        if mm == 60:
            mm = 0
            hh += 1
            if hh == 24:
                hh = 0
    return "%02d:%02d:%09.6f" % (hh, mm, ss)


def decs(dec):
    s = abs(dec) * (60.0 * 60.0 * RAD2DEG)
    dd = int(s / 3600.0)
    mm = int(s / 60.0) - dd * 60
    ss = s - 60 * (mm + 60 * dd)
    if "%8.5f" % ss == '60.00000':
        ss = 0.0
        mm += 1
        if mm == 60:
            mm = 0
            dd += 1
    sign = ' '
    if dec < 0.0:
        sign = '-'
    return "%s%02d:%02d:%08.6f" % (sign, dd, mm, ss)


def _rotate_v_x(vec, a):
    """Return a skypos determined by rotating vec about the X-axis by 
    angle a."""
    ca, sa = math.cos(a), math.sin(a)
    x = vec[0]
    y = vec[1] * ca - vec[2] * sa
    z = vec[1] * sa + vec[2] * ca
    return [x, y, z]


def _rotatev_y(vec, a):
    """Return a skypos determined by rotating vec about the Y-axis by 
    angle a."""
    ca, sa = math.cos(a), math.sin(a)
    x = vec[0] * ca - vec[2] * sa
    y = vec[1]
    z = vec[0] * sa + vec[2] * ca
    return [x, y, z]


def _rotate_v_z(vec, a):
    """Return a skypos determined by rotating vec about the Z-axis by 
    angle a."""
    ca, sa = math.cos(a), math.sin(a)
    x = vec[0] * ca - vec[1] * sa
    y = vec[0] * sa + vec[1] * ca
    z = vec[2]
    return [x, y, z]


def ras_rad(ra_string):
    """
    Convert right ascension string to radians
    :param ra_string: right ascension string (hh:mm:ss.ss)
    :type ra_string: str
    :return: right ascension in radians
    :rtype: float
    """
    if ra_string[0] == '-':
        raise ValueError('Right ascension may not be negative: {}'.format(ra_string))
    (a, b, c) = re.findall("[0-9.]+", ra_string)
    hh, mm = list(map(int, [a, b]))
    ss = float(c)
    return (ss + 60.0 * (mm + 60.0 * hh)) * 2.0 * math.pi / 86400.0


def decs_rad(dec_string):
    """
    Convert declination string to radians
    :param dec_string: declination string (dd:mm:ss.ss)
    :type dec_string: str
    :return: declination in radians
    :rtype: float
    """
    a, b, c = re.findall('[0-9.]+', dec_string)
    dd, mm = list(map(int, [a, b]))
    ss = float(c)
    r = (ss + 60.0 * (mm + 60.0 * dd)) * 2.0 * math.pi / 1296000.0
    if dec_string[0] == '-':
        r = -r
    return r


def dp_to_lm(dp):
    """
    Given a distance, position_angle offset relative to a sky position,
    return the equivalent (l,m)

    :param dp: distance, position_angle offset
    :return: rectangular offset

    """
    #
    x = math.sin(dp[0]) * math.sin(dp[1])
    y = math.sin(dp[0]) * math.cos(dp[1])
    return x, y


def lm_to_dp(lm):
    """
    Given an (l,m) rectangular offset relative to a sky position,
    return the equivalent distance,position angle

    :param lm: rectangular offset
    :return:  distance,position angle
    """
    p = math.atan2(lm[0], lm[1])
    d = math.asin(math.sqrt(lm[0] * lm[0] + lm[1] * lm[1]))
    return d, p


def lm_to_true(lm):
    """
    Given an (l,m) rectangular offset relative to a sky position,
    return the equivalent true angle offsets in the antenna frame
    as used by casa ms

    :param lm: rectangular offset (orthographic projection)
    :return:  true angle offsets a,b (radians)
    """

    # note, adding zero to avoid -0.0000 representation of zero breaking tests
    a = math.asin(lm[0])
    b = math.atan2(lm[1], math.cos(math.asin(math.sqrt(lm[0]**2 + lm[1]**2))))

    return [a+0, b+0]


def rd_xyz(ra, dec):
    """TBD"""
    v = [math.cos(ra) * math.cos(dec), math.sin(ra) * math.cos(dec),
         math.sin(dec)]
    return v


def from_parset(parset, pitch, angle=0.0, lm_boresight=None):
    """

    :param parset: a ParameterSet defining the footprint
    :param pitch: factor by which to scale footprint by in radians
    :param angle: a single optional value (radians) by which to rotate
                  the footprint by
    :param lm_boresight: offset from boresight in radians
    :return: :obj:`Footprint`

    """
    if np.isclose(pitch, 0.) and (parset.n_beams > 1):
        raise ValueError('Footprint beam pitch may not be zero for n_beams={}.'.format(parset.n_beams))

    if parset.type == 'rect':
        # noinspection PyUnresolvedReferences
        pitch_scale = np.sin(pitch)
        offs = [[o[0]*pitch_scale, o[1]*pitch_scale] for o in parset.offsets]
        ioffs = [[i[0]*pitch_scale, i[1]*pitch_scale] for i in parset.interleaves]
        toffs = [[to[0]*pitch_scale, to[1]*pitch_scale] for to in parset.tile_offsets]
    elif parset.type == 'polar':
        pitch_scale = pitch
        offs = [[o[0]*pitch_scale, o[1]] for o in parset.offsets]
        ioffs = [[i[0]*pitch_scale, i[1]] for i in parset.interleaves]
        toffs = [[to[0]*pitch_scale, to[1]] for to in parset.tile_offsets]
    else:
        raise IndexError('Unsupported footprint offset type: {}'.format(parset.type))

    #logger.debug('Initialising {} footprint of type {} from parset'.format(parset.name, parset.type))
    pfp = Footprint(offs, ioffs, parset.type, angle, lm_boresight, tile_offsets=toffs)
    pfp.name = parset.name
    pfp.pitch_scale = pitch_scale
    return pfp


# DMcC: Amended 2015-03-11 to provide the PA correction required to maintain
# the PA in offset beam operations.
# DMcC: Amended 2015-05-05 to add standard_pitch function
class Footprint(object):
    def __init__(self, offsets, inter_offsets, offset_type="polar", angle=0.0,
                 lm_boresight=None, tile_offsets=None):
        """Expects all angular quantities in radians

        :param offsets: n_beam pairs of beam offsets of type offset_type
        :param inter_offsets: Pairs of type offset_type giving interleaving
                              offsets
        :param offset_type: "polar" or "rectangular", or "absolute";
                            if absolute, offsets is a list of positions
        :param angle: a single value (radians) by which to rotate the footprint.

        ASKAP forms up to 36 dual-polarision beams.  Their arrangement on the sky, called the
        "footprint", is arbitrary provided all lie with a six-degree circle.  This program
        allows the defintion of footprints from standard patterns or from customised lists of
        each beam position.

        Beams are defined on the sky relative to the "boresight" - the direction of the
        optical axis of the antenna.  The beam offset must be expressed in a way that is
        independent of the actual boresight direction. Two methods are used:

          1. Polar: distance, position angle (D,PA) where the distance is the great-circle
             distance between the boresight and the offset beam, and the position angle is
             measured as usual from celestial north, counterclockwise looking at the sky,
             that is "through east".
          2. Rectangular: (l,m) are the offset coordintes after an orthographic projection
             referenced to the Boresight.

        For each beam Footprint gives two kinds of output:
          1. A pointing direction for the boresight direction (command to the antenna) to
             place the specified reference direction in the beam.
          2. A position on the sky, given a reference position, that marks the beam's peak
             response.

        These two differ by 180 degrees of position angle.  The natural way to refer to the
        beam's offset is the second of these.  When calculating pointing required for forming
        or calibrating a beam, the offset must be shifted by 180 deg in PA (or by negating (l,m)).

        When operating the 3-axis antenna in equatorial mode (pa_fixed), simply shifting the offset
        by 180 degrees to place the offset beam on the desired point (cal source, or beam-forming
        source) is not correct. This is because the footprint will rotate so as to keep its n-s axis
        aligned north; that is the boresight beam keeps PA = 0.  To correct this, we need to apply a
        small feed rotation at the offset positions so that the offset beam maintains its position
        angle.
        """
        self.refpos = [0.0, 0.0]
        self.n_beams = len(offsets)
        """number of beams"""
        self.pa_corr = []
        self.pa = angle
        self.offsets = offsets
        """dimensionless beam offsets"""
        self.interleaved_offsets = inter_offsets
        """dimensionless tile offsets"""
        self.tile_offsets = tile_offsets

        """pitch in radians applied to parset footprint definition to get this object"""
        self.pitch_scale = None

        self.offset_type = offset_type
        """offset type e.g. rect, polar"""
        self.name = "custom"
        """The given name"""
        if "po" in offset_type:
            a = np.array(offsets)
            b = np.array(inter_offsets)
        elif "rect" in offset_type:
            a = np.array([lm_to_dp(o) for o in offsets])
            b = np.array([lm_to_dp(o) for o in inter_offsets])
        elif "abs" in offset_type:
            self.refpos = offsets[0]
            b0 = Skypos(*offsets[0])
            a = np.array([b0.d_pa(Skypos(p[0], p[1])) for p in offsets])
            b = np.array([b0.d_pa(Skypos(p[0], p[1])) for p in inter_offsets])
        else:
            raise KeyError("Unknown offset_type '{}'".format(offset_type))
        if len(b) == 0:
            b = np.array([[0.0, 0.0]])
        self.offsetsPolar = np.array([a.T[0], a.T[1] + angle]).T
        self.offsetsRect = np.array([dp_to_lm(o) for o in self.offsetsPolar])
        self.interOffsPolar = np.array([b.T[0], b.T[1] + angle]).T
        self.interOffsRect = np.array(
            [dp_to_lm(o) for o in self.interOffsPolar])
        if lm_boresight is not None:
            self.offsetsRect += lm_boresight
            self.offsetsPolar = np.array(
                [lm_to_dp(o) for o in self.offsetsRect])
            self.interOffsPolar = np.array(
                [lm_to_dp(o) for o in self.interOffsRect])
            self.interOffsRect += lm_boresight

        # now compute the reversed set of offsets
        a = self.offsetsPolar
        self.offsetsPolarReverse = np.array([a.T[0], a.T[1] + math.pi]).T
        self.positions = []
        self.positionsReverse = []
        self.n_beams = len(self.offsetsPolar)
        self.interLeaves = []
        self.il_pa_corr = []
        self.trueAngleOffsets = [lm_to_true(lm) for lm in self.offsetsRect]

    def set_refpos(self, refpos):
        """ Expects `refpos` as pair of radians"""
        t = Skypos(*refpos)
        self.refpos = t
        self.positions = [t.offset(dpa) for dpa in self.offsetsPolar]
        self.positionsReverse = [t.offset(dpa) for dpa in
                                 self.offsetsPolarReverse]
        self.pa_corr = [
            pr.d_pa(t)[1] - t.d_pa(p)[1] if abs(t.d_pa(p)[0]) > 1.0e-5 else 0.0
            for pr, p in zip(self.positionsReverse, self.positions)]

        self.interLeaves = [self.refpos.offset(dpa) for dpa in
                            self.interOffsPolar]
        iop_r = [[a[0], a[1] + math.pi] for a in self.interOffsPolar]
        il_r = [self.refpos.offset(dpa) for dpa in iop_r]
        pa_corr = [
            pr.d_pa(t)[1] - t.d_pa(p)[1] if abs(t.d_pa(p)[0]) > 1.0e-5 else 0.0
            for pr, p in zip(self.interLeaves, il_r)]
        self.il_pa_corr = pa_corr

    def get_true_angle_offsets(self):
        """
        Get true angle offsets for MS metadata
        :return: list of offset pairs [[delta_lon1, delta_lat1], [delta_lon2, delta_lat2], ...], one pair per beam
        :rtype: nested list of [float, float]
        """
        return np.array(self.trueAngleOffsets)

    def get_positions(self, reverse=False):
        """
        Get computed beam positions.

        :param reverse: return the reverse list
        :return: list of n_beam direction on the sky pairs
        """
        if reverse:
            return np.array(self.positionsReverse)
        else:
            return np.array(self.positions)

    def get_interleaves(self):
        return np.array(self.interLeaves)

    def get_tile_offsets(self):
        return np.array(self.tile_offsets)

    def get_interleave_pa(self):
        return np.array(self.il_pa_corr)

    def get_interleaved_footprint(self, ileave_num):
        offset_center = self.get_interleaves()[ileave_num]
        offset_pa = self.get_interleave_pa()[ileave_num]
        lm_boresight = (offset_center.ra, offset_center.dec)
        fpb = Footprint(self.offsetsRect, self.interOffsRect,
                        'rectangular', angle=offset_pa)
        fpb.set_refpos(lm_boresight)
        return fpb

    def get_pa_corr(self):
        ret = (3 * math.pi + np.array(self.pa_corr)) % (2 * math.pi) - math.pi
        return ret

    def to_parset(self):
        """
        Export the footprint as a ParameterSet.  Pitch scaling is undone before export.
        """
        out = ParameterSet()
        out.name = self.name
        out.n_beams = len(self.offsets)
        out.url = ''
        out.type = self.offset_type
        if self.offset_type == 'rect':
            out.offsets = [[o[0]/self.pitch_scale, o[1]/self.pitch_scale] for o in self.offsets]
            out.interleaves = [[o[0]/self.pitch_scale, o[1]/self.pitch_scale] for o in self.interleaved_offsets]
            out.tile_offsets = [[o[0]/self.pitch_scale, o[1]/self.pitch_scale] for o in self.tile_offsets]
        elif self.offset_type == 'polar':
            out.offsets = [[o[0]/self.pitch_scale, o[1]] for o in self.offsets]
            out.interleaves = [[o[0]/self.pitch_scale, o[1]] for o in self.interleaved_offsets]
            out.tile_offsets = [[o[0]/self.pitch_scale, o[1]] for o in self.tile_offsets]
        else:
            raise ValueError('Invalid footprint offset type: {}'.format(self.offset_type))

        return out

    @classmethod
    def standard_pitch(cls, band_freq, dish_dia=12.0, spacing_factor=0.75):
        """
        Return standard pitch based on lambda on d.
        Set factor for spacing relative to lambda/D (at band centre)
        Require that the response internal to the footprint does not fall below
        50% at the top of the band. For band 1, FWHM is 1.02*lambda/D
        (full illumination), HMR (radius) is 0.72deg. If the centroid to vertex
        distance of an equilateral triangle is x, the triangle has sides
        root(3).x
        Therefore we want a pitch in band 1 of
        1.25 deg = 0.75 * lambda(midband)/D

        :param band_freq:  frequency in MHz
        :param dish_dia: dish diameter in meters
        :param spacing_factor:
        :return: pitch in radians

        """
        approx_light_speed = 3e8  # m/s
        band_freq_hz = band_freq * 1e6
        lambda_on_d = approx_light_speed / band_freq_hz / dish_dia
        return spacing_factor * lambda_on_d
