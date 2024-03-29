from builtins import str
from builtins import range
from past.builtins import basestring
from builtins import object
# Copyright (c) 2009-2017 CSIRO
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
import os
import re
import warnings
import ast

__all__ = ["ParameterSet", "dict_to_parset", "parset_to_dict", "slice_parset",
           "sub_parset", "merge"]


def to_dict(parmap, typed=False):
    if isinstance(parmap, dict):
        return parmap
    elif isinstance(parmap, ParameterSet):
        return parmap.to_flat_dict(typed)
    else:
        raise TypeError("argument not a ParamterSet (or dict)")


def slice_parset(d, key):
    """Return a subset of the `dict` d  return matching key prefix only"""
    return dict((k, v) for k, v in list(d.items()) if k.startswith(key))


def sub_parset(d, key):
    return dict((k.replace(key + ".", ""), v)
                for k, v in list(d.items()) if k.startswith(key))


def parset_to_dict(st, raw=False):
    """Turn a parameterset string into a python `dict`"""
    d = {}
    for line in st.splitlines():
        pval, comm = extract(line)
        if pval:
            val = raw and pval[1] or decode(pval[1])
            d[pval[0]] = val
    return d


def dict_to_parset(d, sort=False):
    """Turn a python `dict` into a parameterset string.

    param d:    the dictionary object
    param sort: whether to sort by key

    """
    out = ""
    keys = list(d.keys())
    if sort:
        keys = sorted(keys)
    for k in keys:
        line = " = ".join((k, encode(d[k])))
        out = "".join((out, line, "\n"))
    return out


sentinel = object()


class ParameterSet(object):
    """
    The default constructor creates an empty ParameterSet instance.
    ParameterSet keys can be accessed as attributes::

        >>> p = ParameterSet('x.y.z', 1)
        >>> print(p.x.y.z)
        1

    or dictionary keys::

        >>> p = ParameterSet('x.y.z', 1)
        >>> print(p["x"]["y"]["z"])
        1
        >>> print(p["x.y.z"])
        1

    or using :meth:`get_value` which can be used with a default value if the
    key doesn't exists::

        >>> p = ParameterSet('x.y.z', 1)
        >>> print(p.get_value('x.y.z'))
        1
        >>> print(p.get_value('x.y.a', 10))
        10

    To test if a key exists us the `in` keyword. Note that the key must start
    at the root value::

        >>> p = ParameterSet('x.y.z', 1)
        >>> print('x' in p)
        True
        >>> print('y' in p)
        False

    :param args: if only one argument, it has to be a string
                 which is assumed to be the file name of
                 a ParameterSet file or a `dict` object.
                 If two arguments are provided this is
                 assumed to be a key and a value.
    :param kw:   key/value parameters

    Examples::

        >>> p0 = ParameterSet()
        >>> p1 = ParameterSet('x.y', 1)
        >>> x = ParameterSet('a', 3)
        >>> p2 = ParameterSet(x=1, y=2, z=x)
        >>> print(p2.x)
        1
        >>> print(p2["y"])
        2
        >>> print(p2.z.get_value('a'))
        3
        >>> p3 = ParameterSet('xyz.parset')
        >>> p1.set_value('x.a', 2)
        >>> p4 = ParameterSet({'x.y': 1})
        >>> p5 = ParameterSet(**{'x.y': 1})

    """

    def __init__(self, *args, **kw):
        object.__setattr__(self, "_keys", [])
        object.__setattr__(self, "_pdict", {})
        object.__setattr__(self, "_docdict", {})
        # from file
        if len(args) == 1:
            if isinstance(args[0], basestring):
                if args[0].find("=") > -1:
                    pfile = args[0].splitlines()
                elif os.path.exists(args[0]):
                    pfile_file = open(args[0], "r")

                    # split into lines here so that we can close pfile nicely
                    pfile = pfile_file.readlines()

                    pfile_file.close()
                else:
                    raise OSError("Given (single) argument is not filename.")
                i = 1
                pairs = []
                doc = ""
                for line in pfile:
                    pair, comment = extract(line, doc)
                    if pair and pair in pairs:
                        lineno = pairs.index(pair) + 1
                        msg = "Overwriting value for key '%s' first " \
                              "defined in line %d" % (pair[0], lineno)
                        warnings.warn(msg)
                    pairs.append(pair)
                    if pair:
                        try:
                            self.set_value(pair[0], pair[1], comment)
                            doc = ""
                        except ValueError as ex:
                            raise ValueError(
                                "In line %d of %s. %s" % (i, args[0],
                                                          ex.message))
                        except IndexError as ex:
                            raise IndexError(
                                "In line {} of {}: {} {}".format(i, args[0],
                                                                 line, ex.message))
                    else:
                        doc = comment
                    i += 1
            elif isinstance(args[0], dict):
                for k, v in list(args[0].items()):
                    self.set_value(k, v)
            else:
                raise ValueError(
                    "Given (single) argument is not a file name or dict")
        # from key, value
        elif len(args) == 2:
            self.set_value(*args)
        elif len(kw):
            for k, v in list(kw.items()):
                self.set_value(k, v)
        elif len(args) == 0 and len(kw) == 0:
            pass
        else:
            raise ValueError("Incorrect arguments to constructor.")

    # noinspection PyIncorrectDocstring
    def get_value(self, k, default=sentinel, _full=None):
        """Return the value from the ParameterSet using an optional default
        value if the key doesn't exist.

        :param k: the key to get the value for
        :param default: the default value to return if the key doesn't exist

        """
        fullkey = _full or k
        inkey = k
        k, tail = self._split(inkey)
        if k in self._keys:
            child = self._pdict[k]
            if isinstance(child, self.__class__):
                if tail is not None:
                    return child.get_value(tail, default, fullkey)
                else:
                    return decode(child)
            else:
                if tail is not None:
                    if default is sentinel:
                        raise KeyError("Key '%s' not found." % fullkey)
                    else:
                        return default
                return decode(child)
        else:
            if default is sentinel:
                raise KeyError("Key '%s' not found." % fullkey)
            else:
                return default

    def get_doc(self, k):
        """Get the documentation for the specified key `k`"""
        inkey = k
        k, tail = self._split(inkey)
        if k in self._keys:
            child = self._pdict[k]
            if isinstance(child, self.__class__):
                if tail is not None:
                    return child.get_doc(tail)
                else:
                    return self._docdict[k]
            else:
                return self._docdict[k]
        else:
            raise KeyError("Key '%s' not found." % inkey)

    def keys(self):
        """Return the name of the keys in this `ParameterSet`."""
        out = []
        for k in self._keys:
            child = self._pdict[k]
            if isinstance(child, self.__class__):
                for key in list(child.keys()):
                    out.append(".".join([k, key]))
            else:
                out.append(k)
        out.sort()
        return out

    def set_doc(self, k, doc):
        """Set the documenation for the given key
        
        :param string k:    the key
        :param string doc:  the (new) documentation

        """

        v = self.get_value(k)
        self.set_value(k, v, doc)

    def set_value(self, k, v, doc=""):
        """
        Add or replace key/value pair. This will recursively create keys if
        necessary when the key contains '.' notation. This is the only way to
        add keys of this form. To set non-nested attributes one can use
        attribute or item set notation, so that the following are equivalent::

            p = ParameterSet()
            p.x = 1
            p["x"] = 1
            p.set_value("x", 1)
            # to add nested keys use
            p.set_value('x.y', 1)
            # this fails
            p.x.y = 1
            # as does this
            p["x"]["y"]

        :param k: key
        :param v: value
        :param doc: documentation string

        """
        k, tail = self._split(k)
        if k in self._keys:
            child = self._pdict[k]
            if isinstance(child, self.__class__):
                child.set_value(tail, v, doc)
            else:
                if tail:
                    raise ValueError("Leaf node %s can't be extended" % k)
                else:
                    self._pdict[k] = v
                    if k in self._docdict:
                        if doc:
                            self._docdict[k] = doc
                    else:
                        self._docdict[k] = doc
        else:
            if not tail:
                self._pdict[k] = v
                if k in self._docdict:
                    if doc:
                        self._docdict[k] = doc
                else:
                    self._docdict[k] = doc
            else:
                child = self.__class__()
                child.set_value(tail, v, doc)
                self._pdict[k] = child
            self._keys.append(k)

    def __len__(self):
        return len(list(self.keys()))

    def __setitem__(self, k, v):
        self.set_value(k, v)

    def __setattr__(self, k, v):
        self.set_value(k, v)

    def __getitem__(self, k):
        return self.get_value(k)

    def __delitem__(self, k):
        inkey = k
        k, tail = self._split(inkey)
        if k in self._keys:
            child = self._pdict[k]
            if isinstance(child, self.__class__) and tail is not None:
                if tail is not None:
                    del child[tail]
                    return
            else:
                del self._pdict[k]
                if k in self._docdict:
                    del self._docdict[k]
                self._keys.remove(k)
        else:
            raise KeyError("Key '%s' not found." % inkey)

    def __getattr__(self, k):
        return self.get_value(k)

    @staticmethod
    def _split(k):
        keys = k.split(".")
        k = keys[0]
        tail = None
        if len(keys) > 1:
            tail = ".".join(keys[1:])
        return k, tail

    def __contains__(self, k):
        k, tail = self._split(k)
        if k in self._keys:
            child = self._pdict[k]
            if isinstance(child, self.__class__):
                if tail is not None:
                    return child.__contains__(tail)
                else:
                    return True
            else:
                return tail is None
        else:
            return False

    def to_dict(self):
        """
        Returns a python :class:`dict` representation of the `ParameterSet`,
        decoding all values using :func:`decode`
        """
        out = {}
        for k in self._keys:
            if isinstance(self._pdict[k], self.__class__):
                out[k] = self._pdict[k].to_dict()
            else:
                out[k] = decode(self._pdict[k])
        return out

    def to_flat_dict(self, typed=False):
        """
        Returns a python :class:`dict` representation of the `ParameterSet`,
        with flat keys and encoded (string) values when typed==False or python
        typed values when typed==False

        :param bool typed: choosed to have values python typed or ParameterSet
                           encoded strings.
        """
        if typed:
            return dict([(k, v) for k, v in list(self.items())])
        return dict([(k, encode(v)) for k, v in list(self.items())])

    def _get_strings(self):
        """
        Get a list of key=values strings as they appear in ParameterSet files.
        """
        return ["%s = %s" % (k, encode(v)) for k, v in list(self.items())]

    def __str__(self):
        return "\n".join(self._get_strings())

    def __repr__(self):
        return self.__str__()

    def to_file(self, filename):
        f = open(filename, 'w')
        f.write(str(self) + '\n')
        f.close()

    def __iter__(self):
        yield iter(list(self.keys()))

    def pop(self, key):
        out = self[key]
        del self[key]
        return out

    def items(self):
        return [(k, self.get_value(k)) for k in list(self.keys())]


def encode(value):
    """Encode a python value as ParameterSet string.

       Note that python hex values are `int`s, so the only way to encode
       handle hex values is to make them strings.
    """

    def single_str(v):
        if isinstance(v, bool):
            return v and 'true' or 'false'
        if isinstance(v, str):
            return v
        return v

    value = single_str(value)
    # deal with numpy arrays, by converting to lists
    if hasattr(value, 'tolist'):
        value = value.tolist()
    if isinstance(value, list) or isinstance(value, tuple):
        # noinspection PyUnusedLocal
        def to_str(val, islement=False):
            if isinstance(val, list) or isinstance(val, tuple):
                vals = []
                for j in val:
                    # quote text in lists which contains whitespace
                    if isinstance(j, basestring) and len(j.split()) > 1:
                        j = '"' + j + '"'
                    vals.append(to_str(j))
                return "[" + ", ".join(vals) + "]"
            else:
                if isinstance(val, bool):
                    return val and 'true' or 'false'
                return str(val)

        if len(value) > 2 and not isinstance(value[0], list):
            # [ n * val ]
            if value == [value[0]] * len(value):
                val = str(single_str(value[0]))
                return '[' + str(len(value)) + ' * ' + val + ']'
            # n..m
            elif all([isinstance(i, int) for i in value]) \
                    and value == list(range(value[0], value[-1] + 1)):
                return str(value[0]) + '..' + str(value[-1])
        return to_str(value)
    return str(value)


def _fromdotdot(match, value):
    rxisrange = re.compile(r"(\d+)\.{2}(\d+)")
    r = match.groups()
    r0 = int(r[0])
    padded = False
    if r[0].startswith("0") and len(r[0]) > 1:
        padded = True
    r1 = int(r[1])
    if r[1].startswith("0") and len(r[1]) > 1:
        padded = True
    sgn = (r0 < r1) and 1 or -1
    rng = list(range(r0, r1 + sgn, sgn))
    # just numerical range
    if match.span()[0] == 0 and match.span()[1] == len(value):
        return rng
    nwidth = 1
    if padded:
        nwidth = max(len(r[0]), len(r[1]))
    strg = rxisrange.sub("%%0%ii" % nwidth, value)
    return [strg % i for i in rng]


def decode(value):
    """
    This function takes text a string which is using ParameterSet syntax
    and is decoding it into a valid python value, e.g.:: python

        >>> p = ParameterSet('x.y', '[a,b]')
        >>> print(decode(p.x.y))
        ['a', 'b']

    Supported value encodings are:

    * ranges i..j or j..i with padding, e.g. abc00..09.txt
    * lists
    * numerical arrays (lists of lists with numerical values)
    * booleans true/false

    """
    if not isinstance(value, basestring):
        return value
    rxislist = re.compile(r"^\[(.*)\]$")
    rxbool = re.compile(r"([tT]rue|[fF]alse)")
    rxisrange = re.compile(r"(\d+)\.{2}(\d+)")
    rxisnum = re.compile(r"^([+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)$")
    rxishex = re.compile("^0x[0-9a-fA-F]+$")
    rxstr = re.compile('^[\'"](.*)[\'"]$')
    # lists/arrays
    match = rxislist.match(value)
    if match:
        # look for  '01..10' type pattern inside []
        svalue = match.groups()[0]
        if svalue == '':
            return []
        submatch = rxisrange.search(svalue)
        if submatch:
            return _fromdotdot(submatch, svalue)
        # check for [ n * <value> ] and expand
        # doesn't work for vectors elements
        if value.count(",") == 0:
            rxtimes = re.compile(r"\[.*?(\d+)\s*[*](.+)\]")
            mtch = rxtimes.match(value)
            if mtch:
                fact = int(mtch.groups()[0])
                val = mtch.groups()[1].strip()
                if rxbool.match(val):
                    val = ast.literal_eval(val.title())
                elif rxishex.match(val):
                    val = ast.literal_eval(val)
                elif rxisnum.match(val):
                    val = ast.literal_eval(val)
                return fact * [val]

        if value.count("[") > 1:
            inner = value[1:-1]
            rxitem = re.compile(r"(\[.+?\])")
            out = []
            for item in rxitem.findall(inner):
                out.append(decode(item))
            return out
        out = []
        items = match.groups()[0].split(",")
        for i in items:
            i = i.strip()
            if rxisnum.match(i):
                i = ast.literal_eval(i)
            elif rxbool.match(i):
                i = ast.literal_eval(i.title())
            elif rxstr.match(i):
                i = rxstr.match(i).groups()[0]
            out.append(i)
        return out
    # look for  '01..10' type pattern
    match = rxisrange.search(value)
    if match:
        return _fromdotdot(match, value)
    # hex
    if rxishex.match(value):
        return ast.literal_eval(value)
    # int/float
    if rxisnum.match(value):
        return ast.literal_eval(value)
    # true/false
    if rxbool.match(value):
        return ast.literal_eval(value.title())
    return value


def extract(line, comment=""):
    """
    Return a key/value pair from a string. This will most likely be a line in a
    ParameterSet file. It also returns (optionally) a documentation string from
    lines beginning with '##'.
    """
    line = line.lstrip()
    line = line.strip("\n")
    if len(line) == 0:
        return None, comment
    elif line.startswith("#"):
        if line.startswith("##"):
            if len(comment):
                line = line.strip("#")
                if len(line) and line.startswith(" "):
                    line = line[1:]
                return None, "\n".join((comment, line))
            else:
                line = line.strip("#")
                if len(line) and line.startswith(" "):
                    line = line[1:]
                return None, line
        return None, comment
    else:
        line = line.strip()
        kv = line.split("=", 1)
        return [i.strip() for i in kv], comment


# noinspection PyUnusedLocal
def merge(*parsets, **kw):
    """Merge n :class:`ParameterSet` into one. Existing keys will be 
    overwritten, so that the last given argument wins.

    :param parsets: positional argument list of :class:`ParameterSet`.
    
    """
    if len(parsets) < 2:
        raise ValueError("Need to or more ParameterSets")
    out = ParameterSet()
    for p in parsets:
        for k, v in list(p.items()):
            doc = p.get_doc(k)
            out.set_value(k, v, doc)

    return out
