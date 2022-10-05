#!/bin/bash

# Create a set of zip files with an absorption result for a field.
# Usage: prepare_zip.sh sbid [folder] [prefix]
# The default folder is sb<sbid> in the current folder (e.g. ./sb14211) and the default prefix is sb<sbid>- (e.g. sb14211- )
# The zip file produced are:
# * ..abs_spectra.zip - Absorption spectra VOTable files for all sources.
# * ..em_spectra.zip - Emission spectra VOTable files for all sources.
# * ..plots.zip - Web pages and absoprtion spectra plots
# * ..diagnostics.zip - Diagnositc plots including moment 0 maps of emission and absorption data
# * ..ds9regions.zip - DS9 region files for detections and non detections.
# * ..catalogues.zip - The spectra and absorption detections catalogues in VOTable format.

if [ $# -lt 1 ]; then
    echo "Usage: ${0} sbid [folder] [prefix]"
    exit 1
fi

folder="sb$1"
if [ $# -gt 1 ]; then
    folder="$2"
fi
if [ ! -d "$folder" ]; then
    echo "Folder $folder does not exist."
    exit 1
fi
prefix="sb$1-"
if [ $# -gt 2 ]; then
    prefix="$3"
fi
echo "Producing zips for $1 in $folder with prefix $prefix"

# Prepare zips of spectra and plots
cd $folder
zipname="${prefix}abs_spectra.zip"
echo "Building $folder/${zipname}"
rm ${zipname}
zip -9 -q ${zipname} spectra/*_spec.vot

zipname="${prefix}em_spectra.zip"
echo "Building $folder/${zipname}"
rm ${zipname}
zip -9 -q ${zipname} spectra/*_emission.vot

zipname="${prefix}plots.zip"
echo "Building $folder/${zipname}"
rm ${zipname}
zip -9 -q ${zipname} figures/*_combined.png figures/*_loc*.png spectra/*_spec.png *.html 

zipname="${prefix}diagnostics.zip"
echo "Building $folder/${zipname}"
rm ${zipname}
zip -9 -q ${zipname} figures/*_mom0.png figures/*_data.png figures/*_weights.png figures/*_image_wcs.png

zipname="${prefix}ds9regions.zip"
echo "Building $folder/${zipname}"
rm ${zipname}
zip -9 -q ${zipname} *.reg

zipname="${prefix}catalogues.zip"
echo "Building $folder/${zipname}"
rm ${zipname}
zip -9 -q ${zipname} *.vot
