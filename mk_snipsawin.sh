#!/bin/sh

set -e

#
# Win  install:
# https://www.kannon.link/free/2019/06/24/distribution-using-a-python-embeddable-zip-file/
#  - python.exe get-pip.py
#  - python37._pth
#  - python.exe -m pip install bs4
#  - tkinter
#

#python embed and pre-installed tk/other libs as a corresponding zip is needed in the folder
SNIPSAWIN_DIR=snipsawin
PYTHON_EMBED=python-3.8.5-embed-amd64-preinstalled

echo Preparing Snipsa Windows package...

rm -rf $SNIPSAWIN_DIR/ || true
rm snipsawin.zip || true
mkdir $SNIPSAWIN_DIR

#python
echo Installing python environment ...
mkdir $SNIPSAWIN_DIR/$PYTHON_EMBED
unzip -q $PYTHON_EMBED.zip -d $SNIPSAWIN_DIR

echo python38.zip > temp._pth
echo . >> temp._pth
echo .\\Lib >> temp._pth
echo .\\DLLs >> temp._pth
echo ..\\snipsa >> temp._pth
echo import site >> temp._pth
mv temp._pth $SNIPSAWIN_DIR/$PYTHON_EMBED/python38._pth

#snipsa
echo Installing snipsa ...
mkdir $SNIPSAWIN_DIR/snipsa
cp \
	haplomt.py haplomt_find.py haplomt_db_import.py \
	haploy.py haploy_find.py haploy_db_import.py \
	haploy_anno_import.py \
	snpload.py \
	bamload.py \
	snipsa-gui.py \
	$SNIPSAWIN_DIR/snipsa/

#DB licenses included in implementation
#TODO: FTDNA
cp haplomt_map.txt haploy_map2j.txt haploy_annodb_example.txt haploy_annodb_yfull.txt haploy_annodb_ancientdna.txt $SNIPSAWIN_DIR/snipsa/

# Ybrowse (CC BY-NC-SA 3.0)
# https://isogg.org/wiki/ISOGG_Wiki:Terms_of_Use
#cp str_hg19.gff3 str_hg38.gff3 snps_hg38.gff3 $SNIPSAWIN_DIR/snipsa/
cp str_hg19.gff3 str_hg38.gff3 $SNIPSAWIN_DIR/snipsa/

#http://May2021.archive.ensembl.org/info/about/legal/disclaimer.html free use license
cp -a crossmap $SNIPSAWIN_DIR/snipsa/

#starter
cp snipsa.bat $SNIPSAWIN_DIR/

echo Zipping Snipsa Windows package...
zip -q -r snipsawin.zip $SNIPSAWIN_DIR/
rm -rf $SNIPSAWIN_DIR/
echo Done!
