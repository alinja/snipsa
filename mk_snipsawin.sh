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
PYTHON_EMBED=python-3.9.2-embed-amd64

echo Preparing Snipsa Windows package...

rm -rf $SNIPSAWIN_DIR/
rm snipsawin.zip
mkdir $SNIPSAWIN_DIR

#python
echo Installing python environment ...
mkdir $SNIPSAWIN_DIR/$PYTHON_EMBED
unzip -q $PYTHON_EMBED.zip -d $SNIPSAWIN_DIR/$PYTHON_EMBED
unzip -q python-3.9.2-libs.zip -d $SNIPSAWIN_DIR/$PYTHON_EMBED

echo python39.zip > temp._pth
echo . >> temp._pth
echo .\\Lib >> temp._pth
echo .\\DLLs >> temp._pth
echo ..\\snipsa >> temp._pth
echo import site >> temp._pth
mv temp._pth $SNIPSAWIN_DIR/$PYTHON_EMBED/python39._pth

#snipsa
echo Installing snipsa ...
mkdir $SNIPSAWIN_DIR/snipsa
cp \
	haplomt.py haplomt_find.py haplomt_db_import.py \
	haploy.py haploy_find.py haploy_db_import.py \
	haplomt_map.txt haploy_map2.txt \
	snpload.py \
	snipsa-gui.py \
	$SNIPSAWIN_DIR/snipsa/

#starter
cp snipsa.bat $SNIPSAWIN_DIR/

echo Zipping Snipsa Windows package...
zip -q -r snipsawin.zip $SNIPSAWIN_DIR/
rm -rf $SNIPSAWIN_DIR/
echo Done!
