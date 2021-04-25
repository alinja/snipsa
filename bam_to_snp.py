#!/usr/bin/python3
import os
import re
import sys
import csv
import urllib.request
import zipfile
import gzip
from pyliftover import LiftOver

import pysam
import argparse
import snpload
import bamload

#TODO: map quality
#TODO: mt str option

parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='+')
args = parser.parse_args()

samfname=args.file[0]
bamload.get_build(samfname)
snpset = bamload.full_convert(samfname)
snpload.save(samfname+'.snp.txt', snpset, 38)
