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

obuild=38
force_ibuild=None
min_qual=20
min_mapqual=20
min_reads=1
read_max=100

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--obuild', help='Output build')
parser.add_argument('-i', '--ibuild', help='Force input build')
parser.add_argument('-b', '--bqual', help='Min base quality')
parser.add_argument('-q', '--mqual', help='Min map quality')
parser.add_argument('-m', '--minread', help='Min reads to qualify')
parser.add_argument('-x', '--maxread', help='Max reads to analyze per location')
parser.add_argument('--no-mt', action="store_true", help='Skip MT')
parser.add_argument('--no-y', action="store_true", help='Skip Y')
parser.add_argument('--no-ystr', action="store_true", help='Skip Y-STR')
parser.add_argument('--auto', action="store_true", help='Convert autosomal')
parser.add_argument('--ybrowse', action="store_true", help='Use full ybrowse database positions')
parser.add_argument('file', nargs='+')
args = parser.parse_args()

if args.obuild:
    obuild = int(args.obuild)
if args.mqual:
    min_mapqual = int(args.mqual)
if args.bqual:
    min_qual = int(args.bqual)
if args.minread:
    min_reads = int(args.minread)
if args.maxread:
    read_max = int(args.maxread)
if args.no_mt:
    bamload.convert_mt=0
if args.no_y:
    bamload.convert_y=0
if args.no_ystr:
    bamload.convert_ystr=0
if args.auto:
    bamload.convert_snpauto=1
if args.ybrowse:
    bamload.load_ybrowse=1

bamload.min_qual=min_qual
bamload.min_mapqual=min_mapqual
bamload.min_reads=min_reads
bamload.read_max=read_max

for samfname in args.file:
    if args.ibuild:
        in_build = int(args.ibuild)
    else:
        in_build = bamload.get_build(samfname)
    bamload.setup_conv(in_build)
    snpset = bamload.full_convert(samfname)
    snpload.save(samfname+'.snp.txt', snpset, obuild)
