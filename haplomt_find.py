#!/usr/bin/python3
import sys
import os
import snpload
import haplomt
import argparse

#configure how many candidates are shown:
n_single = 1
n_multi = 5
force=''
filt=''
all=False

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--single', help='Analyse a path for single group')
parser.add_argument('-a', '--all', action='store_true', help='Show listing of all found mutations')
parser.add_argument('-n', '--num', help='Show num best matches')
parser.add_argument('file', nargs='+')

args = parser.parse_args()
if args.single:
    force = args.single
    filt = '='+args.single
if args.num:
    n_single = int(args.num)
    n_multi = int(args.num)
if args.all:
    all=True

if len(args.file) < 1:
    print(sys.argv[0]+" <filename>")
    print(sys.argv[0]+" <filter> <filenames>..")
elif len(args.file) < 2:
    print("Loading DB...")
    haplomt.load_db()
    print("DB loaded!")
    
    print("Loading chr data...")
    rep = haplomt.report(args.file[0], n_single, do_all=all, filt=filt, force=force)
    print(rep)
    
else:
    print("Loading DB...")
    haplomt.load_db()
    print("DB loaded!")
    
    lookfor = args.file[0].split(',')
    for fname in args.file[1:]:
        #try:
        snpset, meta = snpload.load(fname, ['MT'])
        
        if 'MT' not in snpset:
            print('%s: no MT data'%fname)
            continue

        best_trees = haplomt.mtfind(snpset, n_multi, filt, force)
        
        found=0
        for bt in best_trees:
            for cmp in lookfor:
                leaf_mut = bt['ut'][len(bt['ut'])-1]
                if leaf_mut['g'].startswith(cmp) or cmp == '0':
                    found=1
                if leaf_mut['raw'].startswith(cmp):
                    found=1
        if found:
            print('%s: match found'%fname)
            for bt in best_trees:
                leaf_mut = bt['ut'][len(bt['ut'])-1]
                print("Result (%-8s %5.3f%% -%d +%d): %-8s"%(leaf_mut['raw'], bt['score'], bt['neg'], len(bt['extras']), leaf_mut['g']))
        else:
            print('%s: no match'%fname)
