#!/usr/bin/python3
import sys
import os
import snpload
import haploy
import argparse

#configure how many candidates are shown:
n_single = 1
n_multi = 5
force=''
filt=''
all=False
min_match_level=0
min_tree_load_level=0
new_yfind=1

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--single', help='Analyse a path for single group')
parser.add_argument('-a', '--all', action='store_true', help='Show listing of all found mutations')
parser.add_argument('-n', '--num', help='Show num best matches')
parser.add_argument('-q', '--quick', help='Quick mode')
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
if args.quick:
    min_match_level=int(args.quick)
    min_tree_load_level=int(args.quick)

if len(args.file) < 2:
    
    if new_yfind:
        print("Loading DB2...")
        haploy.load_db2j(min_tree_load_level=min_tree_load_level)
        print("DB loaded!")
        rep = haploy.report(args.file[0], n_single, do_all=all, filt=filt, force=force, min_match_level=min_match_level)
        print(rep)
    else:
        # keep old one available
        print("Loading DB...")
        haploy.load_db()
        print("DB loaded!")
        print("Loading chr data...")
        fname=args.file[0]
        snpset, meta = snpload.load(fname, ['Y'])

        if 'Y' not in snpset:
            print("No Y data found")
            sys.exit()

        print("%s: Total SNPs: %d"%(fname,meta['total']))

        found_mutations=haploy.yfind(snpset)
        for m in found_mutations:
            print("%-24s %-50s b38:%s %s"%(m['g'], m['mall'], m['b38'], m['rs']))
else:
    lookfor = args.file[0].split(',')
    if new_yfind:
        print("Loading DB...")
        haploy.load_db2j(min_tree_load_level=min_tree_load_level)
        print("DB loaded!")
        for fname in args.file[1:]:
            snpset, meta = snpload.load(fname, ['Y'])

            if 'Y' not in snpset:
                print('%s: no Y data'%fname)
                continue

            #print('Build: %d'%meta['build'])
            b3x='b36'
            if meta['build']==37:
                b3x='b37'
            if meta['build']==38:
                b3x='b38'
            best_trees = haploy.yfind2(snpset, n_multi, filt, force, b3x, min_match_level)

            found=0
            for bt in best_trees:
                for cmp in lookfor:
                    leaf_mut = bt['ut'][len(bt['ut'])-1]
                    path = haploy.path_str(bt['ut'], 15)
                    if leaf_mut['g'].startswith(cmp) or cmp == '0':
                        found=1
                    if cmp in leaf_mut['raw']: ##todo
                        found=1
                    if cmp in path:
                        found=1
            if found:
                print('%s: match found'%fname)
                for bt in best_trees:
                    leaf_mut = bt['ut'][len(bt['ut'])-1]
                    #print("Result (%-8s %5.1f%% -%d +%d): %-8s (ISOGG: %s)"%(leaf_mut['raw'], bt['score'], bt['neg'], len(bt['extras']), haploy.path_str(bt['ut'], 15), leaf_mut['isog']))
                    print("Result (%.1f%% %d -%d +%d): %-8s"%(bt['score'], bt['tot'], bt['neg'], len(bt['extras']), leaf_mut['g']))
                    print("  %s"%(haploy.path_str(bt['ut'], 20)))
                    #haploy.print_extras(snpset, bt, True, b3x)
            else:
                print('%s: no match'%fname)
    else:
        print("Loading DB...")
        haploy.load_db()
        print("DB loaded!")
        for fname in args.file[1:]:
            snpset, meta = snpload.load(fname, ['Y'])
            if 'Y' not in snpset:
                print('%s: no Y data'%fname)
                continue
            found_mutations = haploy.yfind(snpset)
            found_mutations.reverse()
            found=0
            for m in found_mutations:
                for cmp in lookfor:
                    if m['g'].startswith(cmp) or cmp == '0':
                        found=1
                    if cmp in m['mall']:
                        found=1
            if found:
                print('%s: match found'%fname)
                for m in found_mutations[0:9]:
                    print("%-24s %-50s b38:%s %s"%(m['g'], m['mall'], m['b38'], m['rs']))
            else:
                print('%s: no match'%fname)
