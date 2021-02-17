#!/usr/bin/python3
import sys
import os
import snpload
import haplomt

#configure how many candidates are shown:
n_single = 1
n_multi = 5

if len(sys.argv) < 2:
    print(sys.argv[0]+" <filename>")
    print(sys.argv[0]+" <filter> <filenames>..")
elif len(sys.argv) < 3:
    print("Loading DB...")
    haplomt.load_db()
    print("DB loaded!")
    
    print("Loading chr data...")
    fname=sys.argv[1]
    snpset, meta = snpload.load(fname, ['MT'])
    
    if 'MT' not in snpset:
        print("No MT data found")
        sys.exit()
        
    print("%s: Total SNPs: %d"%(fname,meta['total']))

    best_trees = haplomt.mtfind(snpset, n_single)
    
    for bt in best_trees:
        haplomt.print_uptree(snpset, bt['ut'])
        leaf_mut = bt['ut'][len(bt['ut'])-1]
        print("Result (%-8s %5.1f%% -%d +%d): %-8s"%(leaf_mut['raw'], bt['score'], bt['neg'], len(bt['extras']), leaf_mut['g']))
        haplomt.print_extras(snpset, bt)

else:
    print("Loading DB...")
    haplomt.load_db()
    print("DB loaded!")
    
    lookfor = sys.argv[1].split(',')
    for fname in sys.argv[2:]:
        #try:
        snpset, meta = snpload.load(fname, ['MT'])
        
        if 'MT' not in snpset:
            print('%s: no MT data'%fname)
            continue

        best_trees = haplomt.mtfind(snpset, n_multi)
        
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
