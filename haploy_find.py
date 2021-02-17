#!/usr/bin/python3
import sys
import os
import snpload
import haploy


if len(sys.argv) < 2:
    print(sys.argv[0]+" <filename>")
    print(sys.argv[0]+" <filter> <filenames>..")
elif len(sys.argv) < 3:
    print("Loading DB...")
    haploy.load_db()
    print("DB loaded!")
    print("Loading chr data...")
    fname=sys.argv[1]
    snpset, meta = snpload.load(fname, ['Y'])
    
    if 'Y' not in snpset:
        print("No Y data found")
        sys.exit()
        
    print("%s: Total SNPs: %d"%(fname,meta['total']))

    found_mutations=haploy.yfind(snpset)
    for m in found_mutations:
        print("%-24s %-50s b38:%s %s"%(m['g'], m['mall'], m['b38'], m['rs']))
else:
    print("Loading DB...")
    haploy.load_db()
    print("DB loaded!")
    lookfor = sys.argv[1].split(',')
    for fname in sys.argv[2:]:
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
                
       
