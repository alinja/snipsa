#!/usr/bin/python3
import sys
import os
import haploy

do_yfull_snps=1

print("Loading SNP DB from ISOGG csv...")
haploy.load_snp()

if do_yfull_snps:
    print("Loading SNP DB from YFull...")
    haploy.load_yfull_snp(234)

print("Running conversion to support build 36 matching...")
haploy.convert_build38to36()

print("Importing tree DB from YFull...")
haploy.import_yfull_tree()
#haploy.show_db2()

print("Writing local mutation DB (legacy)...")
haploy.save_db()

print("Writing local mutation DB (tree)...")
#haploy.save_yfull_db()
haploy.save_db2()

print("Database import done!")
