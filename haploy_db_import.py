#!/usr/bin/python3
import sys
import os
import haploy

do_yfull_snps=0
include_ftdna_tree=1
#use_ftdna_tree=1

print("Loading SNP DB from ISOGG csv...")
haploy.load_snp()

print("Loading SNP DB from YBrowse...")
haploy.load_ybrowse_snp()
haploy.save_alldbj()

if do_yfull_snps:
    print("Loading SNP DB from YFull...")
    haploy.load_yfull_snp(234)

if include_ftdna_tree:
    print("Loading SNP DB from FTDNA...")
    haploy.import_ftdna_tree()

print("Running conversion to support build 36 matching...")
haploy.convert_build38to36()

print("Running conversion to support build 37 matching...")
haploy.convert_build38to37()

print("Importing tree DB from YFull...")
haploy.import_yfull_tree()
#haploy.show_db2()

print("Writing local mutation DB (legacy)...")
haploy.save_db()

print("Writing local mutation DB (FTDNA tree json)...")
haploy.save_db3j()

print("Writing local mutation DB (Yfull tree json)...")
haploy.save_db2j()

print("Database import done!")
