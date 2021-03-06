#!/usr/bin/python3
import sys
import os
import haploy

print("Loading mutation DB from ISOGG csv...")
haploy.load_snp()

print("Running conversion to support build 36 matching...")
haploy.convert_build38to36()

print("Importing mutation DB from YFull...")
haploy.import_yfull_snp()
#haploy.show_db2()

print("Writing local mutation DB...")
haploy.save_db()

print("Writing local mutation DB2...")
haploy.save_db2()

print("Database import done!")
