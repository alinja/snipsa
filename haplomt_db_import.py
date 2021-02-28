#!/usr/bin/python3
import sys
import os
import haplomt

phylotree=0

if phylotree:
    print('Loading mutation DB from "mtDNA tree Build 17.htm"')
    haplomt.import_snp()
else:
    print("Importing mutation DB from YFull...")
    haplomt.import_yfull_snp()

print("Writing local mutation DB...")
haplomt.save_db()

print("Database import done!")
