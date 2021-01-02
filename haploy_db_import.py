#!/usr/bin/python3
import sys
import os
import haploy

print("Loading mutation DB from csv...")
haploy.load_snp()

print("Running conversion to support build 36 matching...")
haploy.convert_build38to36()

print("Writing local mutation DB...")
haploy.save_db()

print("Database import done!")
