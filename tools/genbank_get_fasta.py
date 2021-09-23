#!/usr/bin/python3
import urllib.request
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='+')
args = parser.parse_args()

code = args.file[0]
def get_genbank_page(code):
    url = 'http://www.ncbi.nlm.nih.gov/nuccore/' + code + '?report=fasta&log$=seqview&format=text'
    #print(url)
    htm = urllib.request.urlopen(url).read()
    return htm.decode('UTF-8')

def get_fasta(gid):
    url = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?id=' + gid + '&db=nuccore&report=fasta&retmode=text&withmarkup=on&tool=portal&log$=seqview'
    #print(url)
    htm = urllib.request.urlopen(url).read()
    return htm.decode('UTF-8')

genbank_page = get_genbank_page(code)
#print(genbank_page)

match = re.search(r'class="seq gbff" val="(\d+)"', genbank_page)
gid = match.group(1)
print('ID: ', gid)

fasta = get_fasta(gid)
print(fasta)

with open(code+'.fasta', "w") as f:
    f.write(fasta)