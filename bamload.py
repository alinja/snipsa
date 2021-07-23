#!/usr/bin/python3
import os
import re
import sys
import csv
import urllib.request
import zipfile
import gzip
import json
from pyliftover import LiftOver
import haploy

import pysam
import argparse
#import snpload

in_build=37
b3x='b37'
str_db_file='str_hg19.gff3'
contig='Y'
contigmt='MT'
pos_triplet_fn = None
lo_37to38 = None
lo_38to37 = None

load_ybrowse=0

snp_by_b37 = {}

def load_ysnp_dbj(min_tree_load_level=0):
    with open('haploy_map2j.txt', 'r') as f:
        jdata = json.load(f)
        for snp in jdata['muts']:
            if b3x in snp:
                snp_by_b37[snp[b3x]] = snp

def load_ysnp_ybrowse_db():
    haploy.load_ybrowse_snp()
    global snp_by_b37
    snp_by_b37 = haploy.haplo_ybrowse_muts_by_b38
    if b3x != 'b38':
        #TODO: build conversion
        raise LookupError


#Genotype calling from pileup column
min_qual=20
min_mapqual=20
min_reads=1
read_max=100
heterozyg_ratio=0.3
ancient_mode_ends=0
def col_to_genotype(col):
    stri=''
    c=0
    for r in col.pileups:
        c+=1
        if c > read_max:
            continue
        if r.alignment.mapq < min_mapqual:
            continue
        #print(r.indel)
        if r.is_del:
            #print('TODO del:', r)
            stri += 'D'
        elif r.is_refskip:
            #print('TODO skip:', r)
            stri += 'I'
        elif r.indel < 0:
            #snp db pos is offset like this
            stri += 'D' 
        elif r.indel > 0:
            #snp db pos is offset like this
            stri += 'I'
        else:
            if r.alignment.qual:
                try:
                    q = ord(r.alignment.qual[r.query_position])
                    #print(r)
                except UnicodeDecodeError:
                    return '--'
            else:
                q = 255
            if q >= min_qual:
                base=r.alignment.query_sequence[r.query_position]
                if r.query_position < ancient_mode_ends or (len(r.alignment.query_sequence)-1) - r.query_position < ancient_mode_ends:
                    if r.query_position < ancient_mode_ends:
                        dist=r.query_position
                    else:
                        dist=(len(r.alignment.query_sequence)-1) - r.query_position
                    #if (base == 'T' and not r.alignment.is_reverse) or (base == 'A' and r.alignment.is_reverse):
                    #if (base == 'T' and r.alignment.is_reverse) or (base == 'A' and not r.alignment.is_reverse):
                    if (base == 'T') or (base == 'A'):
                        print('Ancient discard: ', r.query_position, len(r.alignment.query_sequence), dist, base, r.alignment.is_reverse)
                        continue
                stri += base
    #print(stri)
    na=stri.count('A')
    nc=stri.count('C')
    ng=stri.count('G')
    nt=stri.count('T')
    ni=stri.count('I')
    nd=stri.count('D')
    tot=na+nc+nt+ng+ni+nd
    if tot==0:
        print('No read:', col.pos, stri)
        return '--'
    da={'A':na, 'C':nc, 'G':ng, 'T':nt, 'I':ni, 'D':nd}
    das=sorted(da.items(), key = lambda x:(x[1], x[0]), reverse=True)
    first=das[0][1]
    second=das[1][1]
    gt='--'
    if tot >= min_reads:
        gt=das[0][0]+das[0][0]
        hetness = float(second)/first
        if hetness > heterozyg_ratio:
            gt=das[0][0]+das[1][0]
    #print(col.pos, stri, tot, das, hetness, gt)
    #print(col.pos, tot, das, gt) 
    return gt

#interface 1based, internals 0based
def call_range(samfile, contig, pos, num, iter=None):
    ret=['--']*num
    cols = samfile.pileup(contig=contig, start=pos-1, stop=pos+num-1)
    for col in cols:
        if col.reference_pos >= pos-1 and col.reference_pos < pos+num-1:
            if iter and not col.reference_pos+1 in iter:
                #ret.append('??')
                #print('??', col.pos)
                continue
            #print(iter, col.pos)
            r=col_to_genotype(col)
            #ret.append(r)
            ret[col.reference_pos-(pos-1)]=r
    #print(ret)
    return ret

def pos_triplet_37(p):
    global contig
    global lo_37to38
    #TODO 36
    b36='0'
    b37=p
    c38=lo_37to38.convert_coordinate(contig, int(p))
    if c38:
        b38=str(c38[0][1])
    else:
        b38='0'
    return b36, b37, b38

def pos_triplet_38(p):
    global contig
    global lo_38to37
    b36='0'
    c37=lo_38to37.convert_coordinate(contig, int(p))
    #print(p, c37)
    if c37:
        b37=str(c37[0][1])
    else:
        b37='0'
    b38=p
    return b36, b37, b38

def find_ysnps(snpset, samfile):
    global contig
    global pos_triplet_fn
    #call_range(samfile, contig, 13423425, 3)
    #call_range(samfile, contig, 20350749, 3)
    #call_range(samfile, contig, 19730690, 3)
    #call_range(samfile, contig, 5053735, 3)
    #raise
    num_snps=0
    binsize=10000
    binbase=1
    bin=[]
    keys=snp_by_b37.keys()
    keys=sorted(int(i) for i in keys)
    last=keys[-1]
    for pos in keys:
        #if pos > 2930000:
        #    break
        if pos < binbase+binsize:
            bin.append(pos)
        if not pos < binbase+binsize or pos==last:
            #prev bin is rdy
            print('Y', binbase, bin)
            r=call_range(samfile, contig, binbase, binsize, iter=bin)
            #print(r)
            for p in bin:
                if p-binbase < len(r):
                    b36, b37, b38 = pos_triplet_fn(p)
                    gen = r[p-binbase][0]
                    snp = {'id': 'snipsa_%d'%p,
                        'cr': '',
                        'b36': b36,
                        'b37': b37,
                        'b38': b38,
                        'gen': gen }
                    #print(p, r[p-binbase][0])
                    #print(snp)
                    snpset['Y'][p]=snp
                    if gen[0] != '-':
                        num_snps+=1
            
            while not pos < binbase+binsize:
                binbase+=binsize
            bin=[pos]
    return num_snps

str_by_id={}
# Source: http://www.ybrowse.org/gbrowse2/gff/
def load_ystr_db():
    global b3x
    global str_db_file
    #haplo_muts_by_b36 = pickle.load( open( "haploy_map.txt", "rb" ) )
    
    with open(str_db_file, 'r') as f:
        for line in f:
            if len(line) < 3:
                continue
            if line[0] == '#':
                continue
            row = line.split('\t')
            s=row[3]
            e=row[4]
            id=row[8].split(';')[0].split('=')[1]
            strep={b3x:int(s),
                'l': int(e)-int(s)+1,
                'id':id}
            #print(strep)
            str_by_id[id]=strep

#find at least 4 repetitions of 3...30 characters
#str_re = re.compile(r'(.{2,30}?)\1{3,}')
str_re = re.compile(r'([ACGT]{3,30}?)\1{6,}')
def find_ystrs(snpset, samfile):
    global b3x
    global contig
    global pos_triplet_fn
    num_snps=0
    for s in str_by_id:
        if str_by_id[s]['l'] > 500:
            continue
        p=str_by_id[s][b3x]
        seq_list = call_range(samfile, contig, p, str_by_id[s]['l'])
        seq=''
        for c in seq_list:
            seq+=c[0]
        seq=seq.replace('D','')
        #if 'D' in seq:
        #print(p, s, seq)
        
        b36, b37, b38 = pos_triplet_fn(p)
        snp = {'id': 'snipsa_'+s,
            'cr': 'Y',
            'b36': b36,
            'b37': b37,
            'b38': b38,
            'gen': '=0' }
        #snp[b3x]=p

        for m in str_re.finditer(seq):
            num=int(len(m.group(0))/len(m.group(1)))
            print('M:', p, s, m.span(), num, m.group(1), seq)

        m=str_re.search(seq)
        if m:
            num=int(len(m.group(0))/len(m.group(1)))
            snp['gen']='=%d'%num
            if m.group(0)[0]=='-':
                continue
            #print(s, m.group(1), num)
            #print(s, m.group(0), num, m.group(1))
            print(p, s, num, m.group(1), seq)
        else:
            print(s, seq, '-')
            continue
        #print(snp)
        snpset['Y'][p]=snp
        num_snps+=1
    return num_snps





snp_by_mtpos = {}

def load_mtdb():
    with open('haplomt_map.txt', 'r') as f:
        for line in f:
            mut = eval(line)
            snp_by_mtpos[mut['p']] = mut

def find_mtsnps(snpset, samfile):
    global contigmt
    global pos_triplet_fn
    #TODO:check mt
    num_snps=0
    binsize=100
    binbase=1
    bin=[]
    keys=snp_by_mtpos.keys()
    keys=sorted(int(i) for i in keys)
    last=keys[-1]
    for pos in keys:
        #if pos > 3000:
        #    break
        if pos < binbase+binsize:
            bin.append(pos)
        if not pos < binbase+binsize or pos==last:
            #prev bin is rdy
            print('MT', binbase, bin)
            r=call_range(samfile, contigmt, binbase, binsize, iter=bin)
            #print(r)
            for p in bin:
                if p-binbase < len(r):
                    #b36, b37, b38 = pos_triplet_fn(p)
                    gen = r[p-binbase][0]
                    snp = {'id': 'snipsa_%d'%p,
                        'cr': '',
                        'b36': p,
                        'b37': p,
                        'b38': p,
                        'gen': gen }
                    #print(p, r[p-binbase][0])
                    #print(snp)
                    snpset['MT'][p]=snp
                    if gen[0] != '-':
                        num_snps+=1

            while not pos < binbase+binsize:
                binbase+=binsize
            bin=[pos]
    return num_snps





snpauto_by_b37 = {}

def load_snpauto_db():
    with open('snp_db.txt', 'r') as f:
        c=0
        cr=1
        for line in f:
            snp = eval(line)
            if snp['cr'] not in snpauto_by_b37:
                c=0
                print("CR:", snp['cr'])
                snpauto_by_b37[snp['cr']] = {}
            c+=1
            #if c>100:
            #    continue
            snpauto_by_b37[snp['cr']][snp[b3x]] = snp
            
    #print(snpauto_by_b37)
    for cr in snpauto_by_b37.keys():
        print('%s: %d'%(cr, len(snpauto_by_b37[cr])) )

def find_autosnps(snpset, samfile):
    #global contig
    global pos_triplet_fn
    num_snps=0
    for cr in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        if cr not in snpauto_by_b37:
            continue
        lcontig=cr
        binsize=100000
        binbase=1
        bin=[]
        keys=snpauto_by_b37[cr].keys()
        keys=sorted(int(i) for i in keys)
        last=keys[-1]
        for pos in keys:
            #if pos > 2930000:
            #    break
            if pos < binbase+binsize:
                bin.append(pos)
            if not pos < binbase+binsize or pos==last:
                #prev bin is rdy
                print(cr, binbase, bin)
                r=call_range(samfile, lcontig, binbase, binsize, iter=bin)
                #print(r)
                for p in bin:
                    if p-binbase < len(r):
                        b36, b37, b38 = pos_triplet_fn(p)
                        gen = r[p-binbase][0]
                        snp = {'id': 'snipsa_%d'%p,
                            'cr': '',
                            'b36': b36,
                            'b37': b37,
                            'b38': b38,
                            'gen': gen }
                        #print(p, r[p-binbase][0])
                        #print(snp)
                        if cr not in snpset:
                            snpset[cr] = {}
                        snpset[cr][p]=snp
                        if gen[0] != '-':
                            num_snps+=1
                
                while not pos < binbase+binsize:
                    binbase+=binsize
                bin=[pos]

    return num_snps






def get_build(fname):
    build=''
    samfile = pysam.AlignmentFile(fname, get_rtype(fname))
    #print(samfile.text)
    for l in samfile.text.split('\n'):
        if not l.startswith('@SQ'):
            continue
        sn = l.split('\t')[1][3:]
        #print(l.split('\t'))
        #print(sn)
        if sn == '1': build=37
        if sn == 'Y': build=37
        if sn == 'X': build=37
        if sn == 'MT': build=37
        if sn.startswith('chr'): build=38
    #print(build)
    return build

def get_rtype(fname):
    if fname.endswith('.sam'):
        return 'r'
    if fname.endswith('.bam'):
        return 'rb'
    if fname.endswith('.cram'):
        return 'rc'
    raise

def get_index_fname(fname):
    if fname.endswith('.bam'):
        return fname+'.bai'
    if fname.endswith('.cram'):
        return fname+'.crai'
    raise

def index_if_needed(samfname):
    if not os.path.exists(get_index_fname(samfname)):
        print('Creating index for %s'%samfname)
        pysam.index(samfname)

def setup_conv(in_build):
    global b3x
    global str_db_file
    global contig
    global contigmt
    global pos_triplet_fn
    global lo_37to38
    global lo_38to37
    print("Loading LiftOver conversion chain file for build %d..."%in_build)
    if in_build == 19:
        b3x='b37'
        str_db_file='str_hg19.gff3'
        contig='chrY'
        contigmt='chrM'
        pos_triplet_fn = pos_triplet_37
        lo_37to38 = LiftOver('crossmap/GRCh37_to_GRCh38.chain.gz')
    elif in_build == 37:
        b3x='b37'
        str_db_file='str_hg19.gff3'
        contig='Y'
        contigmt='MT'
        pos_triplet_fn = pos_triplet_37
        lo_37to38 = LiftOver('crossmap/GRCh37_to_GRCh38.chain.gz')
    else:
        b3x='b38'
        str_db_file='str_hg38.gff3'
        contig='chrY'
        contigmt='chrM'
        pos_triplet_fn = pos_triplet_38
        lo_38to37 = LiftOver('crossmap/GRCh38_to_GRCh37.chain.gz')

#TODO:temporary interface, redo
convert_ystr=1
convert_y=1
convert_mt=1
convert_snpauto=0
def full_convert(samfname):
    index_if_needed(samfname)
    samfile = pysam.AlignmentFile(samfname, get_rtype(samfname))
       
    snpset={}
    snpset['Y']={}
    snpset['MT']={}
    asnps=0
    ysnps=0
    ysrts=0
    mtnsps=0
    
    if convert_snpauto:
        print("Loading SNP DB...")
        load_snpauto_db()
        print("Reading SNPs...")
        asnps = find_autosnps(snpset, samfile)
        global snpauto_by_b37
        snpauto_by_b37 = {}

    if convert_ystr:
        print("Loading STR DB...")
        load_ystr_db()
        ysrts = find_ystrs(snpset, samfile)

    if convert_y:
        print("Loading Y SNP DB2...")
        if load_ybrowse:
            load_ysnp_ybrowse_db()
        else:
            load_ysnp_dbj()
        print("Reading Y SNPs...")
        ysnps = find_ysnps(snpset, samfile)

    if convert_mt:
        print("Loading MT DB...")
        load_mtdb()
        print("Reading SNPs...")
        mtnsps = find_mtsnps(snpset, samfile)

    samfile.close()
    
    print("Y strs: %d Y snps: %d MT snps: %d"%(ysrts, ysnps, mtnsps))
    
    return snpset

