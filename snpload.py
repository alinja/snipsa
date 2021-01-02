import re
import zipfile
import gzip
from shutil import copyfile
import shutil

def load(fname, snpset, crs=[]):
    
    tmpfile = preprocess_file(fname)
    (build, fmt) = detect_file_format(tmpfile)
    
    if fmt == '23andme':
        importer = import_line_23andme
    elif fmt == 'ancestry':
        importer = import_line_ancestry
    elif fmt == 'ftdna':
        importer = import_line_ftdna
    else:
        print("Undetected format: build%d, %s"%(build, fmt))
        return 0

    
    n_total = 0
    with open(tmpfile) as f:
        for line in f:
                
            try:
                (snp, pos, cr, gen) = importer(line, build)
            except:
                continue
            
            if gen[0] == "-":
                continue
            if cr == "Y":
                if len(gen) > 1:
                    snp['gen'] = gen[0]
                
            if cr in crs or len(crs) == 0:
                if cr not in snpset:
                    snpset[cr] = {}
                #if pos in snpset[cr]:
                #    if not snp['id'].startswith("i"):
                #        #print("dup ", snpset[cr][pos], snp)
                #        pass
                snpset[cr][pos]=snp
                n_total+=1
                continue

    if build != 36 and build != 37 and build != 38:
        print("BUILD NOT SUPPPORTED!!!!! = %d"%build)
        
    return n_total

def import_line_23andme(line, build):
    if len(line) < 7:
        raise Exception()
    if line.startswith('#'):
        raise Exception()
    sline = line.split()
    if len(sline) < 4:
        raise Exception()
    cr = sline[1];
    pos = sline[2];
    gen = sline[3];
        
    snp = {'id': sline[0],
        'cr': cr,
        'gen': gen }

    if build == 36:
        snp['b36'] = pos
    elif build == 37:
        snp['b37'] = pos
    elif build == 38:
        snp['b38'] = pos
        
    return (snp, pos, cr, gen)

def import_line_ancestry(line, build):
    if 'allele1' in line and 'allele2' in line:
        raise Exception()
    if len(line) < 7:
        raise Exception()
    if line.startswith('#'):
        raise Exception()
    sline = line.split()
    if len(sline) < 4:
        raise Exception()
    cr = sline[1];
    if cr == '23':
        cr = 'X'
    elif cr == '24':
        cr = 'Y'
    elif cr == '25':
        cr = 'YAUTO'
    elif cr == '26':
        cr = 'MT'

    pos = sline[2];
    gen = sline[3]+sline[4];
        
    snp = {'id': sline[0],
        'cr': cr,
        'gen': gen }

    if build == 36:
        snp['b36'] = pos
    elif build == 37:
        snp['b37'] = pos
    elif build == 38:
        snp['b38'] = pos
        
    return (snp, pos, cr, gen)

def import_line_ftdna(line, build):
    if len(line) < 7:
        raise Exception()
    if line.startswith('RSID'):
        raise Exception()
    sline = line.split(',')
    if len(sline) < 4:
        raise Exception()
    cr = sline[1].strip('"');
    pos = sline[2].strip('"');
    gen = sline[3].strip().strip('"');
        
    snp = {'id': sline[0].strip('"'),
        'cr': cr,
        'gen': gen }

    if build == 36:
        snp['b36'] = pos
    elif build == 37:
        snp['b37'] = pos
    elif build == 38:
        snp['b38'] = pos
        
    return (snp, pos, cr, gen)

def detect_file_format(fname):
    lc=0
    build=0
    fmt=''
    with open(fname) as f:
        for line in f:
            lc += 1
            if lc > 30:
                break
            if line.startswith('#'):
                if 'build ' in line:
                    build = int(re.findall(r'build \d+', line)[0].split()[1])
                if "23andMe" in line:
                    fmt = '23andme'
                if "AncestryDNA" in line:
                    fmt = 'ancestry'
                if "##fileformat=VCF" in line:
                    fmt = 'vcf'
                    return (build, fmt)
            if lc == 1 and line.startswith('RSID'):
                     fmt = 'ftdna'
            if fmt == 'ftdna':
                sline=line.split(',')
                if sline[0] == '"rs6681049"' and  sline[2] == '"800007"':
                    build = 37
                if sline[0] == '"rs6681049"' and  sline[2] == '"789870"':
                    build = 36
    #print("Detected format: build%d, %s"%(build, fmt))
    return (build, fmt)

def is_gz_file(fname):
    with open(fname, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def preprocess_file(fname):
    #TODO make real temp file
    tmpfile = 'genome_data.tmp'
    if zipfile.is_zipfile(fname):
        with zipfile.ZipFile(fname) as z:
            zfname = z.namelist()[0]
            #print('ZIP input file: %s'%zfname)
            with z.open(zfname) as zf, open(tmpfile, 'wb') as f:
                shutil.copyfileobj(zf, f)
                
    elif is_gz_file(fname):
        #print('gzip input file: %s'%fname)
        with gzip.open(fname, 'rb') as f_in:
            with open(tmpfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        #print('ASCII input file: %s'%fname)
        shutil.copyfile(fname, tmpfile);
    return tmpfile

def index_by_rs(snpset):
    snpseto = {}

    for cr in snpset:
        if cr not in snpseto:
            snpseto[cr] = {}
        for snp in snpset[cr]:
            snpseto[cr][ snpset[cr][snp]['id'] ] = snpset[cr][snp]
    return snpseto



 
 


def show_stats(snpset):
    n_total = 0
    for cr in snpset:
        n_smps = len(snpset[cr])
        print("Chromosome %s: %s SNPs"%(cr, n_smps))
        n_total += n_smps
    print("Total SNPs: %d"%n_total)

def show_gts(snpset):
    d={}
    for cr in snpset:
        for snp in snpset[cr]:
            gt = snpset[cr][snp]['gen']
            if gt not in d:
                d[gt] = 1
            else:
                d[gt] += 1
            
    print(d)



