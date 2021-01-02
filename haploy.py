import re
import csv
import os
#import pickle

def yfind(snpset):
    y=snpset['Y'];
    found_mutations = []
    for snp in y:
        if 'b36' in y[snp] and y[snp]['b36'] in haplo_muts_by_b36:
            if y[snp]['gen'] == haplo_muts_by_b36[y[snp]['b36']]['p']:
                m = haplo_muts_by_b36[y[snp]['b36']]
                found_mutations.append(m)
                #print("Match loc: %s - %s (b38:%s = %s)"%(m['g'], m['m'], m['b38'], m['p']))
        if 'b37' in y[snp] and y[snp]['b37'] in haplo_muts_by_b37:
            if y[snp]['gen'] == haplo_muts_by_b37[y[snp]['b37']]['p']:
                m = haplo_muts_by_b37[y[snp]['b37']]
                found_mutations.append(m)
                #print("Match loc: %s - %s (b38:%s = %s)"%(m['g'], m['m'], m['b38'], m['p']))
        if 'b38' in y[snp] and y[snp]['b38'] in haplo_muts_by_b38:
            if y[snp]['gen'] == haplo_muts_by_b38[y[snp]['b38']]['p']:
                m = haplo_muts_by_b38[y[snp]['b38']]
                found_mutations.append(m)
                #print("Match loc: %s - %s (b38:%s = %s)"%(m['g'], m['m'], m['b38'], m['p']))
    
    found_mutations.sort(key=_mk_sort_key)

    return found_mutations

def _mk_sort_key(mut):
    #sort by mutation age - alphabetical order is roughly age order
    g = mut['g']
    g = g.rstrip('~')
    #some tweaks to improve ordering
    if g == 'CT':
        pr = 'C2'
    elif g.startswith('GH'):
        pr='G2'
    elif g.startswith('IJK'):
        pr = 'I2'
    elif g.startswith('IJ'):
        pr = 'I3'
    elif g.startswith('NO'):
        pr = 'N2'
    else:
        pr = g[0]+'5'
    return pr + g

haplo_muts_by_b38 = {}
haplo_muts_by_b37 = {}
haplo_muts_by_b36 = {}

#Imports database from ISOGG spreadsheet
def load_snp():
    lnum = 0
    pnum = 0
    # https://isogg.org/tree/ISOGG_YDNA_SNP_Index.html
    # Save as csv:
    # https://docs.google.com/spreadsheets/d/1UY26FvLE3UmEmYFiXgOy0uezJi_wOut-V5TD0a_6-bE/edit#gid=1934392066&fvid=105380649
    # Name,Subgroup Name,Alternate Names,rs numbers,Build 37 Number,Build 38 Number,Mutation Info
    with open('SNP Index - Human.csv') as f:
        for line in f:
            lnum += 1
            if lnum <= 2:
                continue
            sline = [ '{}'.format(x) for x in list(csv.reader([line], delimiter=',', quotechar='"'))[0] ]

            if len(sline) < 7:
                continue
            if len(sline[6]) < 3:
                print(line) #TODO: fix other bugs in db
                continue
            try:
                #int(sline[4]) #may be missing
                int(sline[5])
            except:
                print(line) #TODO: ranges are not SNPs
                #try:
                #    sline[5] = re.findall(r'\d+', sline[5])[0]
                #except:
                continue
            #if sline[0] == "PF6503":
            #    continue #TODO!!
            b37=sline[4]
            b38=sline[5]
            #p=sline[6].replace('-.', '->',1).split('->') 
            p=sline[6].split('->') 
            if len(p) < 2:
                print(line) #TODO: fix other bugs in db
                continue
            else:
                p=p[1][0].upper()
            
            mut = {
                'm': sline[0],
                'mall': sline[0],
                'g': sline[1],
                'rs': sline[3],
                'b37': b37,
                'b38': b38,
                'p': sline[6][3],
            }
            if b38 not in haplo_muts_by_b38:
                haplo_muts_by_b38[b38] = mut
            else:
                haplo_muts_by_b38[b38]['mall'] += ", " + sline[0]
            if b37 not in haplo_muts_by_b37:
                haplo_muts_by_b37[b37] = mut
            pnum+=1
    print("Lines in mut DB: %d (%d)"%(lnum, pnum))

# Convert formats with CrossMap and chain file in crossmap/
# http://crossmap.sourceforge.net/
# ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/
def convert_build38to36():
    haplo_mut = haplo_muts_by_b38
    with open('crossmap/conv_in.bed', 'w') as f:
        for mut in haplo_mut:
            f.write("chrY %d %d\n"%(int(haplo_mut[mut]['b38']), int(haplo_mut[mut]['b38'])+1))
    os.system("cd crossmap; CrossMap.py bed GRCh38_to_NCBI36.chain.gz conv_in.bed > conv_out.bed")
    i=0
    with open('crossmap/conv_out.bed', 'r') as f:
        for line in f:
            con=line.split('->')
            if len(con) < 2:
                print(line)
                continue
            b36 = con[1].split()[1]
            #print("Conv %s to %s"%(con[0].split()[1], con[1].split()[1]))
            haplo_muts_by_b36[b36] = haplo_muts_by_b38[con[0].split()[1]]
            haplo_muts_by_b36[b36]['b36'] = b36
            i+=1
    os.system("cd crossmap; rm conv_in.bed conv_out.bed")
    
#Load and save full converted database in a local cache file for faster loading
def save_db():
    #pickle.dump( haplo_muts_by_b36, open( "haploy_map.txt", "wb" ) )
    
    with open('haploy_map.txt', 'w') as f:
        for mut in haplo_muts_by_b36:
            print(haplo_muts_by_b36[mut], file = f)

def load_db():
    #haplo_muts_by_b36 = pickle.load( open( "haploy_map.txt", "rb" ) )
    
    with open('haploy_map.txt', 'r') as f:
        for line in f:
            mut = eval(line)
            haplo_muts_by_b36[mut['b36']] = mut 
            haplo_muts_by_b37[mut['b37']] = mut 
            haplo_muts_by_b38[mut['b38']] = mut 








