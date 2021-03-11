import re
import csv
import os
import copy
import snpload
from bs4 import BeautifulSoup
import urllib.request

def print_uptree(snpset, ut, do_print=True, b3x='b37'):
    rep=''
    y=snpset['Y'];
    pos=0
    neg=0
    tot=0
    for i, mut in enumerate(ut):
        txt=''
        if 'txt' in mut:
            txt=mut['txt']
        if mut[b3x] in y:
            rep += "%-1s%-9s%s %-22s %-38s %s\n"%(mut['tag'], mut['g'], y[mut[b3x]]['gen'], mut['isog'], mut['raw'], txt) 
        else:
            rep += "%-1s%-9s%s %-22s %-38s %s\n"%(mut['tag'], mut['g'], ' ', mut['isog'], mut['raw'], txt)
            pass
    if do_print:
        print(rep)
    return rep

def print_extras(snpset, bt, do_print=True):
    rep=''
    last=''
    out=[]
    outs=''
    for e in bt['extras']:
        out.append(e['raw'])
        outs += e['raw'] + ' '
        last=e['raw']
    rep += 'Extra: '+outs+'\n'
    if do_print:
        print(rep)
    return rep

#TODO
def print_all(snpset, bt, do_print=True):
    rep=''
    last=''
    out=[]
    outs=''
    c=0
    for e in bt['all']:
        out.append(e['raw'])
        outs += '%-40s'%(e['raw']+',')
        c+=1
        if c == 3:
            outs += '\n'
            c=0
        last=e['raw']
    rep += 'All:\n'+outs+'\n'
    if do_print:
        print(rep)
    return rep

def print_data(do_print=True):
    rep=''

    rep += 'Based on data from yfull.com on 2020-03-06 (CC-BY) and isogg.org/tree Y-DNA Tree 2019-2020 and ybrowse.org (CC BY-NC-SA)\n'

    if do_print:
        print(rep)
    return rep

def path_str(ut, n):
    rep=''

    c=0
    prev=''
    #for i, mut in enumerate(ut):
    for mut in reversed(ut):
        if mut['g'] != prev:
            if c > 0:
                rep = ' -> ' + rep
            else:
                rep = ' (ISOGG: %s)'%mut['isog']
            rep = mut['g'] + rep
            c+=1
        prev = mut['g']
        if c >= n:
            break
        if not '-' in mut['g']:
            break

    return rep

def report(fname, n, do_uptree=True, do_extra=True, do_all=False, filt='', force=''):
    rep=''
    snpset, meta = snpload.load(fname, ['Y'])
    if 'Y' not in snpset:
        return "No Y data found\n"

    rep += print_data(False)

    rep += "%s: Total SNPs: %d\n"%(fname, meta['total'])

    b3x='b36'
    if meta['build']==37:
        b3x='b37'
    if meta['build']==38:
        b3x='b38'
    best_trees = yfind2(snpset, n, filt, force, b3x)

    for bt in best_trees:
        if do_all:
            rep += print_all(snpset, bt, False)

        if do_uptree:
            rep += print_uptree(snpset, bt['ut'], False, b3x)

        leaf_mut = bt['ut'][len(bt['ut'])-1]
        #rep += "Result (%-8s %5.1f%% -%d +%d): %-8s (ISOGG: %s)\n"%(leaf_mut['raw'], bt['score'], bt['neg'], len(bt['extras']), leaf_mut['g'], leaf_mut['isog'])
        rep += "Result (%-8s %5.1f%% -%d +%d): %-8s\n"%(leaf_mut['raw'], bt['score'], bt['neg'], len(bt['extras']), leaf_mut['g'])
        #rep += "%s (ISOGG: %s)\n"%(path_str(bt['ut'], 15), leaf_mut['isog'])
        rep += "%s\n"%(path_str(bt['ut'], 20))

        if do_extra:
            rep += print_extras(snpset, bt, False)
    return rep


#Create a list of mutations on a path upwards from a mutation
def yfind_uptree(snpset, find_g):
    found = 0
    uptree = []
    sames = []
    sameg = ''
    for mut in reversed(haplo_muts_list):
        if found == 0:
            if mut['g'] == sameg:
                sames.append(mut)
            else:
                sameg = mut['g']
                sames = []
                sames.append(mut)
            if mut['g'] == find_g: #mut_leaf['g']: #TODO: g?
                found = 1
                depth = mut['l']
                g = mut['g']
                uptree.extend(reversed(sames))
        else:
            if mut['g'] == g:
                uptree.insert(0, mut)
                continue
            if mut['l'] < depth:
                depth = mut['l']
                g = mut['g']
                uptree.insert(0, mut)
    return uptree

def ymatch(snpset, mut, b3x):
    y=snpset['Y'];
    if mut['t'] == y[mut[b3x]]['gen']:
        return True
    else:
        return False

# Finds n best Y-DNA matches by tracing all possible mutation paths in database
def yfind2(snpset, nbest=5, filt='', force='', b3x='b37'):
    mt=snpset['Y'];
    all_mutations = []
    uptrees=[]
    allmuts=[]
    #TODO: support more exotic mutations?

    #find the route uptree for each mutation in snpset
    for mut in haplo_muts_list:
        pos = mut[b3x]
        if pos in mt:
            #print(mt[pos], mut)
            if mt[pos]['gen'] == mut['t']:
                #print(mt[pos], mut)
                #if mut['!']==0:
                all_mutations.append(mut)
                uptree = yfind_uptree(snpset, mut['g'])
                #print(uptree)
                #print_uptree(snpset, uptree, do_print=True, b3x=b3x)
                uptrees.append(uptree)
                allmuts.append(mut)

    #if there is a forced path set, find a path for it too
    if force:
        uptree = yfind_uptree(snpset, force)
        uptrees.append(uptree)

    #find the uptree route that is most consistent
    #random mutations can have many negative matches above them until the path reaches the common
    #common segment with the correct ones
    best_trees=[]
    for ut in uptrees:
        verbose=0
        pos=0
        neg=0
        tot=0
        bm=[]
        dbm=[]
        tbm=[]
        #prepare list of back mutations that need special handling
        #for mut in ut:
        #    #TODO triple bm
        #    if mut['!'] == 1:
        #        bm.append(mut)
        #    if mut['!'] == 2:
        #        dbm.append(mut)
        #    if mut['!'] == 3:
        #        tbm.append(mut)
        tag=''
        ut_copy=[]
        for i, mut in enumerate(ut):
            if mut[b3x] in mt:
                bm_matched = 0 #bm_match(snpset, i, ut, bm)
                dbm_matched = 0 #bm_match(snpset, i, ut, dbm)
                if ymatch(snpset, mut, b3x):
                    if bm_matched:
                        ##Double back mutation
                        if verbose: print("!!%-12s %s %s"%(mut['g'], mt[mut['p']]['gen'], mut['raw']))
                        tag='!!'
                    else:
                        #TODO: should not match double without a single back mutation
                        if verbose: print("+ %-12s %s %s"%(mut['g'], mt[mut['p']]['gen'], mut['raw']))
                        tag='+'
                        pos+=1
                        tot+=1
                else:
                    if bm_matched:
                        if verbose: print("! %-12s %s %s"%(mut['g'], mt[mut['p']]['gen'], mut['raw']))
                        tag='!'
                    elif dbm_matched:
                        if verbose: print("!!%-12s %s %s"%(mut['g'], mt[mut['p']]['gen'], mut['raw']))
                        tag='!!'
                    else:
                        if verbose: print("- %-12s %s %s"%(mut['g'], mt[mut['p']]['gen'], mut['raw']))
                        tag='-'
                        neg+=1
                        tot+=1
            else:
                tag=''
                #if verbose: print(" (%s)"%mut['g'])
                pass
            mut_copy = mut.copy()
            mut_copy['tag']=tag
            ut_copy.append(mut_copy)

        #find extras: mutations that are found in snpset but not in uptree
        extras=allmuts.copy()
        toremove={}
        for mut in ut_copy:
            toremove[mut['raw']]=1
        extras = filter(lambda e: not e['raw'] in toremove, extras)
        extras = sorted(extras, key=lambda i: int(i[b3x]))
        last = ''
        de_duplicate=[]
        for e in extras:
            if last != e['raw']: # and e['!'] == 0:
                de_duplicate.append(e)
                last = e['raw']
        extras = de_duplicate
        nextras = len(extras)

        #extract unique all mutations
        last = ''
        de_duplicate=[]
        #for e in sorted(all_mutations, key=lambda x: int(re.compile(r'[^\d.]+').sub('', x['raw']))):
        for e in sorted(all_mutations, key=lambda x: x['raw']):
            if last != e['raw']:
                de_duplicate.append(e)
                last = e['raw']
        all_mutations = de_duplicate

        #calculate score giving penalty from negative matches, favoring longest matches
        score = 100.0*(pos - 2.0*neg + 0.00*tot - 0.0*nextras)/(tot + nextras)

        bt={
            'ut': ut_copy,
            'score': score,
            'extras': extras,
            'all': all_mutations,
            'pos': pos,
            'neg': neg,
            'tot': tot,
        }
        best_trees.append(bt)

    if filt:
        if filt[0] == '=':
            best_trees = filter(lambda bt: filt[1:] == bt['ut'][-1]['g'], best_trees)
        else:
            best_trees = filter(lambda bt: filt in bt['ut'][-1]['g'], best_trees)
    best_trees=sorted(best_trees, key=lambda i: -i['score'])

    return best_trees[:nbest]

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
haplo_muts_by_name = {}
haplo_muts_list = []

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
                print('TODO1:', line) #TODO: fix other bugs in db
                continue
            try:
                #int(sline[4]) #may be missing
                int(sline[5])
            except:
                print('TODO2:', line) #TODO: ranges are not SNPs
                #try:
                #    sline[5] = re.findall(r'\d+', sline[5])[0]
                #except:
                continue
            b37=sline[4]
            b38=sline[5]
            mname=sline[0].split('.')[0] ##TODO test multiple hg with same mut
            #p=sline[6].replace('-.', '->',1).split('->') 
            p=sline[6].split('->') 
            if len(p) < 2:
                print('TODO:', line) #TODO: fix other bugs in db
                continue
            else:
                p=p[1][0].upper()
            mut = {
                'm': mname,
                'mall': sline[0],
                'g': sline[1],
                'rs': sline[3],
                'b37': b37,
                'b38': b38,
                'p': p,
            }
            #print(mname)
            if mname not in haplo_muts_by_name:
                haplo_muts_by_name[mname] = mut
            if b38 not in haplo_muts_by_b38:
                haplo_muts_by_b38[b38] = mut
            else:
                haplo_muts_by_b38[b38]['mall'] += ", " + sline[0]
            if b37 not in haplo_muts_by_b37:
                haplo_muts_by_b37[b37] = mut
            pnum+=1
    print("Lines in ISOGG mut DB: %d (%d)"%(lnum, pnum))

haplo_yfull_muts_by_name = {}
haplo_yfull_muts_by_b38 = {}

#Imports database from YFull /snp-list/ URL
def load_yfull_snp_file(fname):
    pnum = 0
    with open(fname) as f:
        print('Importing file: ' +fname)
        soup = BeautifulSoup(f.read(), features="html.parser")
        d = soup.find('div', id='t1')
        t = d.find('table')
        rows = t.find_all('tr')
        for row in rows:
            pnum+=1
            if pnum == 1:
                continue
            cells = row.find_all('td')
            mname=cells[0].text
            b38=cells[3].text
            if b38 == '':
                continue
            mut = {
                'm': mname,
                'mall': '?',
                'g': '?',
                'rs': '?',
                'b37': cells[2].text,
                'b38': b38,
                'p': cells[5].text,
            }
            if mname not in haplo_yfull_muts_by_name:
                haplo_yfull_muts_by_name[mname] = mut
            if b38 not in haplo_yfull_muts_by_b38:
                haplo_yfull_muts_by_b38[b38] = mut

def load_yfull_snp(pages):
    try:
        os.mkdir('yfull')
    except OSError:
        pass

    for page in range(1, pages):
        fname='yfull/yfull-snp-'+str(page)+'.html'
        url='https://www.yfull.com/snp-list/?page='+str(page)
        try:
            with open(fname) as f:
                pass
        except OSError:
            print('File not found: ' +fname)
            print('Downloading ' + url + 'to file: ' + fname)
            urllib.request.urlretrieve(url, fname);

        load_yfull_snp_file(fname)

    print("Lines in YFull mut DB: ", len(haplo_yfull_muts_by_name))

haplo_ybrowse_muts_by_name = {}
haplo_ybrowse_muts_by_b38 = {}

# Source: http://www.ybrowse.org/gbrowse2/gff/
def load_ybrowse_snp():
    ln=0;
    with open('snps_hg38.gff3') as f:
        for line in f:
            ln+=1
            if ln % 100000 == 0:
                print('Line: %d...'%ln)
            if len(line) < 2:
                continue
            if line[0] == '#':
                continue
            tlv=line.split('\t')
            if len(tlv) < 8:
                print('TODO:', tlv)
                continue
            #print(tlv)
            if tlv[0] != 'chrY' or  tlv[2] != 'snp' or (tlv[1] != 'point' and tlv[1] != 'indel'):
                if tlv[1] != 'primate':
                    print(tlv)
                continue
            #TODO quality based filtering
            b38 = tlv[3]
            f=tlv[8].split(';')
            mname=f[0].split('=')[1]
            der=f[3].split('=')[1]
            der=der[0].upper()
            isog='' #f[7].split('=')[1]
            #print(b38, mname, der)
            mut = {
                'm': mname,
                'mall': '?',
                'g': isog,
                'rs': '?',
                'b38': b38,
                'p': der,
            }
            if mname not in haplo_ybrowse_muts_by_name:
                haplo_ybrowse_muts_by_name[mname] = mut
            if b38 not in haplo_ybrowse_muts_by_b38:
                haplo_ybrowse_muts_by_b38[b38] = mut
    print("Lines in YBrowse snp DB: ", len(haplo_ybrowse_muts_by_name))

# Convert formats with CrossMap and chain file in crossmap/
# http://crossmap.sourceforge.net/
# ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/
def convert_build38_mkinput():
    haplo_mut = haplo_muts_by_name
    haplo_mut2 = haplo_yfull_muts_by_name
    haplo_mut3 = haplo_ybrowse_muts_by_name
    with open('crossmap/conv_in.bed', 'w') as f:
        for mut in haplo_mut:
            f.write("chrY %d %d %s\n"%(int(haplo_mut[mut]['b38']), int(haplo_mut[mut]['b38']), haplo_mut[mut]['m']))
        for mut in haplo_mut2:
            f.write("chrY %d %d %s\n"%(int(haplo_mut2[mut]['b38']), int(haplo_mut2[mut]['b38']), haplo_mut2[mut]['m']))
        for mut in haplo_mut3:
            f.write("chrY %d %d %s\n"%(int(haplo_mut3[mut]['b38']), int(haplo_mut3[mut]['b38']), haplo_mut3[mut]['m']))

def convert_build38to36():
    convert_build38_mkinput()
    os.system("cd crossmap; CrossMap.py bed GRCh38_to_NCBI36.chain.gz conv_in.bed > conv_out.bed")
    i=0
    with open('crossmap/conv_out.bed', 'r') as f:
        for line in f:
            con=line.split('->')
            if len(con) < 2:
                #print('TODO:', line)
                continue
            b36 = con[1].split()[1]
            b38 = con[0].split()[1]
            mname=''
            if len(con[0].split()) > 3:
                mname = con[0].split()[3]
            if mname in haplo_muts_by_name:
                haplo_muts_by_name[mname]['b36'] = b36
            if mname in haplo_yfull_muts_by_name:
                haplo_yfull_muts_by_name[mname]['b36'] = b36
            if mname in haplo_ybrowse_muts_by_name:
                haplo_ybrowse_muts_by_name[mname]['b36'] = b36
            i+=1
    os.system("cd crossmap; rm conv_in.bed conv_out.bed")

def convert_build38to37():
    convert_build38_mkinput()
    os.system("cd crossmap; CrossMap.py bed GRCh38_to_GRCh37.chain.gz conv_in.bed > conv_out.bed")
    i=0
    with open('crossmap/conv_out.bed', 'r') as f:
        for line in f:
            con=line.split('->')
            if len(con) < 2:
                #print('TODO:', line)
                continue
            b37 = con[1].split()[1]
            b38 = con[0].split()[1]
            mname=''
            if len(con[0].split()) > 3:
                mname = con[0].split()[3]
            if mname in haplo_muts_by_name:
                haplo_muts_by_name[mname]['b37'] = b37
            if mname in haplo_yfull_muts_by_name:
                haplo_yfull_muts_by_name[mname]['b37'] = b37
            if mname in haplo_ybrowse_muts_by_name:
                haplo_ybrowse_muts_by_name[mname]['b37'] = b37
            i+=1
    os.system("cd crossmap; rm conv_in.bed conv_out.bed")
    
#Load and save full converted database in a local cache file for faster loading
def save_db():
    #pickle.dump( haplo_muts_by_b36, open( "haploy_map.txt", "wb" ) )
    
    with open('haploy_map.txt', 'w') as f:
        for mut in haplo_muts_by_name:
            print(haplo_muts_by_name[mut], file = f)

def load_db():
    #haplo_muts_by_b36 = pickle.load( open( "haploy_map.txt", "rb" ) )
    
    with open('haploy_map.txt', 'r') as f:
        for line in f:
            mut = eval(line)
            haplo_muts_by_b36[mut['b36']] = mut 
            haplo_muts_by_b37[mut['b37']] = mut 
            haplo_muts_by_b38[mut['b38']] = mut 

def save_yfull_db():
    with open('yfull_snp.txt', 'w') as f:
        for mut in haplo_yfull_muts_by_name:
            print(haplo_yfull_muts_by_name[mut], file = f)

def save_ybrowse_db():
    with open('ybrowse_snp.txt', 'w') as f:
        for mut in haplo_ybrowse_muts_by_name:
            print(haplo_ybrowse_muts_by_name[mut], file = f)

def save_db2():
    with open('haploy_map2.txt', 'w') as f:
        for mut in haplo_muts_list:
            print(mut, file = f)

def load_db2():
    with open('haploy_map2.txt', 'r') as f:
        for line in f:
            mut = eval(line)
            haplo_muts_list.append(mut)

def show_db2():
    for m in sorted(haplo_muts_list, key=lambda e: int(e['b37'])):
        print(m)


def decode_entry(e):
    global haplo_muts_by_name
    e2=[]
    for e1 in e.split('/'):
        e2.append(e1.replace('(H)',''))
    m={}
    for e1 in e2:
        if e1 in haplo_ybrowse_muts_by_name:
            mut = haplo_ybrowse_muts_by_name[e1]
            m['f']=''
            if not 'isogg' in m:
                m['isog']=''
            #m['isog']+=mut['g']+'(YB)'
            m['t']=mut['p']
            m['b38']=mut['b38']
            if 'b37' in mut:
                m['b37']=mut['b37']
            if 'b36' in mut:
                m['b36']=mut['b36']
            break
    for e1 in e2:
        if e1 in haplo_yfull_muts_by_name:
            mut = haplo_yfull_muts_by_name[e1]
            m['f']=''
            if not 'isogg' in m:
                m['isog']=''
            m['t']=mut['p']
            m['b38']=mut['b38']
            if 'b37' in mut:
                m['b37']=mut['b37']
            if 'b36' in mut:
                m['b36']=mut['b36']
            break
    for e1 in e2:
        if e1 in haplo_muts_by_name:
            mut = haplo_muts_by_name[e1]
            m['f']=''
            #if not 'isogg' in m:
            #    m['isog']=''
            m['isog']=''
            m['isog']+=mut['g']+''
            m['t']=mut['p']
            m['b38']=mut['b38']
            if 'b37' in mut:
                m['b37']=mut['b37']
            if 'b36' in mut:
                m['b36']=mut['b36']
            break
    return m

def yfull_fname(group):
    if group:
        return 'yfull/yfull-ytree-'+group+'.html'
    else:
        return 'yfull/yfull-ytree.html'

def yfull_url(group):
    if group:
        return 'https://www.yfull.com/tree/' + group + '/'
    else:
        return 'https://www.yfull.com/tree/'

# YFull mtree import (experimental)
def download_yfull_file(group):
    try:
        os.mkdir('yfull')
    except OSError:
        pass
    fname = yfull_fname(group)
    url = yfull_url(group)
    print('Downloading ' + url + 'to file: ' + fname)
    urllib.request.urlretrieve("https://www.yfull.com/tree/"+group+"/", fname);

def yfull_parse_muts(li):
    s=''
    snpforhg=li.find('span', class_='yf-snpforhg', recursive=False)
    if snpforhg:
        s+=snpforhg.text
    plussnps=li.find('span', class_='yf-plus-snps', recursive=False)
    if plussnps:
        s += ' * ' + plussnps['title']
    o=[]
    if len(s) > 0:
        for m in s.split('*'):
            o.append(m.strip())
    return o

def yfull_parse_age(li):
    s=''
    agespan=li.find('span', class_='yf-age', recursive=False)
    if agespan:
        s+=agespan.text
    return s

def yfull_is_tree_quirk(group_name, fileroot):
    if fileroot:
        return False
    if group_name=='R-P312':
        return True
    if group_name=='R-Z2118':
        return True
    return False

def yfull_recurse_list(ul_in, level, fileroot):
    lis = ul_in.find_all('li', recursive=False)
    for li in lis:
        muts={}
        muts['l']=level
        g=li.find('a', recursive=False)
        group_name=''
        if g:
            group_name=g.text
            muts['g']=g.text
        l=li.find('a', href=True, recursive=False)
        if l:
            muts['link']=l['href']
            #print(g.text, g['href'])

        mlist = yfull_parse_muts(li)
        age = yfull_parse_age(li)

        #if 'g' in muts and not muts['g'].endswith('*') and not fileroot:
        if 'g' in muts and not muts['g'].endswith('*'):
            #add separate entries for each mut in same hg
            #print(mlist)
            for m in mlist:
                mutse=dict(muts)
                if 'link' in mutse:
                    del mutse['link']
                dec = decode_entry(m)
                mutse['raw']=m
                if not dec:
                    print('No pos found for', m)
                    global no_pos_counter
                    no_pos_counter+=1
                #    dec['isog']='n/a'
                #    dec['t']='?'
                #    dec['b38']='0'
                #    dec['b37']='0'
                #    dec['b36']='0'
                    continue
                if not 'b37' in dec or not dec['b37']:
                    print('TODO b37:', mutse, dec)
                    dec['b37'] = '0'
                    continue #TODO: include also build38-only snps
                if not 'b38' in dec or not dec['b38']:
                    print('TODO b38:', mutse, dec)
                    dec['b38'] = '0'
                if not 'b36' in dec or not 'b36' in dec:
                    print('TODO b36:', mutse, dec)
                    dec['b36'] = '0'
                #muts['f']=dec['f']
                mutse['t']=dec['t']
                mutse['isog']=dec['isog']
                mutse['b36']=dec['b36']
                mutse['b37']=dec['b37']
                mutse['b38']=dec['b38']
                mutse['raw']=m
                if age:
                    mutse['txt']=age
                #print(mutse)
                haplo_muts_list.append(mutse)

        ul = li.find('ul', recursive=False)
        if ul and not yfull_is_tree_quirk(group_name, fileroot):
            #print('->')
            yfull_recurse_list(ul, level+1, False)
            #print('<-')
        else:
            if 'g' in muts and muts['g'].endswith('*'):
                continue
            if 'link' in muts:
                group=muts['link'].split('/')[-2]
                #print('FILE: ' +fname)
                yfull_recurse_file(group, level)
                #print('END: ' +fname)
    return 0

def yfull_recurse_file(group, level):
    fname = yfull_fname(group)
    try:
        with open(fname) as f:
            pass
    except OSError:
        print('File not found: ' +fname)
        download_yfull_file(group)

    with open(fname) as f:
        print('Importing file: ' +fname)
        soup = BeautifulSoup(f.read(), features="html.parser")
        ul = soup.find('ul', id='tree')
        yfull_recurse_list(ul, level, True)

no_pos_counter=0
def import_yfull_tree():
    #TODO complex muts
    #yfull_recurse_file('', 0)       #cannot use top level: muts missing
    yfull_recurse_file('A00', 0)
    yfull_recurse_file('A0-T', 0)
    print('No pos found for %d tree nodes'%no_pos_counter)
    print('Tree database size is %d nodes'%len(haplo_muts_list))






