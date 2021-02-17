import re
import csv
import os
from bs4 import BeautifulSoup
import urllib.request
import re
import copy

def print_uptree(snpset, ut):
    mt=snpset['MT'];
    pos=0
    neg=0
    tot=0
    for i, mut in enumerate(ut):
        if mut['p'] in mt:
            print("%-3s%-12s %s %s"%(mut['tag'], mut['g'], mt[mut['p']]['gen'], mut['raw'])) 
        else:
            print("%-3s%-12s %s %s"%(mut['tag'], mut['g'], ' ', mut['raw'])) 
            pass

def print_extras(snpset, bt):
    last=''
    out=[]
    outs=''
    for e in bt['extras']:
        out.append(e['raw'])
        outs += e['raw'] + ' '
        last=e['raw']
    print('Extra: '+outs)

#Create a list of mutations on a path upwards from a mutation 
def mtfind_uptree(snpset, mut_leaf):
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
            if mut == mut_leaf:
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

def bm_match(snpset, i, ut, bm):
    mt=snpset['MT'];
    if i == len(ut)-1:
        return False
    mut=ut[i]
    bm_match=0
    for bmc in bm:
        if bmc['p'] in mt:
            #print(mut, bmc)
            if mut['p'] == bmc['p'] and mut['t'] != bmc['t']: 
                return True

    return False
        
def mtmatch(snpset, mut):
    mt=snpset['MT'];
    if mut['t'] == mt[mut['p']]['gen']:
        return True
    else:
        return False

# Finds n best mtDNA matches by tracing all possible mutation paths in database 
def mtfind(snpset, nbest=5):
    mt=snpset['MT'];
    all_mutations = []
    uptrees=[]
    allmuts=[]
    #TODO: support more exotic mutations?
    
    #find the route uptree for each mutation in snpset
    for mut in haplo_muts_list:
        pos = mut['p']
        if pos in mt:
            if mt[pos]['gen'] == mut['t']:
                #print(mt[pos], mut)
                all_mutations.append(mut)
                uptree = mtfind_uptree(snpset, mut)
                uptrees.append(uptree)
                allmuts.append(mut)

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
        for mut in ut:
            #TODO triple bm
            if mut['!'] == 1:
                bm.append(mut)
            if mut['!'] == 2:
                dbm.append(mut)
            if mut['!'] == 3:
                tbm.append(mut)
        tag=''
        ut_copy=[]        
        for i, mut in enumerate(ut):
            if mut['p'] in mt:
                bm_matched = bm_match(snpset, i, ut, bm)
                dbm_matched = bm_match(snpset, i, ut, dbm)
                if mtmatch(snpset, mut):
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
        extras = sorted(extras, key=lambda i: int(i['p']))
        last = ''
        de_duplicate=[]
        for e in extras:
            if last != e['raw'] and e['!'] == 0:
                de_duplicate.append(e)
                last = e['raw']
        extras = de_duplicate
        nextras = len(extras)
        
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
    
    best_trees=sorted(best_trees, key=lambda i: -i['score'])   

    return best_trees[:nbest]

#
# Database importing
#

haplo_muts_list = []

non_decimal_re = re.compile(r'[^\d.]+')
def decode_entry(e):
    m={}
    bm=0
    e=e.strip('(').strip(')')
    if e.endswith('!!!'):
        bm=3
    elif e.endswith('!!'):
        bm=2
    elif e.endswith('!'):
        bm=1
    e=e.strip('!')
    m['f']=e[0].upper()
    m['t']=e[-1].upper()
    m['p']=non_decimal_re.sub('', e)
    m['!']=bm
    return m
    
# Imports database from PhyloTree mtDNA tree
# file source: https://www.phylotree.org/builds/mtDNA_tree_Build_17.zip
def import_snp():
    with open("mtDNA tree Build 17.htm", encoding="windows-1252") as f:
        soup = BeautifulSoup(f.read(), features="html.parser")
        table = soup.find("table")
        output_rows = []
        rc=0
        for table_row in table.findAll('tr'):
            rc+=1
            if rc < 19:
                continue
            columns = table_row.findAll('td')
            st=0
            level=0
            muts={}
            cell=[]
            for column in columns:
                level+=1
                cell=column.text.split()
                if st == 0 and len(cell) != 0:
                    #print(level, cell[0])
                    muts['l']=level
                    muts['g']=cell[0]
                    st=1
                    continue
                if st == 1 and len(cell) != 0:
                    break
            if not 'g' in muts:
                continue

            if len(cell) > 0:
                if cell[0][-1].strip(')').strip('!').isdigit():
                    muts['l']=muts['l']-1
                    cell=[muts['g']]
                    print(muts, cell)

            #print(muts, cell)
            for m in cell:
                if m.endswith(')'):
                    continue #TODO!!!
                mutse=dict(muts)
                dec = decode_entry(m)
                #muts['f']=dec['f']
                mutse['t']=dec['t']
                mutse['p']=dec['p']
                mutse['!']=dec['!']
                mutse['raw']=m
                haplo_muts_list.append(mutse)
    
# Load and save full converted database in a local cache file for faster loading
def save_db():
    with open('haplomt_map.txt', 'w') as f:
        for mut in haplo_muts_list:
            print(mut, file = f)

def load_db():
    with open('haplomt_map.txt', 'r') as f:
        for line in f:
            mut = eval(line)
            haplo_muts_list.append(mut)

def show_db():
    for m in sorted(haplo_muts_list, key=lambda e: int(e['p'])):
        print(m)



# YFull mtree import (experimental)
def download_yfull_file(group):
    try:
        os.mkdir('yfull')
    except OSError:
        pass
    fname='yfull/yfull-mtree-'+group+'.html'
    print('Downloading file: ' + fname)
    urllib.request.urlretrieve("https://www.yfull.com/mtree/"+group+"/", fname);

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

def yfull_recurse_list(ul_in, level, fileroot):
    lis = ul_in.find_all('li', recursive=False)
    for li in lis:
        muts={}
        muts['l']=level
        g=li.find('a', href=True)
        if g:
            muts['g']=g.text
            muts['link']=g['href']
        mlist = yfull_parse_muts(li)
        
        if 'g' in muts and not muts['g'].endswith('*') and not fileroot:
            for m in mlist:
                mutse=dict(muts)
                del mutse['link']
                dec = decode_entry(m)
                #muts['f']=dec['f']
                mutse['t']=dec['t']
                mutse['p']=dec['p']
                mutse['!']=dec['!']
                mutse['raw']=m
                #print(mutse)
                haplo_muts_list.append(mutse)
            #print(muts)
        
        #subtree is at the end of list item
        ul = li.find('ul', recursive=False)
        if ul:
            #print('->')
            yfull_recurse_list(ul, level+1, False)
            #print('<-')
        else:
            if 'g' in muts and muts['g'].endswith('*'):
                continue
            if 'link' in muts:
                group=muts['link'].split('/')[-2]
                fname='yfull/yfull-mtree-'+group+'.html'
                #print('FILE: ' +fname)
                yfull_recurse_file(group, level)
                #print('END: ' +fname)
    return 0
    
def yfull_recurse_file(group, level):
    fname = 'yfull/yfull-mtree-'+group+'.html'
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

def import_yfull_snp():
    #TODO complex muts
    yfull_recurse_file('L', 0)
