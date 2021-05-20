#!/usr/bin/python3
import re
import csv
import os
from bs4 import BeautifulSoup
import urllib.request
import json
import glob




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
    #urllib.request.urlretrieve("https://www.yfull.com/tree/"+group+"/", fname);

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

def yfull_parse_person(li):
    sams=[]
    ul = li.find('ul', recursive=False)
    if ul:
        lis = ul.find_all('li', recursive=False)
    else:
        return sams
    if not lis:
        return sams
    for li in lis:
        has_sample=0
        sam=''
        if li.has_attr('valsampleid'):
            sam+=li['valsampleid']+ ': '
            has_sample=1
        for geo in li.find_all('b', recursive=False):
            if geo.has_attr('class') and 'yf-geo' in geo['class'] and 'fl' in geo['class']:
                if geo.has_attr('title'):
                    sam+=geo['title']
                if geo.has_attr('original-title'):
                    sam+=geo['original-title']
                sam+=' '
            if geo.has_attr('class') and 'yf-geo' in geo['class'] and 'yf-lang' in geo['class']:
                if geo.has_attr('title'):
                    sam+=geo['title']
                if geo.has_attr('original-title'):
                    sam+=geo['original-title']
        for geo in li.find_all('span', recursive=False):
            if geo.has_attr('class') and 'yf-a-age' in geo['class']:
                if geo.has_attr('title'):
                    sam+=geo['title']
                if geo.has_attr('original-title'):
                    sam+=geo['original-title']
        if has_sample:
            sams.append(sam)
    #sam+=' '
    #print(sam)
    return sams

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
        #print(li.get_text())
        muts={}
        muts['l']=level
        g=li.find('a', recursive=False)
        group_name=''
        if g:
            group_name=g.text
            muts['g']=g.text
            txts = yfull_parse_person(li)
            grp=g.text.strip('*')
            for txt in txts:
                print(grp, txt)
                anno = {
                    "g": grp,
                    "txt": 'YFULL: %s'%(txt)
                }
                annos.append(anno)
        l=li.find('a', href=True, recursive=False)
        if l:
            muts['link']=l['href']


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
        #yfull_get_info(soup)

def import_yfull_tree(gr):
    yfull_recurse_file(gr, 0)





#
def import_ftdna_chart(fname, info=''):
    with open(fname) as f:
        print('Importing file: ' +fname)
        soup = BeautifulSoup(f.read(), features="html.parser")
        
        #rows = soup.find('div', id='MainContent_color1_GridView1').find('table').find_all("tr")
        #rows = soup.find('table').find_all("tr")
        rows = soup.find('div', {"id" : re.compile('MainContent.*')}).find('table').find_all("tr")
        
        kiti = -1
        pati = -1
        coui = -1
        gri = -1
        row = rows[0]
        ths = row.find_all("th")
        for i, th in enumerate(ths):
            if 'Kit' in th.get_text():
                kiti = i
            if 'Paternal' in th.get_text():
                pati = i
            if 'Country' in th.get_text():
                coui = i
            if 'Haplogroup' in th.get_text():
                gri = i
        for row in rows:
            tds = row.find_all("td")
            if len(tds)>1:
                kit=''
                pat=''
                cou=''
                gr=''
                kit = tds[kiti].get_text().strip()
                if pati >= 0:
                    pat = tds[pati].get_text().strip()
                if coui >= 0:
                    cou = tds[coui].get_text().strip()
                gr = tds[gri].get_text().strip()
                if not gr:
                    continue
                #if 'MIN' in kit or 'MAX' in kit or 'MODE' in kit:
                #    continue
                print(kit, pat, gr)
                anno = {
                    "g": gr,
                    "txt": 'FTDNA: %s %s %s'%(kit, pat, info)
                }
                annos.append(anno)

def save_anno(fname):
    jroot={
        'info': 'haploy_anno_iport.py',
        'annotation': annos }
    with open(fname, 'w') as f:
        json.dump(jroot, f, indent=1);


# Example annotations - it probably doesn't make sense for everyone to import every project

annos=[]
import_yfull_tree('A00')
import_yfull_tree('A0-T')
#import_yfull_tree('N-FGC28435')
#import_yfull_tree('N')
save_anno('haploy_annodb_yfull.txt')

annos=[]
import_ftdna_chart('ftdna/FamilyTreeDNA - Estonia.htm', '[Estonia]')
import_ftdna_chart('ftdna/FamilyTreeDNA - Saami Project.htm', '[Saami]')
import_ftdna_chart('ftdna/FamilyTreeDNA - I1 Suomi Finland & N-CTS8565 -projekti.htm', '[I1 Suomi]')
import_ftdna_chart('ftdna/FamilyTreeDNA - Finland DNA Project.htm', '[FinlandDNA]')
import_ftdna_chart('ftdna/FamilyTreeDNA - RussiaDNA Project.htm', '[RussiaDNA]')
import_ftdna_chart('ftdna/FamilyTreeDNA - R1a1a and Subclades Y-DNA Project.htm', '[R1a1a]')
save_anno('haploy_annodb_ftdnatest.txt')
