#!/usr/bin/python3
import re
import csv
import os
from bs4 import BeautifulSoup
import urllib.request
import json
import glob
import csv
import argparse


#all-ancient-dna-2-07-06.csv (CC-BY indo-european.eu Carlos Quiles, Jean Manco)
#https://indo-european.eu/ancient-dna/
#csv: https://drive.google.com/drive/folders/1FNQNLQs93tmsX5_G728zE4DTIS0WUsXR
#iconv -f windows-1252 -t utf-8 -c all-ancient-dna-2-07-73.csv >all-ancient-dna-2-07-73b.csv
def import_ancient():
    fname='all-ancient-dna-2-07-73b.csv'
    #annos['info']=fname+' (CC-BY indo-european.eu Carlos Quiles, Jean Manco)'
    with open(fname) as f:
        csv_reader=csv.reader(f, delimiter=',')
        lc = 0
        for row in csv_reader:
            #print(row)
            if lc == 0:
                lc += 1
            else:
                id=row[0]
                mt=row[11].strip('/').split('/')[-1]
                y=row[33].strip('/').split('/')[-1].strip('*')
                fty=row[35].strip('/').split('/')[-1]
                src=row[44]
                date=row[46].strip()
                bc1=row[48]
                bc2=row[49]
                age=row[50]
                simp_cult=row[51].strip()
                cult=row[52].strip()
                loc=row[54].strip()
                cou=row[56].strip()
                #print(id, mt, y, src, date, bc1, bc2, age, simp_cult, cult, loc, cou)
                txt='ANCIENT %s (%s): %s, %s %s, %s %s'%(id, src, date, simp_cult, cult, loc, cou)
                print(txt)
                anno = {
                    "g": y,
                    "txt": txt
                }
                annos.append(anno)


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
    if group_name=='R-M335':
        return True
    if group_name=='R-U106':
        return True
    if group_name=='K-Y28299':
        return True
    if group_name=='H':
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
            age=yfull_parse_age(li)
            if age:
                print(grp, age)
                anno = {
                    "g": grp,
                    "txt": 'YF-AGE: %s'%(age)
                }
                annos.append(anno)
            for txt in txts:
                print(grp, txt)
                anno = {
                    "g": grp,
                    "txt": 'YF %s'%(txt)
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





def ftdna_check_exists(kit):
    if kit in anno_by_kit:
        return True
    return False
#
def import_ftdna_chart(fname, info=''):
    with open(fname, encoding="UTF-8") as f:
        print('Importing file: ' +fname)
        soup = BeautifulSoup(f.read(), features="html.parser")
        
        rows = soup.find('div', {"id" : re.compile('MainContent.*')}).find('table').find_all("tr")
        
        kiti = -1
        pati = -1
        coui = -1
        gri = -1
        row = rows[0]
        ths = row.find_all("th")
        heading_mut = ''
        heading_mut2 = ''
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
            if len(tds)==1:
                hding = tds[0].get_text().strip()
                muts = re.findall('>\s*[a-zA-Z]{1,3}[0-9]{1,7}', hding)
                if muts:
                    heading_mut = muts[-1][1:]
                    heading_mut = muts[-1].strip('>').strip()
                else:
                    heading_mut = ''
                muts2 = re.findall('\s*[a-zA-Z]{1,3}[0-9]{1,7}\+', hding)
                if muts2:
                    heading_mut2 = muts2[-1].strip('+').strip()
                else:
                    heading_mut2 = ''
                print('HEADING:', hding, muts, heading_mut, muts2, heading_mut2)
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
                grclass = tds[gri].find('span')['class'][0]
                gr = tds[gri].get_text().strip()
                if not gr:
                    continue
                print(gr, grclass, heading_mut, heading_mut2, kit, pat)
                anno = {}
                if grclass == 'haplo1':
                    if heading_mut:
                        anno['m'] = heading_mut
                    elif heading_mut2:
                        anno['m'] = heading_mut2
                    else:
                        anno['g'] = gr
                else:
                    anno['g'] = gr
                anno['txt'] = 'FT %s %s %s'%(kit, pat, info)
                if ftdna_check_exists(kit):
                    #TODO: keep best match?
                    print('EXISTS:', anno)
                    continue
                print(anno)
                annos.append(anno)
                if 'g' in anno:
                    anno_by_kit[kit]=anno['g']
                if 'm' in anno:
                    anno_by_kit[kit]=anno['m']

def save_anno(fname):
    jroot={
        'info': 'haploy_anno_import.py',
        'annotation': annos }
    with open(fname, 'w') as f:
        json.dump(jroot, f, indent=1);

def init_annos():
    global annos
    global anno_by_kit
    annos=[]
    anno_by_kit={}


def import_single_ft_project(fname):
    print(fname)
    tag = ''
    tag = fname.split(' - ')[-1]
    tag = tag.split('.')[0]
    tag = tag.replace(' ', '')
    ofname = 'haploy_annodb_' + tag.lower() + '.txt'
    tag = '[' + tag + ']'
    print(tag)

    init_annos()
    import_ftdna_chart(fname, tag)
    save_anno(ofname)

def import_example():
    # Example annotations - it probably doesn't make sense for everyone to import every project


    #annos=[]
    init_annos()
    import_ancient()
    save_anno('haploy_annodb_ancientdna.txt')

    #annos=[]
    init_annos()
    import_yfull_tree('A00')
    import_yfull_tree('A0-T')
    #import_yfull_tree('N-FGC28435')
    #import_yfull_tree('N')
    save_anno('haploy_annodb_yfull.txt')

    #annos=[]
    init_annos()
    #Go to project DNA Results->Classic Chart (e.g. https://www.familytreedna.com/public/Finland?iframe=yresults), set Page Size to a big number, load the new page and Save as
    import_ftdna_chart('ftdna/FamilyTreeDNA - Viitasaari (Savo-Tawastian) DNA project.htm', '[Viitasaari]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - Turun Seudun Sukututkijat Dna.htm', '[Turun Seudun]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - Det skogfinske DNA-prosjektet - Forest Finn DNA project.htm', '[Skogfinske]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - Saami Project.htm', '[Saami]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - Karjala DNA -projekti.htm', '[Karjala]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - Finland DNA Project.htm', '[FinlandDNA]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - Sweden DNA PROJECT - Sverigeprojektet.htm', '[SwedenDNA]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - The Norway DNA Project - Norgesprosjektet.htm', '[NorwayDNA]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - Estonia.htm', '[Estonia]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - RussiaDNA Project.htm', '[RussiaDNA]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - Lituania Propria.htm', '[Lituania]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - My FamilyTree DNA Latvia Project Website.htm', '[LatviaDNA]')
    import_ftdna_chart('ftdna/FamilyTreeDNA -.htm', '[BalticSea]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - I1 Suomi Finland & N-CTS8565 -projekti.htm', '[I1 Suomi]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - R1a1a and Subclades Y-DNA Project.htm', '[R1a1a]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - R1b and Subclades Project for R1b (M343+ and M269+) Y DNA Haplogroup.htm', '[R1b]')
    import_ftdna_chart('ftdna/FamilyTreeDNA - N1c1 Haplogroup Y-DNA Project.htm', '[N1c1]')
    save_anno('haploy_annodb_ftdnatest.txt')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #parser.add_argument('file', nargs='+')
    #parser.add_argument('file')
    parser.add_argument('-f', '--file', help='Output build')
    args = parser.parse_args()

    #if len(args.file) > 0:
    if args.file:
        #for fname in args.file:
        import_single_ft_project(args.file)
    else:
        import_example() 