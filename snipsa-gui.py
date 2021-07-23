#!/usr/bin/python3
import sys
import os
import tkinter as tk
import tkinter.filedialog
from tkinter import scrolledtext 
import haplomt
import haploy
import haploy_anno_import

bam_support=1

if bam_support:
    import bamload
    import snpload

def handle_ftdna_select():
    f = tkinter.filedialog.askopenfile()
    if f:
        fname = f.name
    else:
        return
    haploy_anno_import.import_single_ft_project(fname)


def handle_bam_select():
    f = tkinter.filedialog.askopenfile()
    if f:
        bamname = f.name
    else:
        return
    in_build = bamload.get_build(bamname)
    bamload.setup_conv(in_build)
    snpset = bamload.full_convert(bamname)
    snpload.save(bamname+'.snp.txt', snpset, 38)
    snpname = bamname+'.snp.txt'
    message = 'SNP file written to ' + snpname
    scr.config(state=tk.NORMAL)
    scr.delete('1.0', tk.END)
    scr.insert('1.0', message)
    scr.config(state=tk.DISABLED)
    scr.see("end")
    fnamevar.set(snpname)


def handle_file_select():
    f = tkinter.filedialog.askopenfile()
    if f:
        fname = f.name
    else:
        return
    fnamevar.set(fname)

dbmt_loaded = 0
dby_loaded = 0

def handle_findmt():
    global dbmt_loaded
    if not dbmt_loaded:
        dbmt_loaded = 1
        haplomt.load_db()
    do_uptree = report_snps.get()
    do_all = report_all.get()
    force = pathvar.get()
    vcf_sample = vcfvar.get().strip()
    buildstr = buildvar.get()
    force_build=0
    if  buildstr == 'Build36': force_build = 36
    if  buildstr == 'Build37': force_build = 37
    if  buildstr == 'Build38': force_build = 38
    fname = fnamevar.get()
    print("Reporting file: "+fname)
    num = int(numbox.get())
    rep = haplomt.report(fname, num, do_uptree=do_uptree, do_extra=do_uptree, do_all=do_all, filt=force, force=force, vcf_sample=vcf_sample, force_build=force_build)
    print("Done")
    scr.config(state=tk.NORMAL)
    scr.delete('1.0', tk.END)
    scr.insert('1.0', rep)
    scr.config(state=tk.DISABLED)
    scr.see("end")
   
def handle_findy():
    global dby_loaded
    if not dby_loaded:
        dby_loaded = 1
        haploy.load_db2j()
        haploy.load_annotations('haploy_annodb_*.txt')
    do_uptree = report_snps.get()
    do_all = report_all.get()
    force = pathvar.get()
    vcf_sample = vcfvar.get().strip()
    buildstr = buildvar.get()
    force_build=0
    if  buildstr == 'Build36': force_build = 36
    if  buildstr == 'Build37': force_build = 37
    if  buildstr == 'Build38': force_build = 38
    fname = fnamevar.get()
    print("Reporting file: "+fname)
    num = int(numbox.get())
    rep = haploy.report(fname, num, do_uptree=do_uptree, do_extra=do_uptree, do_all=do_all, filt=force, force=force, vcf_sample=vcf_sample, force_build=force_build)
    print("Done")
    scr.config(state=tk.NORMAL)
    scr.delete('1.0', tk.END)
    scr.insert('1.0', rep)
    scr.config(state=tk.DISABLED)
    scr.see("end")

window = tk.Tk()
window.title("Snipsa GUI")
window.geometry("950x700")
window.minsize(600, 300)

hdrframe=tk.Frame(window)
hdrframe.pack(fill='x')

# File
iframe=tk.Frame(hdrframe)
iframe.pack(fill='both')
cbutton4 = tk.Button(iframe, text="Import FTDNA project html", command=handle_ftdna_select)
cbutton4.pack(side=tk.LEFT, fill='x', expand=True)
if bam_support:
    cbutton3 = tk.Button(iframe, text="Import BAM", command=handle_bam_select)
    cbutton3.pack(side=tk.LEFT, fill='x', expand=True)

# File
fframe=tk.Frame(hdrframe)
fframe.pack(fill='both')
button = tk.Button(fframe, text="Choose file", command=handle_file_select)
button.pack(side=tk.LEFT)
fnamevar = tk.StringVar()
fnamevar.set("No file selected")
fnamelabel=tk.Label(fframe, textvariable=fnamevar,  anchor='w')
fnamelabel.pack(side=tk.LEFT, fill='both', expand=True)
vcflabel=tk.Label(fframe, text='VCF sample:', anchor='w')
vcflabel.pack(side=tk.LEFT, fill='both')
vcfvar = tk.StringVar()
vcfbox = tk.Entry(fframe, textvariable=vcfvar,width=10)
vcfbox.pack(side=tk.LEFT)
buildvar = tk.StringVar()
buildchoices = {'Auto', 'Build36', 'Build37', 'Build38'}
buildvar.set('Auto')
builddropdown = tk.OptionMenu(fframe, buildvar, *buildchoices)
builddropdown.pack(side=tk.LEFT)

# Settings
sframe=tk.Frame(hdrframe)
sframe.pack(fill='both')
nbestlabel=tk.Label(sframe, text='Show n best matches:', anchor='w')
nbestlabel.pack(side=tk.LEFT, fill='both')
numbox = tk.Spinbox(sframe,from_=1, to=1000, width=3, text="aaaa")
numbox.pack(side=tk.LEFT)
report_snps = tk.IntVar(value=1)
snplabel=tk.Label(sframe, text='Show path:', anchor='w')
snplabel.pack(side=tk.LEFT, fill='both')
rscbox = tk.Checkbutton(sframe, variable=report_snps)
rscbox.pack(side=tk.LEFT)
report_all = tk.IntVar(value=0)
alllabel=tk.Label(sframe, text='Show all SNPs:', anchor='w')
alllabel.pack(side=tk.LEFT, fill='both')
allbox = tk.Checkbutton(sframe, variable=report_all)
allbox.pack(side=tk.LEFT)
pathlabel=tk.Label(sframe, text='Exact path to:', anchor='w')
pathlabel.pack(side=tk.LEFT, fill='both')
pathvar = tk.StringVar()
pathbox = tk.Entry(sframe, textvariable=pathvar,width=6)
pathbox.pack(side=tk.LEFT)

# Command
cframe=tk.Frame(hdrframe)
cframe.pack()
cbutton1 = tk.Button(cframe, text="Find MT", command=handle_findmt)
cbutton1.pack(side=tk.LEFT, fill='x', expand=True)
cbutton2 = tk.Button(cframe, text="Find Y", command=handle_findy)
cbutton2.pack(side=tk.LEFT, fill='x', expand=True)

# Result area
scr=scrolledtext.ScrolledText(window, wrap=tk.WORD)  
scr.config(state=tk.DISABLED)
scr.pack(fill='both', expand=True)

window.mainloop()