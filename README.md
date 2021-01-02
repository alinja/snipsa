# SNIPSA

A small project experimenting with SNP genome data and python.

## haploy_find.py

This small tool reads a raw SNP data file and lists Y chromosome haplogroup information. 

```
./haploy_find.py <file>
```

The tool accepts also multiple files, and can be used to search haplogroups by the common haplogroup name (separated by comma) or mutation name.

```
./haploy_find.py N1a1,I2,L287 <files>...
```

## haploy_db_import.py

This script imports the ISOGG database to the internal format that is used by haplay_find.py. It needs CrossMap and the ISOGG spreadsheet in csv format.

