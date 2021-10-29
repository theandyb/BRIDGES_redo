#!/bin/sh

# Generate per-subtype singleton files

awk -F, '{if($9 == "AT_CG")print(substr($8,1,21))}' chr*_annotated.csv > AT_CG.txt
awk -F, '{if($9 == "AT_GC")print(substr($8,1,21))}' chr*_annotated.csv > AT_GC.txt
awk -F, '{if($9 == "AT_TA")print(substr($8,1,21))}' chr*_annotated.csv > AT_TA.txt
awk -F, '{if($9 == "GC_AT")print(substr($8,1,21))}' chr*_annotated.csv > GC_AT.txt
awk -F, '{if($9 == "GC_TA")print(substr($8,1,21))}' chr*_annotated.csv > GC_TA.txt
awk -F, '{if($9 == "GC_CG")print(substr($8,1,21))}' chr*_annotated.csv > GC_CG.txt
awk -F, '{if($9 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_annotated.csv > cpg_GC_AT.txt
awk -F, '{if($9 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_annotated.csv > cpg_GC_TA.txt
awk -F, '{if($9 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_annotated.csv> cpg_GC_CG.txt

# Generate per-subtype control files

awk -F, '{if($4 == "AT_CG")print($8)}' chr*_at.csv | sed 's/"//g' > AT_CG.txt
awk -F, '{if($4 == "AT_GC")print($8)}' chr*_at.csv | sed 's/"//g' > AT_GC.txt
awk -F, '{if($4 == "AT_TA")print($8)}' chr*_at.csv | sed 's/"//g' > AT_TA.txt
awk -F, '{if($4 == "GC_AT")print($8)}' chr*_gc.csv | sed 's/"//g' > GC_AT.txt
awk -F, '{if($4 == "GC_TA")print($8)}' chr*_gc.csv | sed 's/"//g' > GC_TA.txt
awk -F, '{if($4 == "GC_CG")print($8)}' chr*_gc.csv | sed 's/"//g' > GC_CG.txt
awk -F, '{if($4 == "cpg_GC_AT")print($8)}' chr*_gc.csv | sed 's/"//g' > cpg_GC_AT.txt
awk -F, '{if($4 == "cpg_GC_TA")print($8)}' chr*_gc.csv | sed 's/"//g' > cpg_GC_TA.txt
awk -F, '{if($4 == "cpg_GC_CG")print($8)}' chr*_gc.csv | sed 's/"//g' > cpg_GC_CG.txt