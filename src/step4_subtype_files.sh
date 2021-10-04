#!/bin/sh

# Generate per-subtype singleton files

awk -F, '{if($9 == "AT>CG")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > AT_CG.txt
awk -F, '{if($9 == "AT>GC")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > AT_GC.txt
awk -F, '{if($9 == "AT>TA")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > AT_TA.txt
awk -F, '{if($9 == "GC>AT")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > GC_AT.txt
awk -F, '{if($9 == "GC>TA")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > GC_TA.txt
awk -F, '{if($9 == "GC>CG")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > GC_CG.txt
awk -F, '{if($9 == "cpg_GC>AT")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > cpg_GC_AT.txt
awk -F, '{if($9 == "cpg_GC>TA")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > cpg_GC_TA.txt
awk -F, '{if($9 == "cpg_GC>CG")print(substr($10,1,21))}' chr*_filtered_annotated_addl.csv > cpg_GC_CG.txt

# Generate per-subtype control files

awk -F, '{if($4 == "AT>CG")print($3)}' chr*_rc.csv > AT_CG.txt
awk -F, '{if($4 == "AT>GC")print($3)}' chr*_rc.csv > AT_GC.txt
awk -F, '{if($4 == "AT>TA")print($3)}' chr*_rc.csv > AT_TA.txt
awk -F, '{if($4 == "GC>AT")print($3)}' chr*_rc.csv > GC_AT.txt
awk -F, '{if($4 == "GC>TA")print($3)}' chr*_rc.csv > GC_TA.txt
awk -F, '{if($4 == "GC>CG")print($3)}' chr*_rc.csv > GC_CG.txt
awk -F, '{if($4 == "cpg_GC>AT")print($3)}' chr*_rc.csv > cpg_GC_AT.txt
awk -F, '{if($4 == "cpg_GC>TA")print($3)}' chr*_rc.csv > cpg_GC_TA.txt
awk -F, '{if($4 == "cpg_GC>CG")print($3)}' chr*_rc.csv > cpg_GC_CG.txt