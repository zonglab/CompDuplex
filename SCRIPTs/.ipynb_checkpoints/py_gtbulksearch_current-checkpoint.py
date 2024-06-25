#! /usr/local/bin/python2.7

import os, sys, pysam
from datetime import datetime as dt_

def report_progess():
    if (vcf_iter + 1) % r_span == 0:
        elapsed_time = dt_.now() - rt
        pct_finished = (vcf_iter + 1) / total_scope * 100
        est_total_ti = elapsed_time / (vcf_iter + 1) * total_scope
        time__remain = est_total_ti - elapsed_time
        print(f"""({pct_finished:5.2f}%) parsed {vcf_iter + 1:7d} mutations, {
                    abscent:_>6,} abscent, {present:_>6,} present, {complex:_>5,} complex. {''
                   }runtime: {str(elapsed_time)[:-7]}, remaining time: {str(time__remain)[:-7]}.{' ' * 20}""", end="")

PATH_in_vcf =     sys.argv[1]
PATH_in_bam =     sys.argv[2]
PATH_ot_pre =     sys.argv[3]
max_de_novo = int(sys.argv[4])
PATH_o_soma = PATH_ot_pre + "_somatc.bdy"  # non--detected
PATH_o_lowF = PATH_ot_pre + "_low_fq.bdy"  # low frequency
PATH_o_comp = PATH_ot_pre + "_cp_err.bdy"  # complex errer
PATH_o_undr = PATH_ot_pre + "_undcov.bdy"  # under covered
print(PATH_o_soma, PATH_o_lowF, PATH_o_comp, PATH_o_undr)

total_scope = 0
with open(PATH_in_vcf, 'r') as FILE_in_vcf:
    for i, line in enumerate(FILE_in_vcf):
        if not line.startswith("#"): total_scope += 1
print(f"{total_scope} mutations to be parsed")

rt = dt_.now() ; r_span = 1000 ; depth_lw = 10 ; depth_up = 400
abscent, present, complex, undrcov = 0,0,0,0

FILE_in_bam = pysam.AlignmentFile(PATH_in_bam, "rb")
FILE_in_vcf = open(PATH_in_vcf, 'r')
FILE_o_soma = open(PATH_o_soma, "w")
FILE_o_lowF = open(PATH_o_lowF, "w")
FILE_o_comp = open(PATH_o_comp, "w")
FILE_o_undr = open(PATH_o_undr, "w")

for vcf_iter, line in enumerate(FILE_in_vcf):
    if not line.startswith("#"):
        report_progess()
        ### retrieve mutation information
        rows = line.strip().split("\t")
        chrnum, var__pos, var__ref, var__alt = rows[0], int(rows[1]), rows[3], rows[4]
        if (len(var__ref) > 1)|(len(var__alt) > 1): continue
        
        ### mutation stat in bulk bam
        denovo, readnum, readlist = 0, 0, set()
        for j, read_obj in enumerate(FILE_in_bam.fetch(chrnum , var__pos-1, var__pos)):
            ### for over-complicated regions, where there're too many reads
            if readnum > depth_up:
                complex +=1 ; FILE_o_comp.write(line) ; break
            ### combine read pairs
            if read_obj.query_name in readlist: continue
            else: readlist.add(read_obj.query_name)
            
            ### summarize
            if var__pos-1 in read_obj.get_reference_positions():
                read_vloc = read_obj.get_reference_positions().index(var__pos-1)
                if read_obj.query_alignment_qualities[read_vloc] > 0:
                    read_base = read_obj.query_alignment_sequence[ read_vloc]  # mutated base
                    readnum += 1
                    if read_base == var__alt: denovo += 1
        wlne = line.strip()+f"\t{readnum}\t{denovo}\t"
        if (depth_lw <= readnum) & (readnum <= depth_up):
            if (denovo <= max_de_novo):   abscent += 1 ; FILE_o_soma.write(wlne + "somatic\n")
            else:                         present += 1 ; FILE_o_lowF.write(wlne + "lowfreq\n")
        else:                             undrcov += 1 ; FILE_o_undr.write(wlne + "undrcov\n")




FILE_in_bam.close()
FILE_in_vcf.close()
FILE_o_soma.close()
FILE_o_lowF.close()
FILE_o_comp.close()