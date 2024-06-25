#!/usr/local/bin/python3.6
import os, re, sys, gzip, pysam
import numpy as np
from datetime import datetime as dt_


def scSNV_assignment(score_tuple, a_thr, s_thr):
    '''
    a_thr: allele read threshold, min count of allele required
    s_thr: strand read threshold, min count of single strand required
    '''
    F_ref, F_alt, R_ref, R_alt = score_tuple
    F_sum, R_sum = F_ref + F_alt, R_ref + R_alt
    assignment = "*"
    if (F_sum+R_sum>=a_thr) and (min(F_sum,R_sum)>=s_thr):
        if F_ref==0 and R_ref==0:                 assignment="ds_Mut"
        elif F_ref*R_ref*F_alt*R_alt ==0:
            if min(R_ref+F_alt, F_ref+R_alt)==0:
                assignment="Dmg_Tran"
                if max(score_tuple)>=a_thr:       assignment="Dmg____0"
            elif (R_alt == F_alt == 0):           assignment="Ref_Call"
            elif (R_alt * F_alt  == 0):           assignment="Dmg_Miss"
            else:
                assignment="Dmg_Swap"
                if min(F_alt,R_alt)>=s_thr:
                    if max(F_ref, R_ref)>=s_thr:  assignment="Dmg_Swap_sf" # strand filtered
        else:                                     assignment="Swap&Mis"
    elif (F_sum+R_sum>=a_thr) and (min(F_sum,R_sum)==0) and (F_ref==R_ref==0):
        assignment="ss_Var"

    return assignment

### def remove_bad_duplexes(read_set, var__pos, FILE__err):
###     sortread = {}
###     for read_pair in read_set:
###         for read_it in read_set[read_pair].values():
###             ### get read_pair positions
###             # blocks are genomic positions, forward: [0][0] ; reverse: [-1][-1]
###             # isize signed by orientations, forward: ++     ; reverse: --
###             forward = read_it.is_forward
###             stt_pos = read_it.blocks[forward-1][forward-1] +forward  
###             end_pos = stt_pos + read_it.isize +1 -2*forward
###             min_cor, max_cor = min(stt_pos,end_pos), max(stt_pos,end_pos)
###             
###             ### get UMI for molecules
###             tags = dict(read_it.tags).get("MR", "")
###             
###             ### assign read pairs to duplexes
###             # [[0,0]] is a legacy for an early try to allow coordinate
###             # miss-match within 1-3 bps, deprecated since it's inaccurate
###             new = True
###             for i,j in [[0,0]]:
###                 real_UMI = f"{min_cor+i}_{max_cor+j}{tags}"
###                 if real_UMI in sortread:
###                     sortread[real_UMI].append(read_set[read_pair]) ; new = False
###             if new:
###                 sortread[f"{min_cor}_{max_cor}_{tags}"] = [read_set[read_pair]]
###             break # to avoid counting a single strand twice
###     dp_2_del = []
###     for duplex_i in sortread:
###         if len(sortread[duplex_i])<4: dp_2_del.append(duplex_i)
###     for duplex_i in dp_2_del: sortread.pop(duplex_i)
###     return sortread

def remove_bad_duplexes(read_set, var__pos, FILE__err):
    sortread = {}
    for read_pair in read_set:
        for read_it in read_set[read_pair].values():
            ### get read_pair positions
            # blocks are genomic positions, forward: [0][0] ; reverse: [-1][-1]
            # isize signed by orientations, forward: ++     ; reverse: --
            forward = read_it.is_forward
            stt_pos = read_it.blocks[forward-1][forward-1] +forward
            end_pos = stt_pos + read_it.isize +1 -2*forward
            min_cor, max_cor = min(stt_pos,end_pos), max(stt_pos,end_pos)
            
            ### get UMI for molecules
            tags = dict(read_it.tags).get("MR", "coor")
            
            ### assign read pairs to duplexes
            real_UMI = f"{min_cor}_{max_cor}_{tags}"
            sortread.setdefault(real_UMI, []).append(read_set[read_pair])
            break # to avoid counting a single strand twice
    dp_2_del = []
    for duplex_i in sortread:
        if len(sortread[duplex_i])<4: dp_2_del.append(duplex_i)
    for duplex_i in dp_2_del: del sortread[duplex_i]
    return sortread

def Count_Duplex(in_dplxs, var__pos, var__alt, FILE__err):
    duplex__cnt, duplex__pos, duplex_stat = {}, {}, {}
    for_flags = {99, 147, 99+1024, 147+1024}
    for duplex_i in in_dplxs:
        duplex__cnt.setdefault(duplex_i, [0,0,0,0])
        q30_counts=[] ; Err_mALT, Err_diss, Err_lowQ = 0,0,0
        for read_pair in in_dplxs[duplex_i]:
            # get the concensus between 2 reads
            cons__q30, cons_call = [], [] 
            for read_obj in read_pair.values():
                read_vloc = read_obj.get_reference_positions().index(var__pos-1)
                ampl_fowd = read_obj.flag in for_flags
                read_qual = read_obj.query_alignment_qualities[read_vloc]  # base quality
                read_base = read_obj.query_alignment_sequence[ read_vloc]  # mutated base
                
                cons__q30.append(int(read_qual>=30))
                cons_call.append(ampl_fowd*2+int(read_base==var__alt) if read_base in {var__ref, var__alt} else -1)
            
            # when the two agree with each other
            if len(set(cons_call)) == 1:
                if cons_call[0] != -1:
                    q30_counts.append(max(cons__q30))
                    if max(cons__q30)==1: 
                        duplex__cnt[duplex_i][cons_call[0]] += 1
                else:   q30_counts.append(0) ; Err_mALT += 1
            # if the two are different, trust the one with higher Q30 (or none)
            elif len(cons_call) == 2:
                if   sum(cons__q30)  == 2:  q30_counts.append(0) ; Err_diss += 1 # two read disagree 
                elif sum(cons__q30)  == 0:  q30_counts.append(0) ; Err_lowQ += 1 # both reads low Q30
                elif     cons__q30[0]== 1: 
                    if   cons_call[0]!=-1:  q30_counts.append(1) ; duplex__cnt[duplex_i][cons_call[0]] += 1
                    else:                   q30_counts.append(0) ; Err_mALT += 1
                elif     cons__q30[1]== 1: 
                    if   cons_call[1]!=-1:  q30_counts.append(1) ; duplex__cnt[duplex_i][cons_call[1]] += 1
                    else:                   q30_counts.append(0) ; Err_mALT += 1
        
        # finished counting, conduct final screening for base Q30
        if   duplex__cnt[duplex_i][1]+duplex__cnt[duplex_i][3]<2:
            del duplex__cnt[duplex_i] ; continue
        elif q30_counts.count(1)/len(q30_counts)<0.6:
            err_line = "\t".join( vcf_line + [duplex_i, str(duplex__cnt[duplex_i]), "low_q30"]) + "\n"
            FILE__err.write(err_line) ; del duplex__cnt[duplex_i] ; continue
        else:
            duplex__pos[duplex_i] = [abs(var__pos-int(i))+1 for i in duplex_i.split("_")[:2]]
            duplex_stat[duplex_i] = [Err_mALT, Err_diss, Err_lowQ, np.abs(read_obj.isize), 
                                     dict(read_obj.tags).get("AS"), dict(read_obj.tags).get("XS", 0)]
    return duplex__cnt, duplex__pos, duplex_stat


def report_progess(vcf_iter, var__chr, start_tme):
    vcf_iter += 1
    if vcf_iter%2000==0:
        avg_time = (dt_.now()-start_tme)/vcf_iter*1000
        print(f"""\rrun iter {vcf_iter}, now at {var__chr}, time per 1k variant: {
                    avg_time.seconds}.{avg_time.microseconds//1000:0>3} sec, projected remaining time: {
                    str(avg_time/1000*(vcf_totl-vcf_iter))[:-7]}, projected total time:{
                    str(dt_.now()-start_tme+avg_time/1000*(vcf_totl-vcf_iter))[:-7]}{' '*20}""", end="")


#################################################################################################
print("config env")
#################################################################################################
# parse arguments
PATH_work =     sys.argv[1]
PATH__vcf =     sys.argv[2]
PATH__bam =     sys.argv[3]
ale_thres = int(sys.argv[4])
std_thres = int(sys.argv[5])
PATH__out =     sys.argv[6]
read__len = int(sys.argv[7])
PATH__vis = PATH__out.replace(".txt","")+".vis"
PATH__err = PATH__out.replace(".txt","")+".err"
try:    print_process = bool(sys.argv[8])
except: print_process = False
os.chdir(PATH_work)
chr_list = [f"chr{i}" for i in list(range(1,23))+["X","Y"]]
flag_set = {99, 147, 83, 163,99+1024,147+1024,83+1024,163+1024}

#################################################################################################
### get total mutation to iterate
#################################################################################################

header_cnt = 0
FILE__vcf = gzip.open(PATH__vcf, mode="rt") if PATH__vcf[-3:]==".gz" else open(PATH__vcf,"r")
for vcf_totl, vcf_line in enumerate(FILE__vcf):
    if vcf_line[0] == "#": header_cnt += 1
FILE__vcf.close()

print(f"""
working in dir{PATH_work}
vcf file: {PATH__vcf} with {vcf_totl} variants
bam file: {PATH__bam}
coverage threshold: allele: {ale_thres} ; strand: {std_thres}
result saved in: {PATH__out}
""")

#################################################################################################
print("iterating")
#################################################################################################

FILE__vcf = gzip.open(PATH__vcf, mode="rt")
FILE__bam = pysam.AlignmentFile(PATH__bam, "rb")
FILE__out = open(PATH__out,"w")
FILE__vis = open(PATH__vis,"w")
FILE__err = open(PATH__err,"w")

print("start")
var__chr = "chr0" ; start_tme = dt_.now()
for vcf_iter, vcf_line in enumerate(FILE__vcf):
    if vcf_line[0]=="#": continue
    report_progess(vcf_iter, var__chr, start_tme)
    vcf_line = vcf_line.strip().split("\t")
    var__chr, var__pos, var__ref, var__alt = vcf_line[0], int(vcf_line[1]), vcf_line[3], vcf_line[4]
    if var__chr not in chr_list:             FILE__err.write("\t".join(vcf_line + ["unknown_chr"])+"\n"); continue
    if not len(var__ref)==len(var__alt)==1:  FILE__err.write("\t".join(vcf_line + ["not__a__SBS"])+"\n"); continue
    
    ### basic QC for reads
    cov_read = {}
    for i, read_obj in enumerate(FILE__bam.fetch(var__chr, var__pos-1, var__pos)):
        #if np.abs(read_obj.isize)<200:                                                     continue # read pair too short
        if read_obj.flag not in flag_set:                                                  continue # un-proper pairs
        if len(re.findall("[a-z]", read_obj.get_reference_sequence())) > 3:                continue # too many softmasked
        if len(re.findall("[ID]",  read_obj.cigarstring)) > 0:                             continue # contain indels
        if (read_obj.cigar[ 0][0] == 4) & (read_obj.cigar[ 0][1] > 20):                    continue # softclipped
        if (read_obj.cigar[-1][0] == 4) & (read_obj.cigar[-1][1] > 20):                    continue # softclipped
        #try:             qual_dlta = dict(read_obj.tags)["AS"] - dict(read_obj.tags)["XS"]
        #except KeyError: qual_dlta = float("inf")
        #if qual_dlta <=90:                                                                 continue
        cov_read.setdefault(read_obj.qname, {})[read_obj.flag]=read_obj
    if len(cov_read)<4: continue
    
    good_duplex = remove_bad_duplexes(cov_read, var__pos, FILE__err)                      ### reconstruct duplexes
    cnt_results = Count_Duplex(good_duplex, var__pos, var__alt, FILE__err)
    duplex__cnt, duplex__pos, duplex_stat = cnt_results                                   ### count variants
    for duplex_i in duplex__cnt:                                                          ### sort out final calls and output result
        var_call = scSNV_assignment(duplex__cnt[duplex_i], ale_thres, std_thres)
        if min (duplex__pos[duplex_i][0],duplex__pos[duplex_i][1]) > read__len - 8: continue
        aMsNs = ",".join([str(i) for i in duplex__cnt[duplex_i]])
        pos_r = ",".join([str(i) for i in duplex__pos[duplex_i]])
        dup_s = ",".join([str(i) for i in duplex_stat[duplex_i]])
        
        out_line = "\t".join([var_call] + vcf_line     + [          aMsNs, pos_r]) + "\n"
        vis_line = "\t".join(             vcf_line[:5] + [duplex_i, aMsNs, pos_r, var_call, dup_s]) + "\n"
        FILE__out.write(out_line)
        FILE__vis.write(vis_line)