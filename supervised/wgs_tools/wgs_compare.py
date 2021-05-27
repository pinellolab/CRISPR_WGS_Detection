import glob
import pandas as pd
import numpy as np
import re
import os
import glob
from collections import defaultdict
import multiprocessing
from numba import jit

def get_all_seqs(desired_seq_len, prev_seqs=[], bases=['G', 'A', 'T', 'C']):
    if len(prev_seqs) > 0:
        if len(prev_seqs[0]) == desired_seq_len:
            return prev_seqs
    else:
        return get_all_seqs(desired_seq_len, prev_seqs=bases, bases=bases)
    new_seqs = []
    for prev_seq in prev_seqs:
        for base in bases:
            new_seqs.append(prev_seq + base)
    return get_all_seqs(desired_seq_len, prev_seqs=new_seqs, bases=bases)

def get_outcomes(left_bit, right_bit, max_insertions, max_deletions):
    
    # test_f = open("./test.csv", "w+")
    
    lines = []
    columns=['Category', 'Genotype position', 'Inserted Bases',
              'Length', 'Microhomology length', 'Predicted frequency', 'Genotype', 'Name', 'Aligned_sequence']
    # wt = {'Category':'WT', 'Genotype position'}
    lines.append(",".join(columns))
    
    aligned_seqs = []
    bases = ['G', 'A', 'T', 'C']
    # repair_outcomes = defaultdict(int)
    for del_len in range(1, max_deletions + 1):
        for del_left in range(0, del_len + 1):
            del_right = del_len - del_left
            # print(del_left, del_right)
            
            seq = left_bit[:len(left_bit) - del_left] + right_bit[del_right:]
            ref_seq =  left_bit[:len(left_bit) - del_left] + "-"*del_len + right_bit[del_right:]
            if len(ref_seq) > len(left_bit+right_bit):
                continue
            
            # repair_outcomes[seq] += 1
            # print("Deletion", del_len)
            if del_left == 0:
                del_lseq = ""
            else:
                del_lseq = left_bit[-del_left:]
            gt_dict = {"Category":"del", "Genotype position": del_left, "Inserted Bases": '', 'Microhomology length':'',
                       "Length": del_len, 'Predicted frequency':0.0,  "Genotype": seq, "Name": "del" + del_lseq + right_bit[:del_right], 'Aligned_sequence':ref_seq}
            #  if len(del_lseq + right_bit[:del_right]) != del_len:
                # print("WARNING!")
            
            # aligned_seqs.append(ref_seq)
            # genotype_df = genotype_df.append(gt_dict, ignore_index=True)
            # genotype_df.loc[genotype_df.shape[0]+1] = gt_dict
            line_i = [str(gt_dict[key]) for key in columns]
            lines.append(",".join(line_i))
            # test_f.write(",".join(line_i) + "\n")

    for ins_len in range(1, max_insertions + 1):
        ins_seqs = get_all_seqs(ins_len)
        for ins_seq in ins_seqs:
            seq = left_bit + '-' + right_bit
            ref_seq = left_bit + ins_seq + right_bit
            # print(seq)
            # print(ref_seq)
            # repair_outcomes[seq] += 1
            # print("Insertion", ins_seq)
            gt_dict = {"Category": "ins", "Genotype position": ' ', "Inserted Bases": ins_seq,
                       'Microhomology length': ' ',
                       "Length": len(ins_seq), 'Predicted frequency': 0.0, "Genotype": seq,
                       "Name": "ins" + ins_seq, 'Aligned_sequence':ref_seq}
            ## genotype_df.loc[genotype_df.shape[0]+1] = gt_dict
            ## aligned_seqs.append(ref_seq)
            line_i = [str(gt_dict[key]) for key in columns]
            lines.append(",".join(line_i))
            # test_f.write(",".join(line_i) + "\n")
    
    # test_f.close()
    # genotype_df['Aligned_sequence'] = aligned_seqs
    return lines

def build_aligned_sequence_for_indelphi(seq, cut_loc, df):
    prediction_df = df
    aligned_sequences = []
    for idx, row in prediction_df.iterrows():
        # gtype_seq = row['Genotype']
        aligned_seq = row['Genotype']
        if row['Category'] == 'del':
            del_seq = row['Name'][3:]
            del_len = int(row['Length'])
            left_end = int(cut_loc - row['Genotype position'])
            right_start = int(cut_loc + row['Genotype position'] + 1)
            left_aligned_del_seq = seq[left_end + 1:left_end + del_len + 1]
            right_aligned_del_seq = seq[right_start-del_len:right_start]
            if left_aligned_del_seq == del_seq:
                aligned_seq = seq[:left_end + 1] + "-" * del_len + seq[left_end + del_len + 1:]
            elif right_aligned_del_seq == del_seq:
                aligned_seq = seq[:right_start-del_len] + "-"*del_len + seq[right_start:]
            else:
                print('No matched del sequence found!')
                return prediction_df
        aligned_sequences.append(aligned_seq)
        # print(len(seq), len(aligned_seq))
    prediction_df['Aligned sequence'] = aligned_sequences
    return prediction_df

def build_pair_seq(read_df, cutsite = 40, window_size=3):
    # idx start from 0
    cutsite = cutsite-1
    target_site = read_df.iloc[0]['Reference_Sequence']
    # print(target_site)
    cutting_win_seq = []
    cutting_loc = []
    for idx, row in read_df.iterrows():
        aligned_seq = row['Aligned_Sequence']
        win_seq = []
        if row['n_deleted'] > 0:
            del_len = row['n_deleted']
            #  match_start = aligned_seq.find("-"*del_len)
            match_sites = re.finditer("-"*del_len, aligned_seq)
            flag = 0
            for ms in match_sites:
                span = ms.span()
                if span[0] <= cutsite+1 and span[1]-1 >= cutsite-1:
                    loc = str(span[0])+":"+str(span[1])
                    start = span[0]-window_size
                    end = span[1]+window_size
                    if start < 0:
                        start = 0
                    #  print(aligned_seq, aligned_seq[start:end])
                    win_seq = aligned_seq[start:end]
                    flag = 1    
        elif row['n_inserted'] == 1 :
            ins_len = row['n_inserted']
            #########
            ins_aligned_seq = aligned_seq[cutsite-window_size:cutsite+window_size+1]
            loc = str(cutsite) + ':' + str(cutsite + 1)
            win_seq = ins_aligned_seq
        else:
            win_seq = aligned_seq
            loc = "0:"+ str(len(aligned_seq))
        #  print(win_seq)
        if len(win_seq) == 0:
            print(aligned_seq)
            print(row)
        cutting_win_seq.append(win_seq)
        cutting_loc.append(loc)
    read_df['cutting_seq'] = cutting_win_seq
    read_df['indel_location'] = cutting_loc
    # read_df['norm_Reads'] = read_df['#Reads']/sum(read_df['#Reads'])
    return read_df

def map_indephi_with_possibleOut(possible_df, indelphi_df):
    n = 0
    for idx, row in possible_df.iterrows():
        #  genotype_seq = row['Genotype']
        #  match_result = prediction_df[prediction_df['Genotype'] == genotype_seq]
        aligned_seq = row['Aligned_sequence']
        match_result = indelphi_df[indelphi_df['Aligned sequence'] == aligned_seq]
        if len(match_result) > 0:
            freq = sum(list(match_result['Predicted frequency']))
            possible_df.loc[idx, 'Predicted frequency'] = freq
            n += 1
    if len(indelphi_df) != n:
        print("Map indelphi_df with Reads_df ERROR!")
    return possible_df


def mapping_read_with_possible(reads_fn, possible_fn, n_pad = 3, max_mis_n = 1):
    tn = re.split("/", reads_fn)[-1]
    # print(tn)
    reads_df = pd.read_csv(reads_fn)
    reads_df['match_flag'] = 0
    reads_df['match_idx'] = -1
    reads_df['match_len_possible'] = -1
    possible_df = pd.read_csv(possible_fn)
    # possible_df['norm_Reads_combined'] = 0
    possible_df['Reads_combined'] = 0
    possible_df['Matched_reads'] = 0
    possible_df['Mismatched_reads'] = 0
    
    for idx_p, row_p in possible_df.iterrows():
        indel_len = row_p['Length']
        if row_p['Category'] == 'del':
            pt = re.compile("-"*indel_len)
            match_res = list(re.finditer(pt, row_p['Aligned_sequence']))
            """
            if len(match_res) == 0:
                pt = re.compile("-"*indel_len+"[A-Z]")
                match_res = list(re.finditer(pt, row_p['Aligned_sequence']))
                if len(match_res) == 0:
                    pt = re.compile("[A-Z]"+"-"*indel_len)
                    match_res = list(re.finditer(pt, row_p['Aligned_sequence']))
            """
            # print(match_res, indel_len)     
            start_idx, end_idx = match_res[0].span()
            indel_loc_p = str(start_idx)+":"+str(end_idx)
            sl_reads_df = reads_df[reads_df['n_deleted'] == indel_len]
            
        elif row_p['Category'] == 'ins':
            cutsite = int(len(row_p['Aligned_sequence'])/2-1)
            start_idx, end_idx = cutsite, cutsite+1
            indel_loc_p = str(start_idx)+":"+str(end_idx)
            sl_reads_df = reads_df[reads_df['n_inserted'] == 1]
        else:
            continue
            
        for idx_r, row_r in sl_reads_df.iterrows():
            if indel_loc_p == row_r['indel_location']:
                aligned_seq_r = row_r['Aligned_Sequence']
                aligned_seq_p = row_p['Aligned_sequence']
                window_seq_r = aligned_seq_r[start_idx-n_pad:end_idx+n_pad]
                window_seq_p = aligned_seq_p[start_idx-n_pad:end_idx+n_pad]
                mis_n = 0
                for i in range(len(window_seq_p)):
                    if window_seq_r[i] != window_seq_p[i]:
                        mis_n += 1
                if mis_n <= 10:
                    possible_df.loc[idx_p, 'Reads_combined'] += row_r['#Reads']
                    possible_df.loc[idx_p, 'Mismatched_reads'] += row_r['#Reads']
                    reads_df.loc[idx_r, 'match_flag'] += 1
                    reads_df.loc[idx_r, 'match_len_possible'] = row_p['Length']
                    reads_df.loc[idx_r, 'match_idx'] = idx_p
                    # print(window_seq_r, window_seq_p)
                              
    return reads_df, possible_df
def mark_gentype_idx(df):
    # df = pd.read_csv("./CRISPRessoWGS_VEGFA3_win80_possible/offby3_332.csv")
    pr_genotype = df.iloc[0]['Genotype']
    geno_id = 1
    seqTypes = []
    for idx, row in df.iterrows():
        genotype = row['Genotype']
        if pr_genotype == genotype:
            seqTypes.append(geno_id)
        else:
            geno_id+=1
            pr_genotype = genotype
            seqTypes.append(geno_id)
    df['genotype_idx'] = seqTypes
    print(len(df), df.iloc[-1]['genotype_idx'])
    return df
