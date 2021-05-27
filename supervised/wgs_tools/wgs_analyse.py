import mapping_tools as mtools
import numpy as np
import glob
import os
from IPython.display import Image
import re
import tqdm
import mapping_tools as mtool
import multiprocessing
import pandas as pd
import json
import ast 
from scipy import stats
from numba import jit

class WGS_analyser:
    def __init__(self, possible_fn=None, reads_fn=None, indelphi_fn=None, ref_seq=None, reads_type='Treated'):
        self.reads_type = reads_type
        self.possible_fn = possible_fn
        self.reads_fn = reads_fn
        self.indelphi_fn = indelphi_fn
        self.ref_seq = ref_seq
        
    def save_df(self, df, save_path, sep_mark=","):
        df.to_csv(save_path, index=False, sep=sep_mark)
    
    def build_aligned_seq_for_inDelphi(self):
        df = pd.read_csv(self.indelphi_fn)
        # refseq_df = pd.read_csv(self.targetSeq_fn, index_col='target_name')
        # tname = re.split("/|.csv", self.indelphi_fn)[-2]
        # ref_seq = refseq_df.loc[tname, 'target_seq']
        # df['Reference sequence'] = ref_seq
        df = mtool.build_aligned_sequence_for_indelphi(self.ref_seq, 39, df)
        self.save_df(df, self.indelphi_fn)
    
    def map_indelphi_with_possible(self):
        possible_df = pd.read_csv(self.possible_fn)
        indelphi_df = pd.read_csv(self.indelphi_fn)
        df = mtool.map_indephi_with_possibleOut(possible_df, indelphi_df)
        self.save_df(df, self.possible_fn)
    
    def build_noMutated_aligned_seq_readfile(self, cutsite = 40, window_size=3):
        # read_fn, target_site_seq = fname_seq
        if self.reads_type != 'Treated' or self.reads_type != 'VEGFA3':
            if not os.path.exists(self.reads_fn):
                return (self.reads_fn, 'No control')
        
        target_site_seq = self.ref_seq
        read_df = pd.read_csv(self.reads_fn, delimiter="\t")
        # idx start from 0
        cutsite = cutsite-1
        # target_site = read_df.iloc[0]['Reference_Sequence']
        # print(target_site)
        cutting_win_seq = []
        cutting_loc = []
        loc = None
        for idx, row in read_df.iterrows():
            aligned_seq = row['Aligned_Sequence']
            win_seq = ""
            if row['n_deleted'] > 0:
                del_len = row['n_deleted']
                # print(del_len)
                #  match_start = aligned_seq.find("-"*del_len)
                match_sites = re.finditer("-"*del_len, aligned_seq)
                flag = 0
                for ms in match_sites:
                    span = ms.span()
                    # print(span)
                    if span[0] <= cutsite+1 and span[1]-1 >= cutsite-1:
                        loc = str(span[0])+":"+str(span[1])
                        # start = span[0]-window_size
                        # end = span[1]+window_size
                        start, end = span[0], span[1]
                        #  print(aligned_seq, aligned_seq[start:end])
                        win_seq = target_site_seq[:start] + "-"*del_len + target_site_seq[end:]
                        # flag = 1 
            elif row['n_inserted'] == 1 :
                # ins_len = row['n_inserted']
                #########
                ins_aligned_seq = aligned_seq[cutsite-window_size:cutsite+window_size+1]
                ins_ncl = aligned_seq[cutsite]
                loc = str(cutsite+2)+ ':' + str(cutsite+3)
                win_seq = target_site_seq[:cutsite+1] + aligned_seq[cutsite+1] + target_site_seq[cutsite+1:]
                # int(win_seq, aligned_seq[cutsite+1])
            else:
                win_seq = aligned_seq
                loc = None
                # print(win_seq)
            cutting_win_seq.append(win_seq)
            cutting_loc.append(loc)
        read_df['cutting_seq'] = cutting_win_seq
        read_df['indel_location'] = cutting_loc
        # read_df['norm_Reads'] = read_df['#Reads']/sum(read_df['#Reads'])
        # return read_df
        """
        if self.reads_type == 'Treated':
            tem_reads_df_path = "./CRISPRessoWGS_VEGFA3_win80_reads_complement/" + self.reads_fn.split("/")[-1]
        else:
            tem_reads_df_path = "./CRISPRessoWGS_DNMT1_win80_reads_complement/" + self.reads_fn.split("/")[-1]
        """
        # print(tem_reads_df_path)
        self.save_df(read_df, self.reads_fn, sep_mark="\t")
        return self.reads_fn
    
    def mark_possible_gentype_idx(self):
            
        df = pd.read_csv(self.possible_fn)
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
        # print(len(df), df.iloc[-1]['genotype_idx'])
        self.save_df(df, self.possible_fn)
        return self.possible_fn
    
    def mapping_read_with_possible(self):
        read_type = self.reads_type
        possible_df = pd.read_csv(self.possible_fn)
        possible_df['match_idx_' + read_type] = ""
        # possible_df['norm_Reads_combined'] = 0
        possible_df['#Reads_total_' + read_type] = 0
        # possible_df['Matched_reads'] = 0
        # possible_df['Mismatched_reads'] = 0
        # indel_reads = reads_df[(reads_df['n_deleted']>0)|(reads_df['n_inserted'])]
        if self.reads_type != 'Treated':
            if not os.path.exists(self.reads_fn):
                self.save_df(possible_df, self.possible_fn)
                return (self.reads_type, 'No control')
            
        
        # tn = re.split("/", reads_fn)[-1]
        reads_df = pd.read_csv(self.reads_fn, delimiter="\t")
        reads_df['match_flag'] = 0
        # reads_df['match_idx'] = 0
        # reads_df['match_len_possible'] = -1
        
        for idx_p, row_p in possible_df.iterrows():
            indel_len = row_p['Length']
            if row_p['Category'] == 'del':
                pt = re.compile("-"*indel_len)
                match_res = list(re.finditer(pt, row_p['Aligned_sequence']))
                # print(match_res, indel_len)     
                start_idx, end_idx = match_res[0].span()
                indel_loc_p = str(start_idx)+":"+str(end_idx)
                sl_reads_df = reads_df[reads_df['n_deleted'] == indel_len]

            elif row_p['Category'] == 'ins':
                cutsite = int(len(row_p['Aligned_sequence'])/2-1)
                start_idx, end_idx = cutsite, cutsite+1
                indel_loc_p = str(start_idx+2)+":"+str(end_idx+2)
                sl_reads_df = reads_df[reads_df['n_inserted'] == 1]
            else:
                continue

            for idx_r, row_r in sl_reads_df.iterrows():
                if indel_loc_p == row_r['indel_location']:
                    aligned_seq_r = row_r['cutting_seq']
                    aligned_seq_p = row_p['Aligned_sequence']
                    if aligned_seq_p == aligned_seq_r:
                        possible_df.loc[idx_p, '#Reads_total_' + read_type] += row_r['#Reads']
                        # possible_df.loc[idx_p, 'Mismatched_reads'] += row_r['#Reads']
                        reads_df.loc[idx_r, 'match_flag'] += 1
                        # reads_df.loc[idx_r, 'match_len_possible'] = row_p['Length']
                        possible_df.loc[idx_p, 'match_idx_'+read_type] += str(idx_r) + "/"
                        # reads_df.loc[idx_r, 'match_idx'] += 1
                        # print(window_seq_r, window_seq_p)
        self.save_df(reads_df, self.reads_fn, sep_mark='\t')
        self.save_df(possible_df, self.possible_fn)
        return self.check_read_possible_mapping(reads_df, possible_df)
    
    def check_read_possible_mapping(self, reads_df, possible_df):
        rn = len(reads_df[(reads_df['n_deleted']>0)|(reads_df['n_inserted'] == 1)])
        mrn = len(reads_df[reads_df['match_flag']>0])
        pn = len(possible_df[possible_df['#Reads_total_' + self.reads_type]>0])
        # reads_df.to_csv(reads_fn, sep="\t", index=False)
        # possible_df.to_csv(possible_fn, index=False)
        return (self.possible_fn, rn, mrn,  pn)
    
        
class stats_WGS_data:
    def __init__(self, reads_treated_fn=None, reads_control_fn=None, possible_fn=None):
        
        self.reads_treated_fn = reads_treated_fn
        self.reads_control_fn = reads_control_fn
        self.possible_fn = possible_fn
        self.target_name = re.split("/|.csv", possible_fn)[-2]
        if not (reads_control_fn is None):    
            if os.path.exists(reads_control_fn):
                self.withControl = True
            else:
                self.withControl = False
        
    def wgs_analysis(self):
        info_dict = {}
        possible_df = pd.read_csv(self.possible_fn)
        treated_reads_df = pd.read_csv(self.reads_treated_fn, delimiter='\t')
        # Treated
        treated_dict = self.get_corr_and_cutting_rate(treated_reads_df, possible_df, 'Treated')
        if self.withControl:
            control_reads_df = pd.read_csv(self.reads_control_fn, delimiter='\t')
            control_dict = self.get_corr_and_cutting_rate(control_reads_df, possible_df, 'Control')
            return {**treated_dict, **control_dict}
        else:
            return treated_dict

    def get_corr_and_cutting_rate(self, reads_df, possible_df, trg_type):
        res_dict = {}
        len_threshold = 90
        idx_f = (reads_df['n_inserted'] <= 1)|(reads_df['n_deleted'] <= len_threshold)
        f_read_df = reads_df[idx_f]
        sum_reads = sum(f_read_df['#Reads'])
        map_reads = sum(f_read_df[f_read_df['match_flag'] > 0]['#Reads'])

        # print(df)
        # possible_df = possible_df.groupby('genotype_idx')

        # idx_f = (possible_df['Length'] <= len_threshold)&(possible_df['withRightLen'])
        idx_f = possible_df['Length'] <= len_threshold
        f_possible_df = possible_df[idx_f]
        gf_possible_df = f_possible_df.groupby('genotype_idx').sum()

        """
        # remove the alleles of Treated that appears in Control
        try:
            false_indels_idx = gf_possible_df['control_flag']<0
            false_indels_reads = sum(gf_possible_df.loc[false_indels_idx, 'Reads_combined'])
            true_indels_reads = sum(gf_possible_df.loc[gf_possible_df['control_flag']>0, 'Reads_combined'])
            gf_possible_df.loc[false_indels_idx, 'Reads_combined'] = 0
        except:
            print(name, 'no control')
        """
        norm_read = gf_possible_df['#Reads_total_' + trg_type]/sum_reads
        edit_rate = float(map_reads/sum_reads)

        res_dict['cutting_rate_'+ trg_type] = edit_rate
        freq = gf_possible_df['Predicted frequency']/100.0
        # weighted_norm_reads = sum(freq*norm_read)
        corr = stats.spearmanr(norm_read, freq)
        kcorr, kpvalue = stats.kendalltau(norm_read, freq)

        #  print(corr)
        res_dict['sum_indel_reads_' + trg_type] = map_reads
        # refseq_df.loc[idx, 'weighted_norm_reads'] = weighted_norm_reads
        res_dict['sum_reads_' + trg_type] = sum_reads
        res_dict['#indels_' + trg_type] = len(f_read_df[(f_read_df['n_deleted']>0)|(f_read_df['n_inserted']>0)])
        # res_dict['sum_indelphi'] = sum(gf_possible_df[gf_possible_df['Reads_combined']>0]['Predicted frequency'])
        cr = corr.correlation
        if cr == cr:
            res_dict['correlationS_'+ trg_type] = cr
            res_dict['correlationS_pvalue_' + trg_type] = corr.pvalue
        else:
            res_dict['correlationS_'+ trg_type] = 0.0
            res_dict['correlationS_pvalue_' + trg_type] = 1.0
        if kcorr == kcorr:
            res_dict['correlationK_'+ trg_type] = kcorr
            res_dict['correlationK_pvalue_' + trg_type] = kpvalue
        else:
            res_dict['correlationK_'+ trg_type] = 0.0
            res_dict['correlationK_pvalue_' + trg_type] = 1.0
        return res_dict
    
    
    def get_corr_cutting_rate_removing_control_indel(self):
        
        trg_type = 'TreatedNoControl'
        res_dict = {}
        len_threshold = 90
        
        reads_df = pd.read_csv(self.reads_treated_fn, delimiter='\t')
        possible_df = pd.read_csv(self.possible_fn)
        
        idx_f = (reads_df['n_inserted'] <= 1)|(reads_df['n_deleted'] <= len_threshold)
        f_read_df = reads_df[idx_f]
        sum_reads = sum(f_read_df['#Reads'])
        # map_reads = sum(f_read_df[f_read_df['match_flag'] > 0]['#Reads'])
        
        # print(df)
        # possible_df = possible_df.groupby('genotype_idx')

        idx_f = (possible_df['Length'] <= len_threshold)&(possible_df['NotInControl'])
        f_possible_df = possible_df[idx_f]
        gf_possible_df = f_possible_df.groupby('genotype_idx').sum()
        
        map_reads = sum(gf_possible_df['#Reads_total_Treated'])
        
        norm_read = gf_possible_df['#Reads_total_Treated']/sum_reads
        edit_rate = 1.0*map_reads/sum_reads

        res_dict['cutting_rate_'+ trg_type] = edit_rate
        freq = gf_possible_df['Predicted frequency']/100.0
        # weighted_norm_reads = sum(freq*norm_read)
        corr = stats.spearmanr(norm_read, freq)
        kcorr, _ = stats.kendalltau(norm_read, freq)

        res_dict['sum_indel_reads_' + trg_type] = map_reads
        # refseq_df.loc[idx, 'weighted_norm_reads'] = weighted_norm_reads
        # res_dict['sum_reads_' + trg_type] = sum_reads
        # res_dict['sum_indelphi'] = sum(gf_possible_df[gf_possible_df['Reads_combined']>0]['Predicted frequency'])
        cr = corr.correlation
        if cr == cr:
            res_dict['correlationS_'+ trg_type] = cr
        else:
            res_dict['correlationS_'+ trg_type] = 0.0
        if kcorr == kcorr:
            res_dict['correlationK_'+ trg_type] = kcorr
        else:
            res_dict['correlationK_'+ trg_type] = 0.0
        return res_dict
        
    
    def mark_overlapped_indels_treated_ctrl(self):
        if self.possible_fn is None:
            print("Please provide possible allels file!")
            return 0
        possible_df = pd.read_csv(self.possible_fn)
        possible_df['NotInControl'] = 'True'
        
        if "#Reads_total_Control" not in possible_df.columns:
            return self.target_name, "no control!"

        for idx, row in possible_df.iterrows(): 
            if row['#Reads_total_Control'] > 0:
                possible_df.loc[idx, 'NotInControl'] = False
        possible_df.to_csv(self.possible_fn, index=False)
        return self.target_name, "succeed!"
    
    
    def stats_indel_types_and_TCreads_correlation(self):
        res_dict = {}
        
        def get_n_indel(df):
            cut_indel_n = len(df[df['match_flag']>0])
            cut_1indel_n = len(df[(df['match_flag']>0)&(df['#indels'] == 1)])
            indel1_n = len(df[df['#indels'] == 1])
            indelL1_n = len(df[df['#indels'] > 1])
            indelL2_n = len(df[df['#indels'] > 2])
            return indel1_n, indelL1_n, indelL2_n, cut_indel_n, cut_1indel_n
        
        reads_treated_df = pd.read_csv(self.reads_treated_fn, delimiter='\t')
        indel1_n, indelL1_n, indelL2_n, cut_indel_n, cut_1indel_n = get_n_indel(reads_treated_df)
        res_dict['#indels_cutsite_Treated'] = cut_indel_n
        res_dict['#1_indels_cutsite_Treated'] = cut_1indel_n
        res_dict['#1_indels_Treated'] = indel1_n
        res_dict['#G1_indels_Treated'] = indelL1_n
        res_dict['#G2_indels_Treated'] = indelL2_n
        
        res_dict['#Reads_correlation'] = 0.0
        if os.path.exists(self.reads_control_fn):
            reads_control_df = pd.read_csv(self.reads_control_fn, delimiter='\t')
            indel1_n, indelL1_n, indelL2_n, cut_indel_n, cut_1indel_n = get_n_indel(reads_control_df)
            res_dict['#indels_cutsite_Control'] = cut_indel_n
            res_dict['#1_indels_cutsite_Control'] = cut_1indel_n
            res_dict['#1_indels_Control'] = indel1_n
            res_dict['#G1_indels_Control'] = indelL1_n
            res_dict['#G2_indels_Control'] = indelL2_n
            
            possible_df = pd.read_csv(self.possible_fn)
            treated_reads = possible_df['#Reads_total_Treated']
            control_reads = possible_df['#Reads_total_Control']
            
            f_read_df = reads_treated_df[(reads_treated_df['n_inserted'] <= 1)|(reads_treated_df['n_deleted'] <= 90)]
            sum_treated_reads = sum(f_read_df['#Reads'])    
            f_read_df = reads_control_df[(reads_control_df['n_inserted'] <= 1)|(reads_control_df['n_deleted'] <= 90)]
            sum_control_reads = sum(f_read_df['#Reads'])
            
            treated_reads = treated_reads/sum_treated_reads
            control_reads = control_reads/sum_control_reads
            
            corr, pvalue = stats.spearmanr(treated_reads, control_reads)
            if corr == corr:
                res_dict['#Reads_correlation'] = corr
            else:
                res_dict['#Reads_correlation'] = 0.0
            
        return res_dict


    @classmethod
    def get_number_indels_possible_df(self, possible_fn):
        trg_name = re.split("/|.csv", possible_fn)[-2]
        res_dict = {}
        len_threshold = 90
        possible_df = pd.read_csv(possible_fn)
        if not "NotInControl" in possible_df.columns:
            return (trg_name, -1)
        idx_f = (possible_df['Length'] <= len_threshold)&(possible_df['NotInControl'])
        f_possible_df = possible_df[idx_f]
        gf_possible_df = f_possible_df.groupby('genotype_idx').sum()
        number_indels = len(gf_possible_df[gf_possible_df['#Reads_total_Treated']>0])
        return (trg_name, number_indels)
    
    @classmethod
    def get_indels_n_reads_df(self, reads_fn):
        if not os.path.exists(reads_fn):
            return (reads_fn, "not exists!")
        
        reads_df = pd.read_csv(reads_fn, delimiter="\t")
        
        n_del_sites = []
        n_ins_sites = []
        for idx, row in reads_df.iterrows():
            aligned_seq = row['Aligned_Sequence']
            ref_seq = row['Reference_Sequence']
            del_n = count_indels(aligned_seq)
            ins_n = count_indels(ref_seq)
            n_del_sites.append(del_n)
            n_ins_sites.append(ins_n)
            
        reads_df['#deletions'] = n_del_sites
        reads_df['#insertions'] = n_ins_sites
        reads_df['#indels'] = reads_df['#deletions'] + reads_df['#insertions']
        reads_df.to_csv(reads_fn, sep="\t")
        return (reads_fn, "finished!")
    
    @classmethod
    def count_indels(self, seq):
        indel_n = 0
        started = False
        for char in seq:
            if char == '-' and not started:
                indel_n += 1
                started = True
            elif char != '-' and started:
                started = False
        return indel_n
                