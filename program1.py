import sys
import numpy as np
import pandas as pd
from shuffle_opt import *
import Bio
from Bio import AlignIO
from Bio import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from collections import OrderedDict
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

method = sys.argv[2]

gap_open = -12
gap_extend = -3

original_seq_fasta = ''
optimized_seq_fasta = ''

def open_fasta(in_string):
    seq_dict = OrderedDict()
    with open(in_string,'r') as f:
        outfile = f.read().splitlines()
        i = len(outfile)
        j = 0
        aa_seq = ''
        aa_name = ''
        for j in range(i):
            if outfile[j][0] == '>':
                aa_name = outfile[j][1:]
            elif i == j + 1:
                aa_seq += outfile[j]
                seq_dict[aa_name] = aa_seq
            elif outfile[j+1][0] == '>':
                aa_seq += outfile[j]
                seq_dict[aa_name] = aa_seq
                aa_seq = ''
            else:
                aa_seq += outfile[j]
    return seq_dict

def create_tracking_dict(seq_dict):
    count = 0
    tracking_dict = {}
    for i in seq_dict:
        tracking_dict[i] = count
        count += 1
    return tracking_dict

def create_aln_score_matrix(seq_dict, tracking_dict):
    numSeqs = len(seqDict)
    aln_score_df = pd.DataFrame(np.nan, index = range(numSeqs), columns = range(numSeqs))
    for i in seq_dict:
        for j in seq_dict:
            if tracking_dict[i] < tracking_dict[j]:
                aln_score = pairwise2.align.globalds(seq_dict[i], seq_dict[j], matrix, gap_open, gap_extend, penalize_end_gaps=0, score_only=True)
                aln_score_df.iloc[tracking_dict[i], tracking_dict[j]] = aln_score
                aln_score_df.iloc[tracking_dict[j], tracking_dict[i]] = aln_score
    return aln_score_df

def create_codon_usage_dict():
    codon_usage_dict = OrderedDict()
    with open('codon_usage_ecoli.csv','r') as f:
        a = f.read().splitlines()
        for i in a:
            b = i.split(',')
            if b != ['']:
                if b[1] not in codon_usage_dict:
                    codon_usage_dict[b[1]] = OrderedDict()
                    codon_usage_dict[b[1]][float(b[2])] = b[0]
                else:
                    codon_usage_dict[b[1]][float(b[2])] = b[0]
    return codon_usage_dict


def back_translate(seq, codon_usage_dict):
    dna_seq = ''
    for aa in seq:
        randNum = np.random.rand()
        check = 0
        cumFreq = 0
        for codonFreq in codon_usage_dict[aa]:
            cumFreq += codonFreq
            if randNum < cumFreq and check != 1:
                dna_seq += codon_usage_dict[aa][codonFreq]
                check = 1
        if check == 0:
            dna_seq += codon_usage_dict[aa][codonFreq]
    return dna_seq

def generate_opt_seq_df(seq_dict, codon_usage_dict, tracking_dict):
    numSeqs = len(seq_dict)
    opt_seq_df = pd.DataFrame(np.nan, index = range(numSeqs), columns = range(numSeqs))
    for i in seq_dict:
        b_translate_opt_seq = back_translate(seq_dict[i], codon_usage_dict)
        for j in seq_dict:
            if i != j:
                opt_seq = codon_optimize(back_translate(seq_dict[j], codon_usage_dict), b_translate_opt_seq)
                opt_seq_df.iloc[tracking_dict[i], tracking_dict[j]] = opt_seq
            else:
                opt_seq_df.iloc[tracking_dict[i], tracking_dict[j]] = b_translate_opt_seq
    return opt_seq_df

def get_weighted_bases(opt_seq_df, aln_score_df, tracking_dict, seq_dict):
    original_seq_fasta = ''
    optimized_seq_fasta = ''

    for i in range(len(opt_seq_df)):
        opt_nuc_seq = ''
        sumAlignScore = aln_score_df.iloc[i].sum(axis=0)
        weighted_seq_dict = OrderedDict()
        originName = tracking_dict.keys()[tracking_dict.values().index(i)]
        seqLength = len(seq_dict[originName])*3
        nuc_keys = ['A','T','G','C']
        for ii in range(seqLength):
            weighted_seq_dict[ii] = {key: 0 for key in nuc_keys}
        for j in range(len(opt_seq_df.iloc[i])):
            if i != j:
                if method == 'similar':
                    weightScore = aln_score_df.iloc[i][j]/sumAlignScore
                elif method == 'distant':
                    weightScore = 1 - aln_score_df.iloc[i][j]/sumAlignScore
                for ind, nuc in enumerate(opt_seq_df.iloc[i][j].replace('-','')):
                    weighted_seq_dict[ind][nuc] += weightScore
            else:
                original_seq_fasta += '>' + originName + '\n'
                original_seq_fasta += opt_seq_df.iloc[i][j] + '\n'
        for seq_ind in weighted_seq_dict:
            opt_nuc_seq += max(weighted_seq_dict[seq_ind], key=lambda k: weighted_seq_dict[seq_ind][k])
        optimized_seq_fasta += '>' + originName + '\n'
        optimized_seq_fasta += opt_nuc_seq + '\n'
    return original_seq_fasta, optimized_seq_fasta

seqDict = open_fasta(sys.argv[1])

for ind,i in enumerate(seqDict):
    if ind == 0:
        dorn1_name = i
    else:
        print (dorn1_name, i)
        codon_optimize(seqDict[dorn1_name], seqDict[i])


trackingDict = create_tracking_dict(seqDict)

alnScore_df = create_aln_score_matrix(seqDict, trackingDict)

codonDict = create_codon_usage_dict()

optSeq_df = generate_opt_seq_df(seqDict, codonDict, trackingDict)

ori, opt = get_weighted_bases(optSeq_df, alnScore_df, trackingDict, seqDict)

alnScore_df.to_csv(path_or_buf='aln_scores_' + method + '.csv', sep=',', na_rep='na')
file = open('orignal_seq_' + method + '.fa','w')
file.write(ori)
file.close

file = open('optimized_seq_' + method + '.fa','w')
file.write(opt)
file.close




