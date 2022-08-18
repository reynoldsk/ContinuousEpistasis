"""
Identify index, barcoded sgRNA constructs from Illumina HiSeq FASTQ files. The following code takes indexed, paired-end sequencing files, locates two sgRNA sequences and a barcode sequence, and assigns a pairwise construct ID. sgRNA and barcode sequences must match their expected sequence exactly, and criticial sgRNA nucleotides (those that differentiate distinct sgRNAs targeting the same genomic locus) must have a Q-score greater than 30 to be called. Outputs sgRNA counts before and after applying the applicable filter. To get counts without any filter, set qthreshold = 0.

Created on 2/15/21 by Agata Turska-Nowak
Edited: 4/30/21 Ryan Otto
Edited: 8/27/21 Ryan Otto
Edited: 6/29/22 Ryan Otto
"""
# Import packages and libraries
import regex as re
import numpy as np
import pandas as pd
import sys
from Bio import SeqIO
from itertools import islice
from Bio.Seq import Seq
date = '20210914'

pairwise_lib = {}
# From the fasta file take the sequences and their id to create library with sgRNA guides
sgRNA2seq = {}
seq2sgRNA = {}
for record in SeqIO.parse('sgRNAs_seq.fa', 'fasta'):
    seq2sgRNA[record.seq[35:55].reverse_complement()] = record.id
    sgRNA2seq[record.id] = record.seq[35:55].reverse_complement()
guides_list = list(sgRNA2seq.keys())
# Information about the crucial nt in sgRNAs sequences to be check in q30 filter step
nt_sgRNA = {}
with open('nt_sgRNA.txt') as file:
    line = file.readline()
    count = 0
    while line:
        line = line.strip('\n')
        line_split = line.split(':')
        nt_sgRNA[line_split[0]] = []
        number_nt = line_split[1].split(',')
        nt_sgRNA[line_split[0]].append(number_nt[0])
        nt_sgRNA[line_split[0]].append(number_nt[1])
        line = file.readline()
        count += 1
# Read through all files with reads from each timepoint
file_names = open(sys.argv[1]).readlines()
for i in range(len(file_names)):
    file_names[i] = file_names[i].strip('\n')
print(file_names)

def seq_find(seq, qscore, left, right, target_len, hamming_dist):
    """Identifies a target sequence using regular expressions. Sequence must be of specified length, and 
    gaps
    cannot exist in the flanking sequences.
    Arguments:
    seq: Full sequence
    qscore: Full qscore array
    left: Left flanking sequence
    right: Right flanking sequence
    target_len: Length of desired sequence
    hamming_dist: Acceptable number of mismatches in flanking regions (no gaps)
    Returns:
    target_seq: If found, returns the desired sequence.
    target_qscore: If found, returns the desired sequence's Q-scores.
    """
    make_string = '(' + left + '.'*target_len + right + ')' +  '{s<' + str(hamming_dist) + '}'
    matches = re.finditer(make_string, str(seq))
    for i, match in enumerate(matches):
        if i > 0:
            target_seq = 'Not found'
            target_qscore = 'Not found'
            return target_seq, target_qscore
        found_seq = seq[match.start():match.end()]
        found_qscore = qscore[match.start():match.end()]
    try:
        target_seq = found_seq[len(left):-len(right)]
        target_qscore = found_qscore[len(left):-len(right)]
    except:
        target_seq = 'Not found'
        target_qscore = 'Not found'
    return target_seq, target_qscore


def qscore_filter(seq_qscores, qthreshold):
    """Filters sgRNAs based on Q-score confidence.
    Arguments:
    seq_qscores: Sequence's Q-scores
    qthreshold: Confidence threshold
    Returns:
    high_quality: Boolean stating whether the sequence passed the Q-score filter
    """
    high_quality = True
    for base in seq_qscores:
        if (ord(base) -33) <= qthreshold:
            high_quality = False
            break
    return high_quality


def check_qscore(sgRNA1_name, qscore_sgRNA1, sgRNA2_name, qscore_sgRNA2, qthreshold):
    """Second portion of Q-score filter that determines which sgRNA bases must be passed through the 
    filter.
    Arguments:
    sgRNA1_name: Presumed sgRNA1 identity
    qscore_sgRNA1: sgRNA1 Q-scores
    sgRNA2_name: Presumed sgRNA2 identity
    qscore_sgRNA2: sgRNA2 Q-scores
    qthreshold: Confidence threshold
    Returns:
    a: Boolean stating whether the sequence 1 passed the Q-score filter
    b: Boolean stating whether the sequence 2 passed the Q-score filter
    """
    if sgRNA1_name in nt_sgRNA.keys():
        if sgRNA1_name in nt_sgRNA.keys():
            for nt1 in nt_sgRNA[sgRNA1_name]:
                a = qscore_filter(qscore_sgRNA1[int(nt1)-1], qthreshold)
                if not a:
                    break
            for nt2 in nt_sgRNA[sgRNA2_name]:
                b = qscore_filter(qscore_sgRNA2[int(nt2)-1], qthreshold)
                if not b:
                    break
    else:
        a = False
        b = False
    return a, b


sgRNA1_dict = {}
qscore_sgRNA1_dict = {}
sgRNA2_dict = {}
qscore_sgRNA2_dict = {}
BC_dict = {}
qscore_BC_dict = {}
sequence = {}
q_score = {}
read_total = 0
for read_file in file_names:
    read_1_file = f'{read_file}_R1_001.fastq'
    read_2_file = f'{read_file}_R2_001.fastq'
    print(read_file, flush=True)
    sp = read_file.split('/')
    file_id = sp[1] 
    sgRNA1_dict[file_id] = []
    qscore_sgRNA1_dict[file_id] = []
    sgRNA2_dict[file_id] = []
    qscore_sgRNA2_dict[file_id] = []
    BC_dict[file_id] = []
    qscore_BC_dict[file_id] = []
    with open(read_1_file) as fwdFile:
        while True:
            next_n_lines = list(islice(fwdFile,4))
            if not next_n_lines:
                break
            sequence = next_n_lines[1]
            q_score = next_n_lines[3]
            read_total =+ 1
            sgRNA2_seq, qscore_sgRNA2 = seq_find(sequence, q_score, 'CTAGCTCTAAAAC', 'A', 20, 3)
            BC_seq, qscore_BC = seq_find(sequence, q_score, 'GTACAGCGAGGCAAC', 'ACGGATCCCCAC', 6, 3)
            sgRNA2_dict[file_id].append(sgRNA2_seq)
            qscore_sgRNA2_dict[file_id].append(qscore_sgRNA2)
            BC_dict[file_id].append(BC_seq)
            qscore_BC_dict[file_id].append(qscore_BC)
    with open(read_2_file) as fwdFile:
        while True:
            next_n_lines = list(islice(fwdFile,4))
            if not next_n_lines:
                break
            sequence = Seq(next_n_lines[1]).reverse_complement()
            q_score = next_n_lines[3][::-1]
            read_total =+ 1
            sgRNA1_seq, qscore_sgRNA1 = seq_find(sequence, q_score, 'CTAGCTCTAAAAC', 'A', 20, 3)
            sgRNA1_dict[file_id].append(sgRNA1_seq)
            qscore_sgRNA1_dict[file_id].append(qscore_sgRNA1)
all_seqs_total = 0
reads_dict = {}
for read_file in file_names:
    sp = read_file.split('/')
    file_id = sp[1]
    reads_dict[file_id] = []
    sgRNA1_sum = sum([x != 'Not found' for x in sgRNA1_dict[file_id]])
    sgRNA2_sum = sum([x != 'Not found' for x in sgRNA2_dict[file_id]])
    BC_sum = sum([x != 'Not found' for x in BC_dict[file_id]])
    total = len(sgRNA1_dict[file_id])
    reads_dict[file_id].append(total)
    print(sgRNA1_sum)
    print(sgRNA2_sum)
    print(BC_sum)
    print(f'Total: {total}')
    both = [len(x) == 20 and len(y) == 20 and len(z) == 6 for x, y, z in zip(sgRNA1_dict[file_id], sgRNA2_dict[file_id],
                                                                             BC_dict[file_id])]
    reads_dict[file_id].append(sum(both))
    print(f'Success: {sum(both)}')
    print('')
    all_seqs_total += sum(both)
ab = 0
onlya = 0
onlyb = 0
none = 0
failed = 0
qthreshold = 30
for read_file in file_names:
    sp = read_file.split('/')
    file_id = sp[1]
    pairwise_lib[file_id] = {}
    for BC in ['TGAAAG', 'CATGAT', 'CCATGC', 'TCATAC', 'TAGACT', 'ACTAGG', 'Other']:
        pairwise_lib[file_id][BC] = pd.DataFrame(np.zeros([len(guides_list), len(guides_list)]),\
                                               guides_list,guides_list)
    for i, guide1 in enumerate(sgRNA1_dict[file_id]):
        guide2 = sgRNA2_dict[file_id][i]
        BC = BC_dict[file_id][i]
        if len(guide1) == 20 and len(guide2) == 20:
            if BC in ['TGAAAG', 'CATGAT', 'CCATGC', 'TCATAC', 'TAGACT', 'ACTAGG']:
                try:
                    guide1_id = seq2sgRNA[guide1]
                    guide2_id = seq2sgRNA[guide2]
                    a, b = check_qscore(guide1_id, qscore_sgRNA1_dict[file_id][i], guide2_id, qscore_sgRNA2_dict[file_id][i], qthreshold)
                    if a and b:
                        pairwise_lib[file_id][BC].loc[guide1_id][guide2_id] = 1 + \
                            pairwise_lib[file_id][BC].loc[guide1_id][guide2_id]
                        ab += 1
                    elif a and not b:
                        onlya += 1
                        continue
                    elif b and not a:
                        onlyb += 1
                        continue
                    else:
                        none += 1
                        continue
                except:
                    failed += 1
            else:
                try:
                    guide1_id = seq2sgRNA[guide1]
                    guide2_id = seq2sgRNA[guide2]
                    pairwise_lib[file_id]['Other'].loc[guide1_id][guide2_id] = 1 + \
                    pairwise_lib[file_id]['Other'].loc[guide1_id][guide2_id]
                except:
                    failed += 1
    sum_low = onlya + onlyb + none
    sum_all = sum_low + ab
    print('timepoint: ' + str(file_id))
    print('qscore > ' + str(qthreshold) + ' for both sgRNAs: ' + str(ab))
    print('qscore < ' + str(qthreshold) + ' for sgRNA in position 1: ' + str(onlya))
    print('qscore < ' + str(qthreshold) + ' for sgRNA in position 2: ' + str(onlyb))
    print('qscore < ' + str(qthreshold) + ' for both sgRNAs: ' + str(none))
    print('sum < ' + str(qthreshold) + ': ' + str(sum_low))
    print('all successful reads: ' + str(sum_all))
    print('total reads: ' + str(sum_all + failed))
for read_file in file_names:
    sp = read_file.split('/')
    file_id = sp[1]
    for BC in pairwise_lib[file_id].keys():
        guides_list.sort()
        pairwise_lib[file_id][BC] = pairwise_lib[file_id][BC][guides_list].T
        pairwise_lib[file_id][BC] = pairwise_lib[file_id][BC][guides_list].T
for TP in pairwise_lib.keys():
    for BC in pairwise_lib[TP].keys():
        pd.DataFrame(pairwise_lib[TP][BC]).to_csv('BarcodedCounts/q30_filter/' + date + '_' + TP + '_' + BC + '_counts.csv')