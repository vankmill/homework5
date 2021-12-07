###
#   name    | mary lauren benton
#   created | 2021
#   input   | putative sirta sequence (ACGT)
#   output  | score
#
#   this program takes a DNA sequence and optional penalty (default 1) and scores
#   it based on TG composition using Kathy's rules
###

import argparse
import math
import regex as re
import numpy as np
from Bio import SeqIO


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate SiRTA score for sequence.")

arg_parser.add_argument("sequence_file", help='FASTA file of putative SiRTA sequences')
arg_parser.add_argument("-p", "--penalty", type=int, default=1,
                        help='penalty to apply for each GGTGG occurence; default=1')
arg_parser.add_argument("-t", "--ttpenalty", type=int, default=1,
                        help='penalty to apply for each flanking TT occurence; default=1')
arg_parser.add_argument("-w", "--window", type=int, default=100,
                        help='sliding window size; default=100bp')
arg_parser.add_argument("-s", "--step", type=int, default=1,
                        help='step size for sliding windows; default=1bp')
arg_parser.add_argument('--signal', dest='signal', action='store_true', default=False,
                        help='flag to print signal output instead of max score')
arg_parser.add_argument('--thresh', type=float, default=0.0,
                        help='when used with --signal, only return scores above this value; default=0')
arg_parser.add_argument('--strand', dest='strand', type=str, default='+', choices=['+', '-'],
                        help='flag to specify strand with --signal; default +')


###
# functions
###
# test the match to rule out sequences of all gs
def all_gs(seq):
    c = seq[0]
    for s in seq:
        if s != c:
            return False
    return True


# test the candidate sequence to ensure that it has no internal TT
def test_candidate(c):
    result = []
    tt_pos = [s.start(0) for s in re.finditer(r"TT", c)]
    if len(tt_pos) == 0 and len(c) >= 4:
        result.append(c)
    elif len(c) > 4:
        result += test_candidate(c[0:tt_pos[0]+1])
        result += test_candidate(c[tt_pos[0]+1:])
    return result


# find all strings of 4 or more nucleotides containing only consecutive G/T
# Ts must be single, Gs can be 1, 2, or 3 consecutive nucleotides
def find_tg_kmers(seq):
    kmers, pos = [], []

    for match in re.finditer(r'(?!TT)([GT]{4,})', seq):
        if not all_gs(seq[match.start():match.end()]):
            ks = test_candidate(seq[match.start():match.end()])
            for k in ks:
                if k != "":
                    kmers.append(k)
                    if match.end() - match.start() == len(k):
                        pos.append((match.start(), match.end()))
                    else:
                        offset = re.search(k, seq[match.start():match.end()]).start()
                        new_start = match.start()+offset
                        pos.append((new_start, new_start+len(k)))
    return kmers, pos


# count the number of nucleotides found using the find_tg_kmers function
def score_groups(seq_groups):
    return sum([len(x) for x in seq_groups])


# adjust the score by subtracting [penalty p] for each instance of GGTGG
def calculate_penalty(seq_groups, p, tp, pos, seq):
    penalty_sum = 0
    for g in seq_groups:
        # GGTGG penalty
        ggtgg = re.findall(r"(GGTGG)", g, overlapped=True)
        if ggtgg is not None:
            penalty_sum += (len(ggtgg)*p)

    # flanking TT penalty (could be optimized)
    for start, end in pos:
        if start-1 >= 0:
            if seq[start] == "T" and seq[start-1] == "T":
                penalty_sum += tp
        if end < len(seq):
            if seq[end] == "T" and seq[end-1] == "T":
                    penalty_sum += tp
    return penalty_sum


# calculate score across sliding window, return the max score of the windows
def calc_max_score_over_windows(seq, penalty, tt_penalty, window, step):
    num_windows = round(((len(seq)-window)/step)+1)
    max_score = np.nan

    # slide windows along sequence
    for i in range(0, num_windows*step, step):
        seq_len = len(seq[i:i+window])

        # calculate score
        groups, pos = find_tg_kmers(seq[i:i+window])
        score = score_groups(groups) - calculate_penalty(groups, penalty, tt_penalty, pos, seq[i:i+window])
        if np.isnan(max_score):
            max_score = score
        elif score > max_score:
            max_score = score
    return max_score


# parse the fasta file and calculate scores, assumes 1 seq per line
def parse_seqs_from_fasta_max(seq_file, penalty, tt_penalty, win_size, step_size):
    for record in SeqIO.parse(seq_file, 'fasta'):
        print(record.id)
        print(calc_max_score_over_windows(str(record.seq), penalty, tt_penalty, win_size, step_size))


# calculate score over sliding windows, return all scores
def calc_score_over_windows(seq, chrom, start, penalty, tt_penalty, window, step, strand):
    num_windows = math.ceil(((len(seq)-window)/step)+1)
    scores = []

    # slide windows along sequence
    for i in range(0, num_windows*step, step):
        seq_len = len(seq[i:i+window])

        # calculate score
        groups, pos = find_tg_kmers(seq[i:i+window])
        score = score_groups(groups) - calculate_penalty(groups, penalty, tt_penalty, pos, seq[i:i+window])
        if strand == '+':
            scores.append([chrom, start+i, start+i+seq_len, score])
        else:
            scores.append([chrom, start-i-seq_len, start-i, score])
    return scores


# parse the fasta file and calculate scores, assumes 1 seq per line
def parse_seqs_from_fasta_signal(seq_file, penalty, tt_penalty, win_size, step_size, thresh, strand):
    final_scores = []

    # parse fasta, one sequence per line (bedtofasta result default)
    for record in SeqIO.parse(seq_file, 'fasta'):
        print(f'#{record.id}')
        header = re.split(":|-", record.id)
        chrom, start = header[0][0:3].lower()+header[0][3:], header[1]
        scores = calc_score_over_windows(str(record.seq), chrom=chrom, start=int(start), penalty=penalty,
                                         tt_penalty=tt_penalty, window=win_size, step=step_size, strand=strand)
        signal_to_bedgraph(scores, thresh)


def signal_to_bedgraph(signals, threshold):
    for s in signals:
        if s[3] >= threshold:
            print(f'{s[0]}\t{s[1]}\t{s[2]}\t{s[3]}')


# main function
def main():
    args = arg_parser.parse_args()

    # save parameters
    SIRTA_FILE = args.sequence_file
    PENALTY = args.penalty
    TT_PENALTY = args.ttpenalty
    WINDOW = args.window
    STEP = args.step
    SIGNAL = args.signal
    THRESH = args.thresh
    STRAND = args.strand

    # parse fasta and calculate score
    if SIGNAL:
        s = parse_seqs_from_fasta_signal(SIRTA_FILE, PENALTY, TT_PENALTY,
                                         WINDOW, STEP, THRESH, STRAND)
    else:
        parse_seqs_from_fasta_max(SIRTA_FILE, PENALTY, TT_PENALTY, WINDOW, STEP)


###
# main
###
if __name__ == '__main__':
    main()
