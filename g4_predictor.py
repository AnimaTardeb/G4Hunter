"""
Guanine_quadruplex_predictor.py -i <inputFasta> -o <outputDir> -w Window_width -s Threshold_score
Author: J Sietsma Penington
Date: 21 April 2017

Implementation of G4Hunter algorithm from 
Re-evaluation of G-quadruplex propensity with G4Hunter 
Amina Bedrat  Laurent Lacroix  Jean-Louis Mergny
Nucleic Acids Res (2016) 44 (4): 1746-1759. 
DOI:  https://doi.org/10.1093/nar/gkw006

Input Fasta-format file of sequences is scanned for regions with high 'G-richness' + 
'G-skewness', implemented as a score based on the number of consecutive G nucleotides.
To simultaneously assess the complementary strand the same score is calculated as a 
negative number for consecutive C nucleotides. Skewness is implemented by adding the 
positive G scores and negative C scores so that regions with balanced G and C numbers have 
small scores.

Two output files are written, one containing all sliding windows with absolute value of 
mean score above the given threshold, and one containing extended regions of merged 
windows which pass the criterion.

Intervals are written in 0-based, half open format, (]. 'start' is the position 
before the interval and 'end' is the last position in the interval. 
For example the first nucleotide would be described:
start	end	length
    0	  1	     1

Two decisions that may need changing:
	1. contrary to Bedrat-Lacroix-Mergny algorithm, fixed-width window scores are not 
	adjusted for terminal runs that extend outside the window. For example, a final 
	single G adds 1 to the score, regardless of how many Gs may follow the window.
	This can give a difference in scores for a window, to a maximum difference of 
	8/Window_width. Worst case example: if window GGA...TGG was flanked by 2 more Gs at 
	each end, the four initial and final Gs in the window score 2 each in this algorithm, 
	but 4 each in the original algorithm, giving a total difference of 8.
	
	2. A merged region of adjacent windows has the score for the region calculated, but 
	is not filtered by the final score. It is possible a number of windows with 
	scores above threshold may combine to a region with mean score below threshold, which
	would still be written to the output file. This is same as Bedrat-Lacroix-Mergny 
	output, but may be undesirable.
	 
"""

import sys
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq

def ScoreSeq(seq):
    score = 0
    priorGcount = 0
    priorCcount = 0
    for i in range(len(seq)):
        if seq[i] == "G":
            if priorGcount < 4:
                score += priorGcount*2 + 1
                priorGcount += 1
            else:   # if run reaches GGGG, each subsequent G adds 4 regardless of number
                score += 4
            priorCcount = 0
        elif seq[i] == "C":
            if priorCcount < 4:
                score -= priorCcount*2 + 1
                priorCcount += 1
            else:   # if run reaches CCCC, each subsequent C subtracts 4 regardless of number
                score -= 4
            priorGcount = 0
        else:    # nucleotide is A, T or ambiguous: no change to score, any run is terminated
            priorGcount = 0
            priorCcount = 0
    meanScore = float(score) / len(seq)
    return meanScore


def G4predictor(inFasta, out_win_file, out_region_file, win_width, thresh):
    regionCount = 0
    seqRecordCount = 0
    for seq_record in SeqIO.parse(inFasta, "fasta"):
        seqRecordCount += 1
        out_win_file.write(seq_record.id + '\n')
        out_win_file.write("Start\tEnd\tSequence\tLength\tScore\n")
        out_region_file.write(seq_record.id + '\n')
        seq_curr = seq_record.upper()
        prevWinPass = False
        G4predCount = predG4start = predG4end = 0
        predG4seq = ""
        
        # Slide window along sequence, advance 1 nucleotide each loop:
        for i in range(0, len(seq_curr) - win_width + 1):
            win_curr = seq_curr[i : i + win_width]
            score_curr = ScoreSeq(win_curr)
            if abs(score_curr) > thresh :
                out_win_file.write('%i\t%i\t%s\t%i\t%.2f\n'.format( 
                    i, i+win_width, win_curr.seq, win_width, score_curr))
                if i > predG4end:
                    # start of new predicted G4 region : write any previous region to file
                    if G4predCount == 0:
                        out_region_file.write("Start\tEnd\tSequence\tLength\tScore\tNumber\n")
                    else:
                        out_region_file.write('%i\t%i\t%s\t%i\t%.2f\t%i\n'.format(
                           predG4start, predG4end, predG4seq, 
                            len(predG4seq), ScoreSeq(predG4seq), G4predCount) )
                    G4predCount += 1
                    # Adjust start of new region: 
                    if win_curr[0] in ["C", "G"]:
                        #  extend backwards to prior Gs or Cs if they match
                        predG4seq = win_curr.seq
                        predG4start = i
                        while predG4start > 0 and predG4seq[0] == seq_curr.seq[predG4start-1]:  
                            predG4seq = seq_curr.seq[predG4start-1] + predG4seq
                            predG4start -= 1
                    else: # strip off leading A and T nucleotides
                        predG4seq = win_curr.seq.lstrip("AT")
                        predG4start = i + win_width - len(predG4seq)
                else:  # current window extends current G4 region
                    if prevWinPass:  # simple case of adjacent windows
                        predG4seq = predG4seq + win_curr.seq[-1]
                    else:   
                        predG4seq = seq_curr[predG4start:i + win_width].seq
                predG4end = i + win_width
                prevWinPass = True
            else:
                if prevWinPass:
                    # a run of 1 or more windows with scores above threshold has just ended 
                    if predG4seq[-1] in ["C", "G"]:  
                        #  extend forward to next base if it matches
                        if predG4seq[-1] == win_curr.seq[-1]:
                            predG4seq = predG4seq + win_curr.seq[-1]
                    else:  # strip off trailing A and T nucleotides
                        predG4seq = predG4seq.rstrip("AT")
                    predG4end = predG4start + len(predG4seq)
                    
                prevWinPass = False
        
        # Adjust end of, and write, last region
        if G4predCount == 0:
            out_region_file.write("# No predicted quadruplex regions found\n")
        else:
            while (predG4end < len(seq_curr)-1 and predG4seq[-1] in ["C", "G"] and 
                   predG4seq[-1] == seq_curr.seq[predG4end+1]):
                        predG4seq = predG4seq + seq_curr.seq[predG4end+1]
                        predG4end += 1
            predG4seq = predG4seq.rstrip("AT")
            predG4end = predG4start + len(predG4seq)
            out_region_file.write('%i\t%i\t%s\t%i\t%.2f\t%i\n'.format(
                    predG4start, predG4end, predG4seq, len(predG4seq), 
                     ScoreSeq(predG4seq), G4predCount) )
        regionCount += G4predCount
    print "Number of sequences:", seqRecordCount
    print "Number of potential quadruplex regions from given parameters:", regionCount


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=
    'Search a FASTA file for regions likely to form G4 quadruplexes')
    parser.add_argument('-i', dest='fasta_name', metavar='<infile>', 
                       help='the input fasta file')
    parser.add_argument('-o', dest='out_dir',  metavar='<outdir>',
                       help='directory for output files')
    parser.add_argument('-w', dest='window_width', type=int, 
                       help= 'Width of sliding window. Recommended value 25 nts')
    parser.add_argument('-s', dest='threshold', type=float, 
    help= 'Threshold for absolute value of score to predict a quadruplex region. Typical value 1.0 to 1.5')
                       

    args = parser.parse_args()
    basefname = os.path.basename(args.fasta_name).split('.')[0]
    outDirName = os.path.join(args.out_dir, "g4Results_" + str(args.threshold))
    os.mkdir(outDirName)
    outWinName = os.path.join(outDirName, 
        basefname + "_" + str(args.window_width) + "nts.tsv")
    outputWindows = open(outWinName, 'w')
    outRegionName = os.path.join(outDirName, 
        basefname + "_merged.tsv")
    outputRegions = open(outRegionName, 'w')
    
    print "Input file:", os.path.basename(args.fasta_name)
    G4predictor(args.fasta_name, outputWindows, outputRegions, args.window_width, args.threshold)
