#!/usr/bin/python3
# PROJECT: Wavefront Alignments Algorithms 
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Compare alignment (*.alg) files
# USAGE: python3 wfa.alg.cmp.py -h

import sys
import re
import enum
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.constants.constants import alpha

################################################################################
# Distances & penalties
################################################################################
class Distance(enum.Enum):
    edit = 1
    gap_linear = 2
    gap_affine = 3
    gap_affine_2p = 4
    
class Penalties:
    def __init__(self):
      self.distance = Distance.edit
      self.match = 0
      self.mismatch = 0
      self.indel = 0
      self.gap_open1 = 0
      self.gap_extend1 = 0
      self.gap_open2 = 0
      self.gap_extend2 = 0
      
def cigar_compute_score_edit(cigar_vector,penalties,ignore_misms):
  score = 0
  for op in cigar_vector:
    if op[1] in "DI": score += int(op[0])
    if op[1] in "X" and not ignore_misms: score += int(op[0])
  return score

def cigar_compute_score_linear(cigar_vector,penalties,ignore_misms):
  score = 0
  for op in cigar_vector:
    if op[1] == "M": score -= int(op[0]) * penalties.match
    if op[1] == "X" and not ignore_misms: score -= int(op[0]) * penalties.mismatch
    if op[1] in "DI": score -= int(op[0]) * penalties.indel
  return score

def cigar_compute_score_affine(cigar_vector,penalties,ignore_misms):
  score = 0
  for op in cigar_vector:
    if op[1] == "M": score -= int(op[0]) * penalties.match
    if op[1] == "X" and not ignore_misms: score -= int(op[0]) * penalties.mismatch
    if op[1] in "DI": score -= penalties.gap_open1 + int(op[0]) * penalties.gap_extend1
  return score

def cigar_compute_score_affine_2p(cigar_vector,penalties,ignore_misms):
  score = 0
  for op in cigar_vector:
    if op[1] == "M": score -= int(op[0]) * penalties.match
    if op[1] == "X" and not ignore_misms: score -= int(op[0]) * penalties.mismatch
    if op[1] in "DI": 
      score1 = penalties.gap_open1 + int(op[0]) * penalties.gap_extend1
      score2 = penalties.gap_open2 + int(op[0]) * penalties.gap_extend2
      score -= min(score1,score2)
  return score
      
def cigar_compute_score(cigar,penalties,ignore_misms):
  # Parse CIGAR
  #   cigar = "10M3D5M1I10M"
  #   cigar_vector = [('10', 'M'), ('3', 'D'), ('5', 'M'), ('1', 'I'), ('10', 'M')]
  cigar_vector = re.findall(r'(\d+)([MXDI])',cigar)
  # Evaluate CIGAR
  if penalties.distance == Distance.edit:
    return cigar_compute_score_edit(cigar_vector,penalties,ignore_misms)
  elif penalties.distance == Distance.gap_linear:
    return cigar_compute_score_linear(cigar_vector,penalties,ignore_misms)
  elif penalties.distance == Distance.gap_affine:
    return cigar_compute_score_affine(cigar_vector,penalties,ignore_misms)
  elif penalties.distance == Distance.gap_affine_2p:
    return cigar_compute_score_affine_2p(cigar_vector,penalties,ignore_misms)
  
################################################################################
# Alignment stats
################################################################################  
class AlignmentStats:
    def __init__(self):
      self.score_same = 0
      self.score_best1 = 0
      self.score_best2 = 0
      self.scores1 = list()
      self.scores2 = list()

    def __str__(self):
      total_alignments = float(self.score_same + self.score_best1 + self.score_best2)
      str = "[WFACompareAlignments::Stats]\n"
      str += "  => Alignments.common %d (%2.1f %%)\n" % (self.score_same,100.0*float(self.score_same)/total_alignments)
      str += "  => Alignments.best1  %d (%2.1f %%)\n" % (self.score_best1,100.0*float(self.score_best1)/total_alignments)
      str += "  => Alignments.best2  %d (%2.1f %%)"   % (self.score_best2,100.0*float(self.score_best2)/total_alignments)
      return str

    def __repr__(self):
      return __str__(self)

def plot_score_distribution(stats,input_path1,input_path2):
  # Plot
  matplotlib.use('Agg')
  # Draw length histogram
  fig,ax1 = plt.subplots()
  ax1.set_xlabel('Score')
  ax1.set_ylabel('Total Count')
  #   ax1.xaxis.grid(True)
  #   ax1.yaxis.grid(True)
  # Plot score histogram
  range_min = min(min(stats.scores1),min(stats.scores2))
  range_max = max(max(stats.scores1),max(stats.scores2))
  n, bins, patches = ax1.hist(stats.scores1,50,range=[range_min,range_max],color="royalblue",edgecolor='black',alpha=0.5) 
  n, bins, patches = ax1.hist(stats.scores2,50,range=[range_min,range_max],color="darkorange",edgecolor='black',alpha=0.5)
  start, end = ax1.get_xlim()
  ax1.set_xticks(np.arange(start,end,(end-start)/5))
  # Leyend
  handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in ["royalblue","darkorange"]]
  labels= [input_path1,input_path2]
  plt.legend(handles,labels,loc="upper left")
  # Plot
  plt.title("Score distribution for '%s' and '%s'" % (input_path1,input_path2))
  fig.savefig("wfaCmpAlg." + input_path1 + "." + input_path2 + ".png",
              format='png',dpi=100,bbox_inches='tight')
  #plt.show()

################################################################################
# Compare both files
################################################################################
def compare_alignments(input_path1,input_path2,penalties,use_score,ignore_misms,verbose):
  # Alignment stats
  stats = AlignmentStats()  
  # Read both files and compare  
  input_file1 = open(input_path1,"rt")
  input_file2 = open(input_path2,"rt")
  # Read lines
  line_no = 0
  while True:
      # Read lines
      line_no += 1
      line1 = input_file1.readline()
      if not line1: break # End of file
      line2 = input_file2.readline()
      if not line2: 
        print("[WFACompareAlignments] Files unsynch")
        exit(-1)
      # Extract alignments
      fields1 = line1.split()
      fields2 = line2.split()
      if use_score:
        cigar1 = None
        cigar2 = None
        score1 = int(fields1[0]) if len(fields1)<=2 else int(fields1[2])
        score2 = int(fields2[0]) if len(fields2)<=2 else int(fields2[2])
      else:
        cigar1 = fields1[1] if len(fields1)<=2 else fields1[5]
        cigar2 = fields2[1] if len(fields2)<=2 else fields2[5]
        # Evaluate CIGAR's score
        score1 = cigar_compute_score(cigar1,penalties,ignore_misms)
        score2 = cigar_compute_score(cigar2,penalties,ignore_misms)
      # Update stats
      stats.scores1.append(score1)
      stats.scores2.append(score2)
      if score1 == score2:
        stats.score_same += 1
      elif score1 > score2:
        stats.score_best1 += 1
      else:
        stats.score_best2 += 1
      # Verbose
      if verbose and score1 != score2:
        print(">Failed::%d" % line_no)
        print("  s1=%d\t%s" % (score1,cigar1))
        print("  s2=%d\t%s" % (score2,cigar2))
  # Close files
  input_file1.close()
  input_file2.close()
  # return stats
  return stats
  
################################################################################
# Main
################################################################################
# Configure arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i1','--input1',required=True,help='Input file1 (*.alg)')
parser.add_argument('-i2','--input2',required=True,help='Input file2 (*.alg)')
parser.add_argument('-e','--edit',
                    help='Score alignments using edit distance (--edit)')
parser.add_argument('-l','--gap-linear',
                    help='Score alignments using gap-linear distance (--gap-linear=M,X,I)')
parser.add_argument('-g','--gap-affine',
                    help='Score alignments using gap-affine distance (--gap-affine=M,X,O,E)')
parser.add_argument('-G','--gap-affine-2p',
                    help='Score alignments using gap-affine-2p distance (--gap-affine-2p=M,X,O1,E1,O2,E2)')
parser.add_argument('--use-score', action='store_true',default=False,
                    help='Compares the score provided in the file (default=use-cigar)')
parser.add_argument('--ignore-misms', action='store_true',default=False,
                    help='Ignores mismatches when computing the score (default=disable)')
parser.add_argument('-p','--plot', action='store_true',default=False,
                    help='Plots score distribution of both inputs (default=disable)')
parser.add_argument('-v','--verbose', action='store_true',default=False,
                    help='Display comparison differences (default=disable)')
parser.add_argument('-H',action='store_true',dest="human_readable",default=False)

# Parse arguments
args = parser.parse_args()

# Select distance
penalties = Penalties()
if args.edit is not None:
  pass # Default
elif args.gap_linear is not None:
  values = args.gap_linear.split(',')
  penalties.distance = Distance.gap_linear
  penalties.match = int(values[0])
  penalties.mismatch = int(values[1])
  penalties.indel = int(values[2])
elif args.gap_affine is not None:
  values = args.gap_affine.split(',')
  penalties.distance = Distance.gap_affine
  penalties.match = int(values[0])
  penalties.mismatch = int(values[1])
  penalties.gap_open1 = int(values[2])
  penalties.gap_extend1 = int(values[3])
elif args.gap_affine_2p is not None:
  values = args.gap_affine_2p.split(',')
  penalties.distance = Distance.gap_affine_2p
  penalties.match = int(values[0])
  penalties.mismatch = int(values[1])
  penalties.gap_open1 = int(values[2])
  penalties.gap_extend1 = int(values[3])
  penalties.gap_open2 = int(values[4])
  penalties.gap_extend2 = int(values[5])
else:
  print("[WFACompareAlignments] No distance provided. Using edit-distance (default)")
  
# Check penalties
if penalties.match > 0: 
  print("[WFACompareAlignments] Match penalty must be negative or zero")
  exit(-1)
if penalties.mismatch < 0 or \
   penalties.mismatch < 0 or \
   penalties.gap_open1 < 0 or \
   penalties.gap_extend1 < 0 or \
   penalties.gap_open2 < 0 or \
   penalties.gap_extend2 < 0:
  print("[WFACompareAlignments] Penalties must be positive or zero")
  exit(-1)
  
# Compare alignments from both files
stats = compare_alignments(
  args.input1,args.input2,penalties,
  args.use_score,args.ignore_misms,args.verbose)
print(stats) # Display stats

# Plot score distribution
if args.plot:
  plot_score_distribution(stats,args.input1,args.input2)
  
  