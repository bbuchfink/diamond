#!/usr/bin/python3
# PROJECT: Wavefront Alignments Algorithms 
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Rescore alignment (*.alg) files
# USAGE: python3 wfa.alg.rescore.py -h

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
      
def cigar_compute_score_edit(cigar_vector,penalties):
  score = 0
  for op in cigar_vector:
    if op[1] in "DI": score += int(op[0])
    if op[1] in "X": score += int(op[0])
  return score

def cigar_compute_score_linear(cigar_vector,penalties):
  score = 0
  for op in cigar_vector:
    if op[1] == "M": score -= int(op[0]) * penalties.match
    if op[1] == "X": score -= int(op[0]) * penalties.mismatch
    if op[1] in "DI": score -= int(op[0]) * penalties.indel
  return score

def cigar_compute_score_affine(cigar_vector,penalties):
  score = 0
  for op in cigar_vector:
    if op[1] == "M": score -= int(op[0]) * penalties.match
    if op[1] == "X": score -= int(op[0]) * penalties.mismatch
    if op[1] in "DI": score -= penalties.gap_open1 + int(op[0]) * penalties.gap_extend1
  return score

def cigar_compute_score_affine_2p(cigar_vector,penalties):
  score = 0
  for op in cigar_vector:
    if op[1] == "M": score -= int(op[0]) * penalties.match
    if op[1] == "X": score -= int(op[0]) * penalties.mismatch
    if op[1] in "DI": 
      score1 = penalties.gap_open1 + int(op[0]) * penalties.gap_extend1
      score2 = penalties.gap_open2 + int(op[0]) * penalties.gap_extend2
      score -= min(score1,score2)
  return score
      
def cigar_compute_score(cigar,penalties):
  # Parse CIGAR
  #   cigar = "10M3D5M1I10M"
  #   cigar_vector = [('10', 'M'), ('3', 'D'), ('5', 'M'), ('1', 'I'), ('10', 'M')]
  cigar_vector = re.findall(r'(\d+)([MXDI])',cigar)
  # Evaluate CIGAR
  if penalties.distance == Distance.edit:
    return cigar_compute_score_edit(cigar_vector,penalties)
  elif penalties.distance == Distance.gap_linear:
    return cigar_compute_score_linear(cigar_vector,penalties)
  elif penalties.distance == Distance.gap_affine:
    return cigar_compute_score_affine(cigar_vector,penalties)
  elif penalties.distance == Distance.gap_affine_2p:
    return cigar_compute_score_affine_2p(cigar_vector,penalties)
  
################################################################################
# Compare both files
################################################################################
def rescore_alignments(input_path,output_path,penalties,verbose):
  # Open files
  input_file = open(input_path,"rt")
  output_file = open(output_path,"wt")
  # Read lines
  while True:
      # Read lines
      line = input_file.readline()
      if not line: break # End of file
      # Extract alignments
      fields = line.split()
      cigar = fields[1] if len(fields)<=2 else fields[5]
      # Evaluate CIGAR's score
      score = cigar_compute_score(cigar,penalties)
      # Print rescore and CIGAR 
      output_file.write("%d\t%s\n" % (score,cigar))
  # Close files
  input_file.close()
  output_file.close()
  
################################################################################
# Main
################################################################################
# Configure arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True,help='Input file (*.alg)')
parser.add_argument('-o','--output',required=True,help='Output file (*.alg)')
parser.add_argument('-e','--edit',
                    help='Score alignments using edit distance (--edit)')
parser.add_argument('-l','--gap-linear',
                    help='Score alignments using gap-linear distance (--gap-linear=M,X,I)')
parser.add_argument('-g','--gap-affine',
                    help='Score alignments using gap-affine distance (--gap-affine=M,X,O,E)')
parser.add_argument('-G','--gap-affine-2p',
                    help='Score alignments using gap-affine-2p distance (--gap-affine=M,X,O1,E1,O2,E2)')
parser.add_argument('-v','--verbose',action='store_true',default=False,
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
  print("[WFAReScoreAlignments] No distance provided. Using edit-distance (default)");
  
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
  
# Rescore alignments
rescore_alignments(args.input,args.output,penalties,args.verbose)

  
  