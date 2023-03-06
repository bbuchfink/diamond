#!/usr/bin/python3
# PROJECT: Wavefront Alignments Algorithms 
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Plot WFA alignment matrices
# USAGE: python3 wfa.plot.py -h

import sys
import copy
import glob
import os.path
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import seaborn as sns

################################################################################
# WFA-file parsing
################################################################################
def wfa_parse_file(filename):
  # Log
  print('[Parsing]',end='',flush=True)
  # Parse
  wfa_info = {}
  with open(filename) as fp:
    line = fp.readline() # Read line
    while line:
      # Parse data
      if line[0] == '#':
        fields = line.split();
        if fields[1] == "Heatmap":
          matrix_name = fields[2]
          matrix_list = []
          matrix_list.append(fp.readline().strip().split(',')) 
          line = fp.readline()
          while line and line[0] != '#':
            matrix_list.append(line.strip().split(',')) 
            line = fp.readline()
          M = np.matrix(matrix_list).astype(int)
          M = np.where(M==-1,np.nan,M)
          wfa_info[matrix_name] = M
          # Parse next line
          if not line: break
          else: continue
        elif fields[1] == "List":  
          if len(fields) <= 3: 
            wfa_info[fields[2]] = None
          else: 
            CIGAR = np.matrix(fields[3].replace(","," ")).astype(int)
            wfa_info[fields[2]] = CIGAR
        else:    
          wfa_info[fields[1]] = fields[2]
      # Read next line
      line = fp.readline()
  # Process some values
  wfa_info["WFAMode"] = wfa_info["WFAMode"][1:-1] 
  wfa_info["Distance"] = wfa_info["Penalties"][1:-1].split(',')[0]
  # Return
  return wfa_info

################################################################################
# WFA-Heatmap plotting
################################################################################
def wfa_plot_ticks(ax,data,detailed):
  # Compute dimensions
  ylen = data.shape[0]
  xlen = data.shape[1]
  plen = int(wfa_info["PatternLength"])
  tlen = int(wfa_info["TextLength"])
  # Detailed mode
  if detailed:
    # Parameters
    pattern = wfa_info["Pattern"]
    text = wfa_info["Text"]
    # Hacky font size estimation
    label_size = 80.0/float(max(xlen,ylen))
    
    # Set Y-axis (major)
    y_ticks_loc = np.arange(-0.5,ylen,1)
    ax.yaxis.set_minor_locator(mticker.FixedLocator(y_ticks_loc))
    ax.set_yticks(y_ticks_loc)
    ax.set_yticklabels([])
    # Set Y-axis (minor)
    ax.set_yticks(np.arange(0,ylen,1),minor=True)
    ax.set_yticklabels([y for y in pattern],fontsize=label_size,fontweight='bold',minor=True)
    
    # Set X-axis (major)
    x_ticks_loc = np.arange(-0.5,xlen,1)
    ax.xaxis.set_minor_locator(mticker.FixedLocator(x_ticks_loc))
    ax.set_xticks(x_ticks_loc)
    ax.set_xticklabels([])
    # Set X-axis (minor)
    ax.set_xticks(np.arange(0,xlen,1),minor=True)
    ax.set_xticklabels([x for x in text],fontsize=label_size,fontweight='bold',minor=True)
    # Set X on top
    ax.xaxis.tick_top()
    
    # Hide minor dashes
    ax.tick_params(axis='both',which='minor',length=0) 
    # Gridlines (based on major)
    ax.grid(which='major',color="black",linestyle='-',linewidth=0.25)
  else:
    # X-ticks
    x_ticks = [i for i in range(0,xlen-1,xlen//5)]
    if x_ticks[-1]+5 < xlen-1: x_ticks.append(xlen-1)  
    else: x_ticks[-1] = xlen-1
    x_ratio = tlen/xlen
    x_ticks_labels = [int(i*x_ratio) for i in x_ticks]
    x_ticks_labels[-1] = tlen-1
    # Y-ticks
    y_ticks = [i for i in range(0,ylen-1,ylen//5)]
    if y_ticks[-1]+5 < ylen-1: y_ticks.append(ylen-1)  
    else: y_ticks[-1] = ylen-1
    y_ratio = plen/ylen
    y_ticks_labels = [int(i*y_ratio) for i in y_ticks]
    y_ticks_labels[-1] = plen-1
    # Set ticks
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticks_labels,fontsize=6)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks_labels,fontsize=6)
    # Set grid
    ax.grid(which='major',color="black",linestyle='-',linewidth=0.25)

def wfa_plot_wavefront(wfa_info,title,ax,data,detailed,
                       xlabel=False,ylabel=False):
  # Compute dimensions
  ylen = data.shape[0]
  xlen = data.shape[1]
  # Title
  ax.set_title(title,fontsize=10) 
  # Set ticks & grid
  wfa_plot_ticks(ax,data,detailed)
  # Set labels
  if xlabel: ax.set_xlabel('Text',fontsize=8)
  if ylabel: ax.set_ylabel('Pattern',fontsize=8)
  # Create colorbar
  #cmap = LinearSegmentedColormap.from_list('rg',['tab:green','yellow','tab:red'],N=256)
  cmap = copy.copy(plt.cm.jet)
  cmap.set_bad('whitesmoke')
  # Heatmap
  im = ax.imshow(data,cmap=cmap)
  # Detailed mode (Loop over data and create text annotations)
  if detailed:
    fontsize = 80.0/float(max(xlen,ylen))
    for i in range(len(data)):
      for j in range(len(data[i])):
        if data[i,j]>=0:
          ax.text(j,i,int(data[i,j]),ha="center",va="center",color='w',fontsize=fontsize,fontweight='bold')
  # Return
  return im

def wfa_plot_cigar(ax,wfa_info):
  # Parameters
  plen = int(wfa_info["PatternLength"])
  tlen = int(wfa_info["TextLength"])
  ylen = wfa_info["M"].shape[0]
  xlen = wfa_info["M"].shape[1]
  marker_size = 80.0/float(max(xlen,ylen))
  # Plot CIGAR series
  def wfa_plot_cigar_scaled(cigar,x_ratio,y_ratio,marker,color,label):
    x = np.floor(cigar[:,0].astype('float64') * x_ratio).astype('int') if x_ratio < 1.0 else cigar[:,0]
    y = np.floor(cigar[:,1].astype('float64') * y_ratio).astype('int') if y_ratio < 1.0 else cigar[:,1]
    ax.scatter([x],[y],marker=marker,color=color,s=marker_size,linewidths=0,label=label)
  # Compute dims
  x_ratio = xlen/tlen  
  y_ratio = ylen/plen
  # Fetch CIGAR
  cigar_m = wfa_info["CIGAR-M"]
  cigar_x = wfa_info["CIGAR-X"]
  cigar_i = wfa_info["CIGAR-I"]
  cigar_d = wfa_info["CIGAR-D"]
  # Plot CIGAR
  if (cigar_m is not None): wfa_plot_cigar_scaled(cigar_m,x_ratio,y_ratio,',','limegreen','match')
  if (cigar_x is not None): wfa_plot_cigar_scaled(cigar_x,x_ratio,y_ratio,'o','red','misms')
  if (cigar_i is not None): wfa_plot_cigar_scaled(cigar_i,x_ratio,y_ratio,'>','orange','ins')
  if (cigar_d is not None): wfa_plot_cigar_scaled(cigar_d,x_ratio,y_ratio,'v','blue','del')
  ax.legend(loc="upper right",prop={'size': 3}) 

def wfa_plot_xtra(wfa_info,title,ax,data,detailed, 
                  xlabel=False,ylabel=False): 
  # Title
  ax.set_title(title,fontsize=10) 
  # Set ticks & grid
  wfa_plot_ticks(ax,data,detailed)
  # Set labels
  if xlabel: ax.set_xlabel('Text',fontsize=8)
  if ylabel: ax.set_ylabel('Pattern',fontsize=8)
  # Create colorbar
  cmap_ext = colors.ListedColormap(['dodgerblue','darkred'])
  cmap_ext_bounds = [10,20]
  cmap_ext_norm = colors.BoundaryNorm(cmap_ext_bounds,cmap_ext.N)
  cmap_ext.set_bad('whitesmoke')
  # Heatmap
  im = ax.imshow(data,cmap=cmap_ext)
  # Plot CIGAR
  if 'A' in wfa_info["WFAMode"]: wfa_plot_cigar(ax,wfa_info)
  # Return
  return im
 
def wfa_plot(filename,wfa_info,dpi,mode,detailed):
  # Log
  print('[Plotting]',end='',flush=True)
  # Parameters
  plen = int(wfa_info["PatternLength"])
  tlen = int(wfa_info["TextLength"])
  ylen = wfa_info["M"].shape[0]
  xlen = wfa_info["M"].shape[1]
  # Create plot
  compact = (mode=="compact");
  extended = (mode=="extended"); 
  if compact: 
    fig, ax1 = plt.subplots(nrows=1,ncols=1,dpi=dpi,sharex=True)
    im1 = wfa_plot_wavefront(wfa_info,'M-Wavefront',ax1,wfa_info["M"],detailed,xlabel=True,ylabel=True)
  elif extended or wfa_info["Distance"]=="Edit":
    fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,dpi=dpi,sharex=True)
    im1 = wfa_plot_wavefront(wfa_info,'M-Wavefront',ax1,wfa_info["M"],detailed,xlabel=True,ylabel=True)
    im3 = wfa_plot_xtra(wfa_info,'CIGAR',ax2,wfa_info["Extend"],detailed,xlabel=True)
  elif wfa_info["Distance"]=="GapLineal" or wfa_info["Distance"]=="GapAffine":
    fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,dpi=dpi,sharex=True)
    im1 = wfa_plot_wavefront(wfa_info,'M-Wavefront',ax1,wfa_info["M"],detailed,xlabel=True,ylabel=True)
    im2 = wfa_plot_wavefront(wfa_info,'I1-Wavefront',ax2,wfa_info["I1"],detailed,xlabel=True)
    im3 = wfa_plot_wavefront(wfa_info,'D1-Wavefront',ax3,wfa_info["D1"],detailed,xlabel=True)
  elif wfa_info["Distance"]=="GapAffine2p":    
    fig, ((ax1,ax2,ax5),(ax3,ax4,ax6)) = plt.subplots(nrows=2,ncols=3,dpi=dpi,sharex=True)
    im1 = wfa_plot_wavefront(wfa_info,'M-Wavefront',ax1,wfa_info["M"],detailed,xlabel=True,ylabel=True)
    im2 = wfa_plot_wavefront(wfa_info,'I1-Wavefront',ax2,wfa_info["I1"],detailed,xlabel=True)
    im4 = wfa_plot_wavefront(wfa_info,'D1-Wavefront',ax4,wfa_info["D1"],detailed,xlabel=True)
    im5 = wfa_plot_wavefront(wfa_info,'I2-Wavefront',ax5,wfa_info["I2"],detailed,ylabel=True)
    im6 = wfa_plot_wavefront(wfa_info,'D2-Wavefront',ax6,wfa_info["D2"],detailed)
    im3 = wfa_plot_xtra(wfa_info,'CIGAR',ax3,wfa_info["Extend"],detailed)
  # Color bar
  if compact:
    p0 = ax1.get_position().get_points().flatten()
    p1 = p0
    pass
  elif extended or wfa_info["Distance"]=="Edit":
    p0 = ax1.get_position().get_points().flatten()
    p1 = ax2.get_position().get_points().flatten()
    pass
  elif wfa_info["Distance"]=="GapLineal" or wfa_info["Distance"]=="GapAffine":
    p0 = ax1.get_position().get_points().flatten()
    p1 = ax3.get_position().get_points().flatten()
  elif wfa_info["Distance"]=="GapAffine2p":    
    p0 = ax3.get_position().get_points().flatten()
    p1 = ax6.get_position().get_points().flatten()
  ax_cbar = fig.add_axes([p0[0],0,p1[2]-p0[0],0.025])
  ax_cbar.tick_params(labelsize=6) 
  plt.colorbar(im1,cax=ax_cbar,orientation='horizontal')
  # Title
  file = os.path.basename(filename).replace('.plot','')
  title = "WFA-Plot(%s) %s[%s]" % (file,wfa_info["Penalties"],wfa_info["WFAMode"])
  fig.suptitle(title,fontsize=10)
  plt.subplots_adjust(top=0.85)
  # Plot
  plt.savefig(filename.replace('.plot','.png'),bbox_inches='tight')

################################################################################
# Main
################################################################################
# Configure arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',action='store',help='Input file (*.plot)')
parser.add_argument('--dpi',type=int,action='store',default=1000,help='Plot resolution (default=1000)') # More than 2000 is hard to handle
parser.add_argument('--mode',action='store',default='compact',help='Plot mode in {compact,extended,full}')
parser.add_argument('-d','--detailed',action='store_true',default=False,help='Plot score values and sequences')
parser.add_argument('-H',action='store_true',dest="human_readable",default=False)

# Parse arguments
args = parser.parse_args()

# Open input
if args.input:
  input_files = [args.input]
else:
  input_files = glob.glob("*.plot")
  print('[WFA2png] Searching all *.plot ( Found %d file%c)' % (len(input_files),'s' if len(input_files)>1 else ' '))
  
# Plot each WFA file
print('[WFA2png] Plotting at %d dpi' % (args.dpi))
idx = 0
for filename in input_files:
  print('[WFA2png] [#%d] Generating \'%s\' ' % (idx,filename),end='',flush=True)
  wfa_info = wfa_parse_file(filename)
  wfa_plot(filename,wfa_info,args.dpi,args.mode,args.detailed)
  print('[Done!]',)
  idx += 1


  
  
