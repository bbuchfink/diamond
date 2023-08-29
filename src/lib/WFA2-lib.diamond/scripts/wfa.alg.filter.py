################################################################################
# PROJECT: SEQAlignment tools
# NAME: Alignment filtering tool
# AUTHOR: Santiago Marco Sola <santiagomsola@gmail.com>
################################################################################
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

#
# Main
#
# Configure arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True,
                    help="Input alignment file (*.alg)")
parser.add_argument('-o','--output',default=None,
                    help="Output sequence file (*.seq)")
parser.add_argument('-d','--max-length-diff',type=float,default=-1.0,
                    help="Filter sequences by length difference percentage")
parser.add_argument('-e','--max-error',type=float,default=-1.0, 
                    help="Filter sequences by alignment error")
# Parse arguments
args = parser.parse_args()
max_ldiff = args.max_length_diff
max_error = args.max_error

# Length/Error distributions
dist_length = []
dist_error = []

# Read/Write files
input_file = open(args.input, "rt")
output_file = open(args.output, "wt") if args.output != None else None
while True:
  # Read line
  line = input_file.readline()
  if not line: break
  # Split
  fields = line.split()
  # Check relative length difference
  plen = int(fields[0])
  tlen = int(fields[1])
  if max_ldiff >= 0.0 and abs(plen-tlen) > plen*max_ldiff: continue
  # Check error
  error = - int(fields[2])
  if max_error >= 0.0 and error > plen*max_error: continue
  # Print output
  if output_file != None:
    output_file.write(">%s\n" % fields[3])
    output_file.write("<%s\n" % fields[4])
  # Store stats
  dist_length.append(plen/1000)
  dist_length.append(tlen/1000)
  dist_error.append(error)
  
# Close files
input_file.close()
if output_file != None: output_file.close()

# Plot
matplotlib.use('Agg')

# Draw length histogram
ax1_color = "royalblue"
fig,ax1 = plt.subplots()
ax1.set_xlabel('Length (Kbp)')
ax1.xaxis.grid(True)
ax1.yaxis.grid(True)
ax1.set_ylabel('Sequence Count',color=ax1_color)

# Draw length histogram (cumulative)
ax3_color = "darkorange"
ax3 = ax1.twinx()
values, base, _ = ax1.hist(dist_length,bins=100,color=ax1_color)
values = np.append(values,0)
ax3.plot(base,np.cumsum(values)/np.cumsum(values)[-1],color=ax3_color, 
         marker='o',linestyle='-',markersize=1)
ax3.set_ylabel('Sequence Count (%)',color=ax3_color)

# Draw error histogram
ax2_color = "firebrick"
y, bin_length = np.histogram(dist_error,bins=100)
bin_centers = 0.5 * (bin_length[1:] + bin_length[:-1])
ax2 = ax1.twiny()
ax2.set_xlabel('Nominal error',color=ax2_color)
ax2.plot(bin_centers,y,'-',c=ax2_color)

# Leyend
handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in [ax1_color,ax2_color,ax3_color]]
labels= ["Seq. Length Histogram","Alignment Error Distribution","Seq. Length Cumulative Histogram"]
plt.legend(handles,labels,loc="upper left")

# Plot
plt.title("Sequence statistics for '%s'" % args.input)
fig.savefig(args.input + ".png",
            format='png',
            dpi=100,
            bbox_inches='tight')
#plt.show()
