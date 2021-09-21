import random
import numpy as np
from optparse import OptionParser
import scipy.stats
import pandas

# get command line input from user using optparse.
parser = OptionParser()
parser.add_option("-i","--inputfile", help="gCF output from IQTree.",
                    action="store", type="string", dest="inputfile")
parser.add_option("-o","--outputfile", help="CSV file for results.",
                    action="store", type="string", dest="outputfile")
(options, args) = parser.parse_args()

if not options.inputfile:   # if name of gcf file not given
    parser.error('Inputfile file not given.')
df = pandas.DataFrame(columns = ['delta', 'stdev', 'p_values', 'branch'])

infile = open(options.inputfile, 'r').readlines()
countofcomparisons=0
for line in infile:
    resampled_gCFN = 0
    resampled_gDF1N = 0
    resampled_gDF2N = 0
    resampled_gDFPN = 0
    if not line.startswith('#'):
        if line.startswith('ID'):
            header = line.split()
        else:
            splitline=line.split()
            if splitline[0] == '31':
                continue
            else:        
                """Get the information from the table."""   
                branch=float(splitline[0])
                gCF=float(splitline[1])
                gCF_N=int(splitline[2])
                gDF1=float(splitline[3])
                gDF1_N=int(splitline[4])
                gDF2=float(splitline[5])
                gDF2_N=int(splitline[6])
                gDFP=float(splitline[7])
                gDFP_N=int(splitline[8])
                gN=int(splitline[9])
                treeDF1=splitline[10]
                treeDF2=splitline[11]
                """Calculate delta"""
                
                if (gDF1+gDF2>5):
                    delta = abs(gDF1_N- gDF2_N)/(gDF1_N+gDF2_N)

                    countofcomparisons+=1
                    """Do bootstrap replicates."""
                    deltaresample = []
                    for rep in range(1, 2001):
                        resampled_gCFN = 0
                        resampled_gDF1N = 0
                        resampled_gDF2N = 0
                        resampled_gDFPN = 0
                        # for each replicate draw a list of gene trees to sample with replacement
                        tosample = np.random.choice(gN, gN)
                        for sample in tosample:
                            if sample <= gCF_N:
                                resampled_gCFN+=1
                            elif sample <= gDF1_N+gCF_N:
                                resampled_gDF1N+=1
                            elif sample <= gDF2_N+gCF_N+gDF1_N:
                                resampled_gDF2N+=1
                            elif sample <= gDFP_N+gCF_N+gDF1_N+gDF2_N:
                                resampled_gDFPN+=1
                        if resampled_gDF1N+resampled_gDF2N==0:
                            thisdelta=np.nan
                        else:
                            thisdelta = abs(resampled_gDF1N- resampled_gDF2N)/(resampled_gDF1N+resampled_gDF2N)
                        deltaresample.append(thisdelta)
                    stdev = np.nanstd(deltaresample)
                    mean = np.nanmean(deltaresample)
                    Z=(delta-0)/stdev
                    p_values = scipy.stats.norm.sf(abs(Z))*2 #twosided
                    df = df.append({'delta' : delta, 'stdev' : stdev, 'p_values' : p_values, 'branch' : branch}, ignore_index=True)
duncansidol = 1 - (1-0.05)**(1/countofcomparisons)
print('Duncan sidol p value is %r ' %duncansidol)

df.to_csv(path_or_buf=options.outputfile)
#                
                

