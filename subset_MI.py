"""This script takes as input a folder of alignments,
and a folder of trees under the MI criterion. It randomly selects
one alignment and tree per ortholog, and outputs them into a new
folder."""

# Load a tree and link it to an alignment.
from optparse import OptionParser
import os
import random

# get command line input from user using optparse.
parser = OptionParser()
parser.add_option("-i","--inputfolder", help="Input MI folder.",
                    action="store", type="string", dest="inputfolder")
parser.add_option("-o","--outputfolder", help="Output folder.",
                    action="store", type="string", dest="outputfolder")


(options, args) = parser.parse_args()

if not options.inputfolder:   # if slim script name is not given
    parser.error('Input folder not given')
if not options.outputfolder:   # if divergence time is not given
    parser.error('Output folder not given')

os.mkdir(options.outputfolder)

myinput = sorted(os.listdir(options.inputfolder))
ORTHO = None
currentorthos = []
keeping = []
count = 1
for file in myinput:
    if file.endswith('.tre'):
        ORTHO_New = file.split('_')[0]
        if ORTHO_New == ORTHO:
            currentorthos.append(file)
        elif ORTHO_New!= ORTHO:
            if count > 1:
                # randomly select an ortholog from the list
                selectedortho = random.sample(currentorthos, 1)
                keeping.append(selectedortho)
                ORTHO = ORTHO_New
                currentorthos = [file]
            else:
                ORTHO = ORTHO_New
                currentorthos = [file]
                count += 1

for ortholog in keeping:
    os.system('cp %s/%s %s/%s' % (options.inputfolder, ortholog[0], options.outputfolder, ortholog[0]))
   
        

