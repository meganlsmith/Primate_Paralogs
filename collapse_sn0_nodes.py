"""This script will take as input a folder of gene trees.
It will loop through the gene trees. For each gene tree, it will consider each node. For each node,
It will match it with a node in the species tree. Then, it will check the number of informative sites
for that node. If the number of sites is zero, then the node will be collapsed. We'll put these into a new folder.

input:  gene tree folder
        alignment folder
actions: collapse nodes for which there are no informative sites
output: gene tree folder with zero site nodes collapsed.


NOTE: PLEASE PROVIDE FULL PATHS BECAUSE COMPUTATIONS WILL HAPPEN IN SPECIFIED OUTPUT DIRECTORY.
"""

from optparse import OptionParser
from collections import OrderedDict
import ete3
import os
import pandas as pd

# get command line input from user using optparse.
parser = OptionParser()
parser.add_option("-i","--genetreefolder", help="Folder with all gene trees.",
                    action="store", type="string", dest="genetreefolder")
parser.add_option("-a","--alignmentfolder", help="Folder with all alignments.",
                    action="store", type="string", dest="alignmentfolder")
parser.add_option("-o","--outputfolder", help="Folder for output.",
                    action="store", type="string", dest="outputfolder")
parser.add_option("-q","--iqtreepath", help="Path to IQTree.",
                    action="store", type="string", dest="iqtreepath")
parser.add_option("-c","--cores", help="Number of cores.",
                    action="store", type="string", dest="cores")
(options, args) = parser.parse_args()

# make output Folder
os.mkdir(options.outputfolder)
os.chdir(options.outputfolder)

# loop through gene trees
genetreelist = os.listdir(options.genetreefolder)
genetreelist = [x for x in genetreelist if x.endswith('tre')]

for genetree in genetreelist:
    genename = genetree.split('.tre')[0]
    if '1to1ortho' in genename:
        genename = genename.split('_1to1ortho')[0]+'.Noambig'

    # get path to tree and alignment
    treename = options.genetreefolder+'/'+genetree
    alignment = options.alignmentfolder + genename + '.fa'

    # midpoint root tree and write temporary file
    otree = ete3.PhyloTree(treename)
    # Calculate the midpoint node
    R = otree.get_midpoint_outgroup()
    
    # and set it as tree outgroup
    if len(R) != len(otree):
        # and set it as tree outgroup
        otree.set_outgroup(R)
    otree.write(outfile='temp.tre')

    # run IQTree
    os.system('%s -t %s  -s %s  --scf 100 --prefix %s -T %s' % (options.iqtreepath, 'temp.tre', alignment, genename, options.cores))

    # read original tree and tree with branch labels
    ltree = ete3.Tree(genename+'.cf.branch', format=1)

    # create dictionry of labels and taxa
    labeldict = OrderedDict()
    # get the leaf node to leaf labels
    node2labels = ltree.get_cached_content(store_attr="name")
    # traverse nodes and create dictionary associating node IDs to a list of species below that node
    for node in ltree.traverse("postorder"):
        # get node ID
        name = node.get_cached_content(store_attr="name")
        names = list(name.keys())
        namesstring = str(names)
        nodeid = namesstring.split(', Tree node ')[-1].split()[0].strip("'")
        inds = node2labels[node]
        newinds = []
        for x in inds:
            new_x = x.replace("__", "@_")
            newinds.append(new_x)
        if nodeid != '[Tree':
            labeldict.update({nodeid: newinds})

    # read in stat table
    stats = pd.read_csv(genename + '.cf.stat', skiprows=[0,1,2,3,4,5,6,7,8,9,10,11,12,13], sep = '\t')

    # get rows where sN == 0 and record node IDs to be collapsed
    zero = stats.loc[stats['sN'] == 0]
    nodestocollapse = zero['ID'].tolist()

    # traverse the tree associate node labels with taxa using dictonary
    for node in otree.traverse("postorder"):

        # get node ID
        name = node.get_cached_content(store_attr="name")
        node2labels2 = node.get_cached_content(store_attr="name")
        if len(node2labels2) > 1:

            # get list of individuals
            inds = list(node2labels2[node])
            # iterate over dictionary and find FIRST item that includes all of those species
            previous = False
            for key in labeldict:
                indsforkey = list(labeldict[key])
                listdiff = list(set(inds) - set(indsforkey))
                if indsforkey == inds:
                    match = key
                    previous = True
                    exact = True

                elif len(listdiff) == 0 and previous == False:
                    match = key
                    previous = True
                    exact = False

            # check to see if match is in tocollapselist
            if int(match) in nodestocollapse:
                node.delete(preserve_branch_length=True)
    otree.unroot()
    outname = genetree
    otree.write(outfile=outname)
    os.system('rm *cf*')
    os.system('rm *.log*')
    os.system('rm temp.tre')
