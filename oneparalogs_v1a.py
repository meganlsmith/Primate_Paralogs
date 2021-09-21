"""This script takes as input a folder with gene trees and a folder with alignments.
The user must supply a minimum number of species for which sequences are required.
Then, the script will keep a single sample per species.
Finally, the user must supply an output file for pruned trees (default ./OneParalogs_trees).
 Sequences must be named
Genus_species_OTHERINFO, as the first two underscore-separated strings will be extracted
as the species name."""

# Load a tree and link it to an alignment.
import ete3
from statistics import median
from optparse import OptionParser
#import multiprocessing as mp
import os
import random


# get command line input from user using optparse.
parser = OptionParser()
parser.add_option("-t","--treesfolder", help="File with gene trees in newick format.",
                    action="store", type="string", dest="treesfolder")
parser.add_option("-a","--alignmentfolder", help="File with alignments in fasta format.",
                    action="store", type="string", dest="alignmentfolder")
parser.add_option("-m","--minimum", help="Minimum number of species required in alignment",
                    action="store", type="int", dest="minimum")
parser.add_option("-o","--outputtrees", help="Folder in which to store pruned trees.",
                    default="OneParalogs_trees",
                    action="store", type="string", dest="outputtrees")


(options, args) = parser.parse_args()

if not options.treesfolder: 
    parser.error('Trees folder not given')
if not options.alignmentfolder:
    parser.error('Alignment location not given')
if not options.minimum:
    parser.error('Minimum number of species not supplied.')


def get_seq_length(tree):
    """This function will get the median length of all sequences for a tree."""
    lengths = []
    for leaf in tree.iter_leaves():
        length = len(leaf.sequence.replace('-',''))
        lengths.append(length)
    return(median(lengths))

def get_species_name(node_name_string):
    """this function will get species names for the leaves of a tree."""
    # Species code is the first part of leaf name (separated by an
    #  underscore character)
    spcode = node_name_string.split("_")[0:2]
    spcode = spcode[0]+'_'+spcode[1]
    # We could even translate the code to complete names
    return spcode

def downsamplealignment(keeplist, alignment, newalignmentname):
    """this function will take a list of individuals to keep and remove all other
    samples from the alignment, creating a new alignment in the output folder."""
    newalignmentwrite = open(newalignmentname, 'w')
    with open(alignment, 'r') as alignmentfile:
        WRITENOW=False
        for line in alignmentfile.readlines():
            if WRITENOW==False and '>' in line and line.split('>')[1].strip() in keeplist:
                newalignmentwrite.write(line)
                WRITENOW=True
            elif WRITENOW==True:
                newalignmentwrite.write(line)
                WRITENOW=False



# make the directory four output gene tree
os.mkdir(options.outputtrees)

# get list of trees
treestoparse = os.listdir(options.treesfolder)

for tree in treestoparse: # for each tree
    full_path = os.path.join(options.treesfolder, tree)

    #read the tree
    t = ete3.PhyloTree(full_path)

    # find the alignment and associate it with the tree.
    alignmentname = tree.split('.treefile')[0]
    alignment = os.path.join(options.alignmentfolder, alignmentname)
    t.link_to_alignment(alignment=alignment, alg_format="fasta")

    # assign species names to tips
    t.set_species_naming_function(get_species_name)

    # get the median sequence length
    medianlength = get_seq_length(t)

    # get the leaf node to leaf labels
    node2labels = t.get_cached_content(store_attr="name")

    # get the node to species labels
    node2species = t.get_cached_content(store_attr="species")
    
    # find cases where there are multiple sequences for a species and downsample the alignment
    species = t.get_species() # get names of species
    leaves = t.get_leaves() # get sample names
    keep = [] # start list of sequences to discard
    # loop through species and randomly sample when necessary
    for thisspecies in species:
        samplesforspecies = []
        for thisleaf in leaves:
            sample = str(thisleaf).split('--')[1]
            if thisspecies in sample:
                samplesforspecies.append(sample)
        if len(samplesforspecies) > 1:
            random.sample(samplesforspecies, 1)
            keep.append(random.sample(samplesforspecies, 1)[0])
        else:
            keep.append(samplesforspecies[0])
    
            
    # create a list to store taxa to drop from tree
    tokeep=[]
        
    # traverse the tree to find duplicates
    for node in t.traverse("postorder"):
        inds = node2labels[node] # get the leaf labels associated with a node
        spec = node2species[node] # get the species labels associated with a node
        for individual in inds:
            if individual in keep: # if there is more than one leaf but only one species, then it is a lineage-specific duplicate
                tokeep.append(individual) # add these to our list of leaves to drop.
    
    t.prune(tokeep, preserve_branch_length = True)

#    check to see if the tree contained only enough species
    node2labels = t.get_cached_content(store_attr="name")
    leaves = t.get_leaves()
    species = t.get_species()
    if len(species)>=options.minimum:
        alignmentname = alignment.split('/')[-1]
        newalignmentname = options.outputtrees+'/'+alignmentname
        treefiledest = options.outputtrees+'/'+alignmentname.split('.fa')[0]+'.tre'
        t.write(outfile=treefiledest)
        downsamplealignment(keep, alignment, newalignmentname)            


myoutputtrees = os.listdir(options.outputtrees)
with open('%s.tre.tmp' % options.outputtrees, 'w') as outfile:
    for fname in myoutputtrees:
        if fname.endswith('.tre'):
            with open(os.path.join(options.outputtrees,fname)) as infile:
                outfile.write(infile.read())
                outfile.write('\n')


os.system('mv %s.tre.tmp %s.tre' % (options.outputtrees, options.outputtrees))

