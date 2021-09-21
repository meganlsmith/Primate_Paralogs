"""This script takes as input a folder with gene trees and a folder with alignments.
The user must supply a minimum number of species for which sequences are required.
Then, the script will prune two-lineage-specific duplicates (and lineage-specific duplicates)
and save trees that contained no duplications
except two-lineage-specific duplicates and lineage-specific duplicates.
Finally, the user must supply an output file for pruned trees. Sequences must be named
Genus_species_OTHERINFO, as the first two underscore-separated strings will be extracted
as the species name."""

# Load a tree and link it to an alignment.
import ete3
from statistics import median
from optparse import OptionParser
#import multiprocessing as mp
import os


# get command line input from user using optparse.
parser = OptionParser()
parser.add_option("-t","--treesfolder", help="File with gene trees in newick format.",
                    action="store", type="string", dest="treesfolder")
parser.add_option("-a","--alignmentfolder", help="File with alignments in fasta format.",
                    action="store", type="string", dest="alignmentfolder")
parser.add_option("-m","--minimum", help="Minimum number of species required in alignment",
                    action="store", type="int", dest="minimum")
parser.add_option("-o","--outputtrees", help="Folder in which to store pruned trees.",
                    action="store", type="string", dest="outputtrees")


(options, args) = parser.parse_args()

if not options.treesfolder: 
    parser.error('Trees folder not given')
if not options.alignmentfolder: 
    parser.error('Alignment location not given')
if not options.minimum:   # if dfe not supported
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
    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)

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

    # create a list to store taxa to drop from tree
    droplist=[]

    # traverse the tree to find duplicates
    for node in t.traverse("postorder"):
        inds = node2labels[node] # get the leaf labels associated with a node
        inds = [x for x in inds if x not in droplist]
        spec = node2species[node] # get the species labels associated with a node
        if len(inds)>1 and len(spec) == 1: # if there is more than one leaf but only one species, then it is a lineage-specific duplicate
            count2=0 # counter so we know if we're on the first one
            for ind in inds: # for each individual sampled
                indseq = ((t&ind).sequence) # get the sequence
                indseqnogap = indseq.replace('-','') # get length excluding missing data/ gaps
                seqdif = abs(len(indseqnogap)-medianlength) # calculate difference in sequence length and median sequence length in the alignment
                if count2 == 0: # if it's the first one we're looking at
                    mindif=seqdif # get differences as minimum current difference
                    tokeep=ind # mark individual as current 'to keep'
                    count2+=1 # increment counter
                else: # if it's not the first one we're evaluating
                    if seqdif < mindif: # if this sequence is closer to the median length than the previous 'to keep' sequence
                        tokeep=ind # then make this the new 'to keep' sequence
                        mindif=seqdif # and update the minimum difference
            todrop=[x for x in inds if x not in tokeep] # all individuals except to keep need to be dropped
            for itemtodrop in todrop:
                nodedrop=t.search_nodes(name=itemtodrop)[0]
                nodedrop.delete(preserve_branch_length=True)
            droplist=droplist + todrop # add these to our list of leaves to drop.


     # drop two lineage duplicates
    droplist2 = []
    for node in t.traverse("postorder"):
        spec = node2species[node] # get the species labels associated with a node
        indivlist = list(node.get_leaf_names())
        indivlist = [x for x in indivlist if x not in droplist2]
        if len(indivlist)>2 and len(spec) == 2: # if there is a duplication event specific to two species:
            count3=1
            speclist = list(spec)
            speciesA = speclist[0]
            speciesB = speclist[1]
            for inda in range(0, len(indivlist)-1):
                for indb in range(inda+1, len(indivlist)):
                    leafa = t&indivlist[inda]
                    leafb = t&indivlist[indb]
                    distance = t.get_distance(leafa, leafb)
                    if (speciesA in indivlist[inda] and speciesB in indivlist[indb]) or (speciesB in indivlist[inda] and speciesA in indivlist[indb]):
                        if count3 == 1:
                            mindist = distance
                            chosen = [indivlist[inda], indivlist[indb]]
                            count3+=1
                        else:
                            if distance < mindist:
                                mindist = distance
                                chosen = [indivlist[inda], indivlist[indb]]
            todrop2=[x for x in indivlist if x not in chosen]
            for item in range(0, len(todrop2)):
                nodedrop = t.search_nodes(name=todrop2[item])[0]
                nodedrop.delete(preserve_branch_length=True)
            droplist2 = droplist2+todrop2

    node2labels = t.get_cached_content(store_attr="name")
    leaves = t.get_leaves()
    species = t.get_species()
    if len(leaves) == len(species) and len(species)>=options.minimum:
        alignmentname = alignment.split('/')[-1]
        newalignmentname = options.outputtrees+'/'+alignmentname
        treefiledest = options.outputtrees+'/'+alignmentname.split('.fa')[0]+'.tre'
        t.unroot()
        t.write(outfile=treefiledest)
        allleaves = []
        for leaf in t.iter_leaves():
            allleaves.append(str(leaf).split('--')[1].strip())
        downsamplealignment(allleaves, alignment, newalignmentname)


myoutputtrees = os.listdir(options.outputtrees)
with open('%s.tre.tmp' % options.outputtrees, 'w') as outfile:
    for fname in myoutputtrees:
        if fname.endswith('.tre'):
            with open(os.path.join(options.outputtrees,fname)) as infile:
                outfile.write(infile.read())
                outfile.write('\n')


os.system('mv %s.tre.tmp %s.tre' % (options.outputtrees, options.outputtrees))
