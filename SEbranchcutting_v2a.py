"""This script takes as input a folder with gene trees and a folder with alignments.
The user must supply a minimum number of species for which sequences are required.
The script will trim lineage-specific duplicates, two species duplicates, and perform the new
branch-cutting. It will then output any tree where there is one sample per species.
Finally, the user must supply an output file for pruned trees. Sequences must be named
Genus_species_OTHERINFO, as the first two underscore-separated strings will be extracted
as the species name."""

# Load a tree and link it to an alignment.
import ete3
from statistics import median
from optparse import OptionParser
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


    # traverse the tree to find duplicates
    for node in t.traverse("postorder"):
        count2 = 0
        inds = list(node.get_leaf_names())
        spec = node2species[node] # get the species labels associated with a node
        if len(inds) > len(spec):
            if len(inds) > 1 and len(spec) == 1:
                for ind in inds:
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
                for i in todrop:
                    nodedrop =t.search_nodes(name=i)[0]
                    nodedrop.delete(preserve_branch_length=True)

    for node in t.traverse("postorder"):
        inds = list(node.get_leaf_names())
        spec = node2species[node] # get the species labels associated with a node
        if len(inds) > len(spec):

            if len(inds)>2 and len(spec) == 2: # if there is more than two leaves but two species, then it is a two lineage-specific duplicate
                count3=1
                speclist = list(spec)
                speciesA = speclist[0]
                speciesB = speclist[1]
                for inda in range(0, len(inds)-1):
                    for indb in range(inda+1, len(inds)):
                        leafa = t&inds[inda]
                        leafb = t&inds[indb]
                        distance = t.get_distance(leafa, leafb)
                        if (speciesA in inds[inda] and speciesB in inds[indb]) or (speciesB in inds[inda] and speciesA in inds[indb]):
                            if count3 == 1:
                                mindist = distance
                                chosen = [inds[inda], inds[indb]]
                                count3+=1
                            else:
                                if distance < mindist:
                                    mindist = distance
                                    chosen = [inds[inda], inds[indb]]
                todrop=[x for x in inds if x not in chosen]
                for i in todrop:
                    nodedrop =t.search_nodes(name=i)[0]
                    nodedrop.delete(preserve_branch_length=True)

    count5=1
    for node in t.traverse("postorder"):
        inds = list(node.get_leaf_names())
        spec = node.iter_species()
        speclist = []
        for specname in spec:
            speclist.append(specname)
        spec = speclist
        if len(inds) > len(spec):
            count4=1
            childrentostore = []
            found = False
            for child in node.children:
                childinds = list(child.get_leaf_names())
                childspecies = node2species[child]
                if count4 == 1:
                    maxspecies = len(childspecies)
                    childtokeep = child
                    childrentostore.append(child)
                    count4+=1
                    found=True
                else:
                    childrentostore.append(child)
                    if len(childspecies) > maxspecies:
                        maxspecies = len(childspecies)
                        childtokeep = child
            for child in childrentostore:
                childinds = list(child.get_leaf_names())
                childspecies = node2species[child]
                if child != childtokeep:
                    s = t.copy()
                    s.prune(childinds, preserve_branch_length=True)
                    if len(childinds) == len(childspecies) and len(childspecies) > options.minimum:
                        tonametree = options.outputtrees + '/'+ alignmentname.split('.fa')[0]+ '_' + str(count5) + '.tre'
                        count5+=1
                        s.unroot()
                        s.write(format=1, outfile=tonametree)
            if found == True:
                allindividuals = list(t.get_leaf_names())
                childinds = list(childtokeep.get_leaf_names())
                todrop=[y for y in inds if y not in childinds]
                tokeep=[x for x in allindividuals if x not in todrop]
                t.prune(tokeep, preserve_branch_length=True)


    # check to see if the tree contained only lineage-specific uplicates
    node2labels = t.get_cached_content(store_attr="name")
    leaves = t.get_leaves()
    spec = t.iter_species()
    speclist = []
    for specname in spec:
        speclist.append(specname)
    species = speclist
    if len(leaves) == len(species) and len(species)>=options.minimum:
        treefiledest = options.outputtrees+'/'+alignmentname.split('.fa')[0]+ '_' + str(count5) +'.tre'
        t.unroot()
        t.write(outfile=treefiledest)
