"""This script takes as input a folder with gene trees and a folder with alignments.
The user must supply a minimum number of species for which sequences are required.
Then, the script save trees that contained at least the specified number of species.
Finally, the user must supply an output file for filtered trees. Sequences must be named
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


def get_species_name(node_name_string):
    """this function will get species names for the leaves of a tree."""
    # Species code is the first part of leaf name (separated by an
    #  underscore character)
    spcode = node_name_string.split("_")[0:2]
    spcode = spcode[0]+'_'+spcode[1]
    # We could even translate the code to complete names
    return spcode


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

    # get the leaf node to leaf labels
    node2labels = t.get_cached_content(store_attr="name")

    # get the node to species labels
    node2species = t.get_cached_content(store_attr="species")

    
    # check to see if the tree contains only single-copy orthologs and filter for missing data
    node2labels = t.get_cached_content(store_attr="name")
    leaves = t.get_leaves()
    species = t.get_species()
    if len(species)>=options.minimum:
        alignmentname = alignment.split('/')[-1]
        treefiledest = options.outputtrees+'/'+alignmentname.split('.fa')[0]+'.tre'
        t.write(outfile=treefiledest)
        os.system('cp %s/%s %s/%s' % (options.alignmentfolder, alignmentname, options.outputtrees, alignmentname))


myoutputtrees = os.listdir(options.outputtrees)
with open('%s.tre.tmp' % options.outputtrees, 'w') as outfile:
    for fname in myoutputtrees:
        if fname.endswith('.tre'):
            with open(os.path.join(options.outputtrees,fname)) as infile:
                outfile.write(infile.read())
                outfile.write('\n')


os.system('mv %s.tre.tmp %s.tre' % (options.outputtrees, options.outputtrees))

