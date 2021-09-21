# Primate_Paralogs

Scripts for 'Using all gene families vastly expands data available for phylogenomic inference in primates'

For all scripts, sequences must be names Species_genus_otherinfo, as the info separated by the first underscore will be used to assign gene copies to species.

singlecopyorthologs_v1a.py: Script for filtering only single copy orthologs with some sampling threshold (SCOs).

allparalogs_v1a.py: Script for filtering gene trees to keep all paralogs with some minimum taxon sampling threshold (ALL PARALOGS)

oneparalogs_v1a.py: Script to sample one paralog per species (ONE PARALOGS)

lineagespecificdups_v3.py: Script to perform Lineage Specific Duplicate filtering (LSDs)

twolineagedups_v3.py: Script to perform Two-lineage and lineage-specific duplicate filtering (TSDs)

SEbranchcutting_v2a.py: Script to perform Subtree Extraction filtering (SE)

subset_MI.py: Scripts to subsample one ortholog per orthogroup for MI filtering.

collapse_sn0_nodes.py: Script to collapse gene tree nodes with no decisive sites.

resampling.py: Script to test for introgression using GCFs.

cliptsds_v1a.py: Script to prune lsds and tsds from gene trees which can then be used in further filtering. DOES NOT FILTER FOR ONLY THOSE TREES WITH TSDS, LSDS, or SCOs (see twolineagedups_v3.py for this funciton.
