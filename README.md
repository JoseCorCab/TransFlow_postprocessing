# TransFlow_postprocessing
Annotation and post-processing of de novo transcriptome assemblies.

This workflow is launched after transcriptome assembly. In this order it performs annottation using Full-Lengther Next, 
mapping onto genome, obtaining of reliable transcripts and generation of all graphs showed in the article.

USAGE: 
1.- Link all data: de novo assembled transcriptomes, genome, and reads used for assembly.
2.- Launch annotation_wf.sh, and you need to add transcriptome name (without '.fasta' extension) as first argument.
3.- When execution has finished, you can launch make_red_graph.sh for generating redundancy graph.
