#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem='4gb'
#SBATCH --constraint=cal
#SBATCH --time='7-00:00:00'


transcriptome_name='' ### SET TRANSCRIPTOME NAME WITHOUT FILE EXTENSION
genome='./genome.fasta.fasta' ###SET GENOME PATH
files_illumina=`pwd`'/paired_1_.fastq,'`pwd`'/paired_2_.fastq' ### SET ILLUMINA FILES PATH
files_454=`pwd`'/sequences_.fastq' ###SET ROCHE/454 FILES PATH

####### CORE SETTINGS FOR FULL-LENGTHER NEXT
srun hostname -s > workers_$transcriptome_name

####### DEPENDENCIES
module load minimap2/2.9
module load samtools/1.8
. ~soft_bio_267/initializes/init_R
. ~soft_bio_267/initializes/init_ruby
 . ~soft_cvi_114/initializes/init_fln
export PATH=`pwd`/scripts:$PATH

###### INITIALIZING WORKFLOW
transcriptome='./'$transcriptome_name'.fasta' #transcriptome path
mkdir $transcriptome_name
cd $transcriptome_name
ln -s $genome ./genome.fasta
ln -s $transcriptome ./transcriptome.fasta
report_path='../report'
if [ ! -d $report_path ]; then
	mkdir $report_path
fi 

###### ANNOTATION
full_lengther_next -s 10.243 -f transcriptome.fasta -z -u `pwd`/Actinopterygii_db/Actinopterygii_db -g vertebrates -c 1500 -q 'd' -w ../workers_$1 -M "$files_illumina;$files_454"

###### MAPPING the transcriptome against the genome (-N0 + --secondary=no = only primary; -u f force minimap2 to consider the forward transcript strand only)
 minimap2 --cs -ax splice -uf --secondary=no -C5 genome.fasta transcriptome.fasta > raw_alignment.sam

samtools view -h -F 2052 raw_alignment.sam > alignment.sam
samtools view -f 4 raw_alignment.sam > unmapped.sam

paftools.js sam2paf alignment.sam > alignment.paf
merge_paf_cigar.rb alignment.paf alignment.sam alignment_cigar.paf

if [ ! -s alignment_cigar.paf ] ; then
	echo "ERROR: minimap2 has failed"
	exit 1
else 

###### METRIC FILES GENERATION
	rm alignment.paf alignment.sam
	mkdir reports
	paf_report.rb -i alignment_cigar.paf > alignment.txt
	echo -e "length\ttag" > length.list
	echo -e "query\tIdentity\tCoverage\tExons\tCategory" > tagged_transcripts.txt

	# #### Raw transcriptome
	echo "Raw transcriptome"
	grep ">" $transcriptome | tr -d ">" > Raw.list
	fln_protein_class.rb -O -i Raw.list  -f 'c' -p "`pwd`/fln_results" -H -n "Raw" 1>Protein_stats 2>>ort_rep

	# #### Mapped transcripts
	echo "Mapped transcriptome"
 	cut -f 1 alignment.txt > mapped.list

	# #### Unmapped transcripts
	echo "Unmapped transcriptome"

	cut -f 1 unmapped.sam > Unmapped.list
	fasta_editor.rb -i transcriptome.fasta -l Unmapped.list -L unmapped.length
		
	# #### Typical transcripts
	echo "Typical transcriptome"
	cat alignment.txt | awk 'BEGIN {OFS="\t"}{if (!($2 > 0.9 && $4 == 1)) print $1,$2,$3,$4,"Typical TTs"}' > typical.txt
	cut -f 1 typical.txt > Typical.list
	fln_protein_class.rb -O -i Typical.list  -f 'c' -p "`pwd`/fln_results" -n "Typical" 1>> Protein_stats 2>>ort_rep

	# #### Atypical transcripts
	echo "Atypical transcriptome"
	cat alignment.txt | awk 'BEGIN {OFS="\t"}{if ($2 > 0.9 && $4 == 1) print $1,$2,$3,$4,"Atypical TTs"}' >> tagged_transcripts.txt
	grep "Atypical" tagged_transcripts.txt | cut -f 1 > Atypical.list
		
	# #### Low quality transcripts
	echo "Low quality transcriptome"
	cat typical.txt | awk 'BEGIN {OFS="\t"}{if (!($2 > 0.7 && $3 > 0.7)) print $1,$2,$3,$4,"Low Quality TTs"}' >> tagged_transcripts.txt
	grep "Low" tagged_transcripts.txt | cut -f 1 > Low.list
		
	# #### Reliable transcripts
	echo "Reliable transcriptome"
	cat typical.txt | awk 'BEGIN {OFS="\t"}{if (($2 > 0.7 && $3 > 0.7)) print $1,$2,$3,$4,"Reliable TTs"}' >> tagged_transcripts.txt
	grep "Reliable" tagged_transcripts.txt | cut -f 1 > pre_reliable.list
	cut -f 1 fln_results/misassembled.txt > artifacts.list
	grep -v "Query_id" fln_results/unmapped.txt | cut -f 1 >> artifacts.list
	fasta_editor.rb -i $transcriptome -l pre_reliable.list -d artifacts.list -c a -o reliable_transcriptome.fasta
	grep ">" reliable_transcriptome.fasta | tr -d ">" > Reliable.list
	fln_protein_class.rb -O -i Reliable.list  -f 'c' -p "`pwd`/fln_results" -n "Reliable" 1>> Protein_stats 2>>ort_rep
		
###### PLOTS GENERATION
	echo "Generating Plots"
	prot_ylim=155000
	if [ $transcriptome_name == "v4_transcriptome" ]; then
		prot_ylim=`grep -c ">" $transcriptome`
		transcriptome_name="v4.0"
	fi
	cat Protein_stats | sed 's/C_terminal/C-terminal/' | sed 's/N_terminal/N-terminal/' |sed 's/nc_rnas/ncRNAs/' | sed 's/new_coding/New coding/' | sed 's/unknown/Unknown/' | sed 's/category/Category/' | sed 's/sequences/Sequences/' | sed 's/Raw/Raw TTs/'| sed 's/Typical/Typical TTs/'| sed 's/Reliable/Reliable TTs/'| sed 's/v4_transcriptome/v4.0/'> renamed_protein_stats
	point_density.R -d renamed_protein_stats -Y $prot_ylim -n $transcriptome_name -o $report_path/prot_density_$transcriptome_name -x sample -y Sequences -c Category -l
	point_density.R -d renamed_protein_stats -n $transcriptome_name -o $report_path/prot_legend -x sample -y Sequences -c Category 
	xyplot_graph.R -d tagged_transcripts.txt -x Identity -y Coverage -g Category -o $report_path/scatter_plot_$transcriptome_name
	if [ $transcriptome_name != "v4_transcriptome" ]; then
		echo "test"
		 cat unmapped.length | awk 'BEGIN {OFS="\t"}{print $2,"Unmapped TTs"}' >> length.list
		 extract_length_dstr.rb alignment_cigar.paf Atypical.list "Atypical TTs" >> length.list
		 extract_length_dstr.rb alignment_cigar.paf Low.list "Low Quality TTs" >> length.list
		 extract_length_dstr.rb alignment_cigar.paf Reliable.list "Reliable TTs" >> length.list
		 boxplo1D.R -d length.list -o $report_path/length_dist_$transcriptome_name -x tag -y length -t $transcriptome_name >& Length.rep
	fi

###### REDUNDANCY CALCULATION
cat *_orth | awk -vtranscriptome="$transcriptome_name" 'BEGIN {OFS="\t"}{print $1,$2,transcriptome}' >> ../redundancy.txt
###### AFTER THIS LAST COMMAND, make_red_graph.sh MUST BE LAUNCHED






