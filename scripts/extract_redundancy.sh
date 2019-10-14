#!/usr/bin/env bash

for file in Typical.list Raw.list Reliable.list; do 
	type=`echo $file |cut -f 1 -d "."`
	prot_list=`extract_pt_seqs.rb ./fln_results/pt_seqs $file| cut -f 3`
	all_prots=`echo $prot_list| wc -l`
	uniq_prots=`echo "prot_list" |sort | uniq |wc -l`
	redundance=$((all_prots - uniq_prots))
	echo -e $redundance"\t"$type"\t"$1
done
