#!/usr/bin/env ruby

#return a file which contains the clasification of proteins from the input list according to a previous full_lengther_next's annotation

require 'optparse'
#############################################################################################################################
## METHODS
#############################################################################################################################
def load_sequences_list(filename)
	sequences_list = []
	File.open(filename).each do |protein|
		sequences_list << protein.chomp
	end
	
	return sequences_list
end

def load_all_sequences_classified(fln_path)
	all_seqs = {}
	line_number = 1
	File.open("#{fln_path}pt_seqs").each do |line|
		if line_number > 1
			line.chomp!
			fields = line.split("\t")
			all_seqs[fields[0]] = "#{fields[4].split(" ").first}"
			
		end
		line_number += 1
	end

	%w(unmapped misassembled nc_rnas new_coding unknown).each do |file|
		line_number = 0 
		File.open(File.join(fln_path, file + '.txt')).each do |line|
			line_number += 1
			next if file != 'unknown' && file != 'misassembled'  && line_number == 1
			line.chomp!
			seq_name = line.split("\t").first
			if ["misassembled", "unmapped"].include?(file)
				all_seqs[seq_name] = "Artifacts"

			else
				all_seqs[seq_name] = file
			end
		end
	end
	return all_seqs
end

def count_orthologs(fln_path, sequences_list, tag)
	red_outfile = "#{tag}_orth"
	all_orthologs = []
	data_orthologs = load_orthologs(fln_path)
	sequences_list.each do |sequence|
		orthologs_name = data_orthologs[sequence]
		all_orthologs << orthologs_name if !orthologs_name.nil?
	end
	File.open(red_outfile, 'w') do |outfile|
		STDERR.puts "Total orthologs of #{tag}\t#{all_orthologs.length}"
		STDERR.puts "Different orthologs of #{tag}\t#{all_orthologs.uniq.length}"
		redundancy = all_orthologs.length - all_orthologs.uniq.length
		outfile.puts "#{redundancy}\t#{tag}"
	end
end

def load_orthologs(fln_path)
	orthologs = {}
	all_orthologs = []
	line_number = 1
	File.open("#{fln_path}pt_seqs").each do |line|
		if line_number > 1
			line.chomp!
			fields = line.split("\t")
			orthologs[fields[0]] = fields[2] 
			all_orthologs << fields[2]
		end
		line_number += 1
	end
	return orthologs
end

def choose_and_save_proteins(sequences_list, all_sequences_classified)
	sequences_classified = []
	sequences_list.each do |sequence|
		query = all_sequences_classified[sequence]
		query = "unknown" if query.nil?
		sequences_classified << [sequence, query]
		
	end
	return sequences_classified
end

def get_stats(sequences_list_classified, format)
	stats = []
	all_sequences_count = sequences_list_classified.count
	counts = Hash.new{0}
	sequences_list_classified.each do |seq, type|
		counts[type] += 1
	end
	counts.each do |category, count|

		stats << [category, count.fdiv(all_sequences_count).round(4)] if format == "p"
		stats << [category, count] if format == "c"
	end
	return stats
end

def save_stats(stats, sample_name)
	stats.sort.each do |category, percentage|
		puts "#{sample_name}\t#{category.gsub(/\s/, "_")}\t#{percentage}"
	end
end

#############################################################################################################################
## INPUT PARSING
#############################################################################################################################

options = {}

OptionParser.new do  |opts|
	options[:input] = ''
	opts.on("-i FILE", "--input_file", "Select the list of proteins to clasificate.") do |i|
		options[:input] = i
	end

	options[:path] = ''
	opts.on("-p PATH", "--fln-path", 'Set the path of the full_lengther_next results EXAMPLE: -p "~/full_lengther_execution/fln_results"') do |path|
		path << "/" if path !~ /\/\z/
		options[:path] = path
	end
	
	options[:header] = false
	opts.on("-H", "--header", "Print header") do
		options[:header] = true
	end
	options[:format] = "p"
	opts.on("-f STRING", "--format STRING", "Print output in percentage 'p' or proteins count 'c'") do |s|
		options[:format] = s
	end
	options[:name] = ''
	opts.on("-n STRING", "--sample-name STRING", "Set the sample name. DEFAULT = filename") do |name|
		options[:name] = name
	end

	options[:find_orthologs] = false
	opts.on("-O", "--orthologs", "Count different ortologous") do
		options[:find_orthologs] = true
	end

	opts.on("-h", "--help", "Displays helps") do 
		puts opts
		abort
	end
end.parse!

#############################################################################################################################
## MAIN PROGRAM
#############################################################################################################################
abort("ERROR: The specified files not exist") if options[:input].empty?
sequences_list = load_sequences_list(options[:input])
fln_data = load_all_sequences_classified(options[:path])
count_orthologs(options[:path], sequences_list,options[:name]) if options[:find_orthologs]
sequences_list_classified = choose_and_save_proteins(sequences_list, fln_data)
# sequences_list_classified.each do |seq, category|
# 	puts "#{seq}\t#{category.split(";").join("\t")}\n"
# end
stats = get_stats(sequences_list_classified, options[:format])
options[:name] = options[:input].split(".").first if options[:name].empty?
if options[:header]
	puts "sample\tcategory\tpercentage" if options[:format] == "p"
	puts "sample\tcategory\tsequences" if options[:format] == "c"
end
#
save_stats(stats, options[:name])
