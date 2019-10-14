#!/usr/bin/env ruby


require 'bio'

def load_list(filename)
	list = File.readlines(filename)
	return list
end

def analyse_input_file(input_file, list, outfile, print_length = false)
  biofastafile = Bio::FlatFile.open(Bio::FastaFormat, input_file)
  File.open(outfile, 'w') do |out_file|
	  biofastafile.each_entry do |entry|
	  	id = entry.entry_id
	  	seq = entry.seq
		if list.include?(entry.entry_id)
			outfile.puts ">#{id}", seq
			puts "#{id}\t#{seq.length}" if print_length
		end
	  end
	end
end

# def save_file(filtered_fasta, outfile, print_length = false)
# 	File.open(outfile, 'w') do |file|
# 		filtered_fasta.each do |id, seq|
# 			file.puts ">#{id}", seq
# 			puts "#{id}\t#{seq.length}" if print_length
# 		end
# 	end
# end



input_fasta = ARGV[0]
list_file = ARGV[1]
outfile = ARGV[2]
length_opt = false
length_opt = true if ARGV[3] == "-l"

list = load_list(list_file)
analyse_input_file(input_fasta, list, outfile, print_length = length_opt)
