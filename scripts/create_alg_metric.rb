#!/usr/bin/env ruby

require 'optparse'


V4_SIZE = 697125

##################################################################
### METHODS
##################################################################


def load_sample_list(file)
	sample_list = File.readlines(file).map{ |line| line.chomp!}
	return sample_list
end

def count_transcripts(file)
	all_transcripts = 0
	transcripts_with_reads = 0
	count_table = File.readlines(file).map do |line| 
		all_transcripts += 1
		line = line.chomp.split("\t")
		transcripts_with_reads += 1 if line[1].to_i > 0 
	end
	return all_transcripts, transcripts_with_reads
end


def load_reads_stats(file)
	all_reads = nil
	not_mapped_reads = nil
	mapped_reads = nil
	File.readlines(file).map do |line| 
		line.chomp!.gsub!(/^ */,'')
		all_reads = line.split(" ").first.to_i if line.include?("reads; of these:")
		not_mapped_reads = line.split(" ").first.to_i if line.include?("aligned concordantly 0 times")
		break if !all_reads.nil? && !not_mapped_reads.nil?
	end
	abort("#{file} not found") if all_reads.nil? || not_mapped_reads.nil?
	mapped_reads = all_reads - not_mapped_reads
	return all_reads, mapped_reads
end

def get_information(transcriptomes, samples, actual_path)
	
	sample_info  ={}
	samples.each do |sample|
		
		transcriptomes.each do |transcriptome|
			sample_info[transcriptome] = [["Mapped reads"],["Reads per TTs"],["Mapped TTs"],["Mapped TTs [v4.0]"]] if !sample_info[transcriptome]
			path = "#{actual_path}/#{transcriptome}/#{sample}"
			bowtie_log = "#{path}/bowtie2_0000/bowtie_log"
			count_matrix = "#{path}/sam2counts_all.py_0000/matrix.txt"
			abort("#{count_matrix} not founded") if !(File.file?(count_matrix) && !File.zero?(count_matrix))
			all_reads, mapped_reads = load_reads_stats(bowtie_log)
			all_transcripts, transcripts_with_reads = count_transcripts(count_matrix)
			mapped_reads_ratio = mapped_reads.fdiv(all_reads)
			mapping_ratio = mapped_reads.fdiv(transcripts_with_reads)
			transcripts_with_reads_ratio = transcripts_with_reads.fdiv(all_transcripts)
			validated_contigs_ratio = transcripts_with_reads.fdiv(V4_SIZE)
			#puts [sample,transcriptome,mapped_reads_ratio,mapping_ratio,transcripts_with_reads,validated_contigs_ratio].join("\t")
			sample_info[transcriptome][0] <<mapped_reads_ratio
			sample_info[transcriptome][1] <<mapping_ratio
			sample_info[transcriptome][2] <<transcripts_with_reads_ratio
			sample_info[transcriptome][3] <<validated_contigs_ratio
		end
	end
	return sample_info
end

def calculate_stats(all_samples)
	all_samples.each do |transcriptome, data|
		data.each do |values|
			metric = values.shift
			# p values
			mean, se = mean_se(values)
			puts [transcriptome, metric,mean,se].join("\t")
		end
	end
end




def mean_se(array)
	total = array.inject(0, :+)
	size = array.length
    mean  = total.fdiv(size)
    variance = array.inject(0) { |variance, value| variance + (value - mean) ** 2 }.fdiv(size - 1)

    se = Math.sqrt(variance).fdiv(Math.sqrt(size))
	return mean, se
end

##################################################################
### OPTS
##################################################################

options = {}

OptionParser.new do  |opts|
	options[:samples] = ''
	opts.on("-s FILE", "--sample_list", "File including all samples listed.") do |i|
		options[:samples] = i
	end

	options[:path] = Dir.pwd
	opts.on("-p PATH", "--path", 'Set the execution path') do |path|
		path.gsub!(/\/\z/, "")
		options[:path] = path
	end
	
	options[:header] = false
	opts.on("-H", "--header", "Print header") do
		options[:header] = true
	end
	
	options[:transcriptomes] = ''
	opts.on("-t TRANS1,TRANS2,TRANS3", "--transcriptomes TRANS1,TRANS2,TRANS3", "Transcriptomes to process 'TRANS1,TRANS2,TRANS3'") do |i|
		options[:transcriptomes] = i.split(",")
	end

	opts.on("-h", "--help", "Displays helps") do 
		puts opts
		abort
	end
end.parse!


##################################################################
### MAIN
##################################################################

samples = load_sample_list(options[:samples])
puts ["transcriptome","metric", "mean", "se"].join("\t") if options[:header]
all_samples = get_information(options[:transcriptomes], samples, options[:path])
calculate_stats(all_samples)