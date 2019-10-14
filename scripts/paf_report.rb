#!/usr/bin/env ruby

require 'optparse'
QUERY_NAME = 0
QUERY_LENGTH = 1
REGION_SIZE = 0
TAG = 1
VALIDATION = 2
MIN = 0
MAX = 1


##############################################################################################################################
## METHODS
##############################################################################################################################

def load_file(file) #files is an array of files
	queries = {}
	cigar_field = nil
	File.open(file).each do |line|
		line.chomp!
		fields = line.split("\t")
		fields[QUERY_LENGTH] = fields[QUERY_LENGTH].to_i
		queries[fields[QUERY_NAME]] = [nil, nil, nil, nil, fields]
	end
	return queries
end

def find_cigar_index(hit) #find the default field of cigar in paf; it returns field index
	cigar_index = nil
	hit.each_with_index do |field, i|
		if field.to_s.include?('cg:Z:')
			cigar_index = i
			break
		end
	end
	p cigar_index if cigar_index.nil?
	return cigar_index
end

def curate_all_cigar(queries, min_match) #format all cigar strings into bidimensional array and curate usin minimun match size threshold
	cigar_field = find_cigar_index(queries.first[1][4])
	queries.each do |query_name, attributes|
		identity, coverage, exons, cigar, hit = attributes
		raw_cigar = hit[cigar_field]
		formatted_cigar = format_cigar(raw_cigar)
		curated_cigar = curate_cigar(formatted_cigar, min_match)
		curated_cigar = curate_cigar(curated_cigar.reverse, min_match).reverse
		queries[query_name][3] = curated_cigar
	end
	return queries
end


def format_cigar(cigar) #format a cigar string and it returns a bidimensional array
	formatted_cigar = []
	temp_size = []
	cigar = cigar.split("cg:Z:").last.split('')
	cigar.each do |letter|
		if letter !~ /[A-Z]/
			temp_size << letter
		elsif letter != "S"
			formatted_cigar << [temp_size.join("").to_i, letter]
			temp_size = []
		elsif letter == "S"
			temp_size = []
		end
	end
	return formatted_cigar
end

def curate_cigar(cigar, threshold)
	regions_to_delete = []
	cigar.each_with_index do |region, i|
		region_size, tag = region
		if (tag == "M" && region_size >= threshold)
			break
		else
			cigar[i] = nil
		end
	end
	cigar.reject! do |position|
		position.nil?
	end

	return cigar
end

def deformat_cigar(cigar_array)
	cigar_string = "cg:Z:"
	cigar_array.each do |region_size, tag|
		cigar_string << "#{region_size.to_s}#{tag}"
	end
	return cigar_string
end

def get_alignment_size(cigar)
	alignment_size = 0
	cigar.each do |region_size, tag|
		next if tag == "N"
		alignment_size += region_size
	end
	return alignment_size
end

def get_matches(cigar)
	matches = 0
	cigar.each do |region_size, tag|
		next if tag != "M"
		matches += region_size
	end
	return matches
end

def get_exon_number(cigar)
	count = 1
	cigar.each do |region_size, tag|
		count += 1 if tag == "N"
	end
	return count
end

def calculate_stats(queries) 
	queries.each do |query_name, attributes|
		identity, coverage, exons, cigar, hit = attributes
		nt_matches = get_matches(cigar)
		query_length = hit[QUERY_LENGTH]
		aligment_length = get_alignment_size(cigar)
		attributes[0] = nt_matches.fdiv(aligment_length)
		attributes[1] = nt_matches.fdiv(query_length)
		attributes[2] = get_exon_number(cigar)
	end
end

def print_paf(queries, filename)
	filename = "#{filename}.paf" if !filename.include?('.paf')
	File.open("#{filename}", 'w') do |file|
		queries.each do |query_name, attributes|
			hit = attributes.last
			file << "#{hit.join("\t")}\n"
		end
	end
end

def filter_queries(queries, options)
	queries.reject! {|query, attributes| 
		identity, coverage, exons, cigar, hit = attributes
		check_query(identity, coverage, exons, options)
	} if options[:select_mode] == 'k' 


	queries.select! {|query, attributes| 
		identity, coverage, exons, cigar, hit = attributes
		check_query(identity, coverage, exons, options)
	} if options[:select_mode] == 'r'
	return queries
end

def check_query(identity, coverage, exons, options)
	return ((( coverage < options[:coverage][MIN] || coverage > options[:coverage][MAX] ) ||   
		( identity < options[:identity][MIN] ||  identity > options[:identity][MAX] ) ) ||
		( exons < options[:exons][MIN] || exons > options[:exons][MAX] ))
end
	
def print_results(queries, print_header, tag)
	if print_header
		if ! tag.nil?
			puts %w(query Identity Coverage Exons tag).join("\t")
		else 	
			puts %w(query Identity Coverage Exons).join("\t")
		end
	end
	queries.each do |query, attributes|
		identity, coverage, exons, cigar, hit = attributes
		if !tag.nil?
			puts "#{query}\t#{[identity, coverage, exons, tag].join("\t")}" if !cigar.empty?
		else
			puts "#{query}\t#{[identity, coverage, exons].join("\t")}" if !cigar.empty?
		end
	end
end

#############################################################################################################################
## INPUT PARSING
##################################################################################################################################3
options = {}

OptionParser.new do  |opts|
	options[:input] = nil 
	opts.on("-i FILE", "--file", "Select the file to extract info.") do |i|
		options[:input] = i	
	end

	options[:min_match] = 50
	opts.on("-m INT", "--min_full_match INT", "Set the minimun size of local full match region. DEFAULT 50") do |i|
		options[:min_match] = i.to_i
	end

	options[:identity] = [0, 1]
	opts.on("-I 'MIN,MAX'", "--filter-identity 'MIN,MAX'", "Set identity conditions. DEFAULT '0,1'") do |i|
		options[:identity] = i.split(',').map {|coord| coord.to_f}
	end

	options[:coverage] = [0, 1]
	opts.on("-C 'MIN,MAX'", "--filter-coverage 'MIN,MAX'", "Set coverage conditions. DEFAULT '0,1'") do |c|
		options[:coverage] = c.split(',').map {|coord| coord.to_f}
	end

	options[:header] = false
	opts.on("-H", "--header", "Print header") do
		options[:header] = true
	end

	options[:select_mode] = "k"
	opts.on("-S STRING", "--selection_mode STRING", "Set the selection mode for remove or keep the alignments that meet the conditions ('r' for remove and 'k' for keep). DEFAULT 'k'") do |mode|
		options[:select_mode] = mode
	end

	options[:exons] = [0,100000]
	opts.on("-E 'MIN,MAX'", "--exons 'MIN,MAX'", "Set exon numbers conditions. CIGAR 'cg' is needed. DEFAULT '0,1'") do |e|
		options[:exons] = e.split(',').map {|coord| coord.to_i}
	end

	options[:paf] = nil
	opts.on("-p STRING", "--paf  STRING", "Return filtered PAF") do |p| 
		options[:paf] = p
	end
	
	options[:tag] = nil
	opts.on("-t string", "--tag string", "Set tag for all entries") do |tag|
		options[:tag] = tag
	end

	# options[:overlap] = nil
	# opts.on("-o INTEGER", "--overlap INTEGER", "Report groups of overlaping hits in INTEGER positions. '0' for any overlapping") do |i| 
	# 	options[:overlap] = i.to_i
	# end
	
	opts.on("-h", "--help", "Displays helps") do 
		puts opts
		abort()
	end
end.parse!

#############################################################################################################################
## MAIN PROGRAM
##################################################################################################################################3

abort("ERROR: You must load a more recent version o f Ruby [ > 1.9.3 ]") if RUBY_VERSION.to_f <= 2
queries = load_file(options[:input])
queries = curate_all_cigar(queries, options[:min_match])
#p queries

# queries.each do |query_name, attributes| 
# 	identity, coverage, exons, cigar, hit = attributes
# 	p "#{hit[0,17].join("\t")}\t#{deformat_cigar(cigar)}"
# end
# queries = format_queries(queries_array)



calculate_stats(queries)


#filtered_queries = filter_queries(queries, options)

# # if !options[:overlap].nil?
# # 	scaffolds = find_overlaping_queries(filtered_queries, options[:overlap])
# # 	report_overlaping(scaffolds, "#{options[:input].first.split(".").first}_overlapping")
# # end


#print_paf(filtered_queries, options[:paf]) if !options[:paf].nil?

print_results(queries, options[:header], options[:tag])




























# def did_overlap?(control_hsp, target_hsp, start_position, end_position, return_positions = false)
# 	overlap = false
# 	overlapping_positions = nil
# 	if control_hsp[start_position] <= target_hsp[start_position] && 
# 		control_hsp[end_position] > target_hsp[start_position] 
		
# 		if control_hsp[end_position] < target_hsp[end_position]
# 			overlapping_positions = [target_hsp[start_position], control_hsp[end_position]]
# 		else 
# 			overlapping_positions = [target_hsp[start_position], target_hsp[end_position]]
# 		end

# 	elsif control_hsp[start_position] >= target_hsp[start_position] &&
# 		control_hsp[start_position] < target_hsp[end_position] 

# 		if control_hsp[end_position] > target_hsp[end_position]
# 			overlapping_positions = [control_hsp[start_position], target_hsp[end_position]] 
# 		else
# 			overlapping_positions = [control_hsp[start_position], control_hsp[end_position]]
# 		end
# 	end
	
# 	if return_positions
# 		overlap = overlapping_positions
# 	else
# 		overlap = true if !overlapping_positions.nil?
# 	end

# 	return overlap
# end


# def find_overlaping_queries(queries, min_overlapping_length)
# 	#this method is for group the queries by alignment overlapping 
# 	scaffolds = {}
# 	query_name = 0
# 	start_position = 1
# 	end_position = 2
# 	all_targets = group_by_target(queries)
# 	all_targets.each do |target_name, formatted_queries|
# 		formatted_loci = {}
# 		raw_loci = find_overlapping(formatted_queries, min_overlapping_length)
		
# 		locus_number = 1
# 		raw_loci.each do |locus|
# 			locus_name = "#{target_name}:#{locus_number}"
# 			formatted_loci[locus_name] = locus
# 			locus_number += 1
# 		end

# 		scaffolds[target_name] = formatted_loci 	
# 	end
# 	return scaffolds
# end
	
# def find_overlapping(queries, min_overlapping_length)
# 	query_name = 0
# 	start_position = 1
# 	end_position = 2
# 	queries.sort_by! { |query_name, start_position, end_position, attributes| 
# 		[start_position, end_position]
# 	}
# 	last_query = queries.shift
# 	loci = [[last_query]]
# 	queries_left = queries.length

# 	queries_left.times do
# 		last_query = queries.shift
# 		overlapping_positions = nil
# 		loci.each do |locus|
# 			locus.each do |saved_query|
# 				overlapping_coords = did_overlap?(last_query, saved_query, start_position, end_position, true)
# 				next if overlapping_coords.nil?
# 				overlapping_positions = overlapping_coords[1] - overlapping_coords[0]
# 				if overlapping_positions < min_overlapping_length
# 					overlapping_positions = nil
# 					next
# 				else 
# 					break
# 				end
# 			end
# 			if !overlapping_positions.nil?
# 				locus << last_query
# 				break
# 			end
# 		end
# 		loci << [last_query] if overlapping_positions.nil?
# 	end
# 	return loci
# end

# def group_by_target(queries)
# 	targets = {}
# 	queries.each do |query_name, hits|
# 		hits.each do |target_name, attributes|
# 			ident, cover, hsps = attributes
# 			sorted_hsps = hsps.sort_by { |hsp| hsp[TARGET_START] }
# 			hit_start = sorted_hsps.first[TARGET_START]
# 			hit_end = sorted_hsps.last[TARGET_END]
# 			hit_name = target_name.split(":").first
# 			if !targets[hit_name].nil?
# 				targets[hit_name] << [query_name, hit_start, hit_end, attributes]
# 			else
# 				targets[hit_name] = [[query_name, hit_start, hit_end, attributes]]
# 			end
# 		end	

# 	end
# 	return targets
# end

# def report_overlaping(scaffolds, filename)
# 	File.open(filename, 'w') do |file|
# 		scaffolds.each do |scaffold_name, loci|
# 			loci.each do |locus_name, hits|
# 				hits.sort_by! { |hit_name, target_start, target_end, attributes| target_start }
# 				locus_start = hits.first[1]
# 				locus_end = hits.last[2]
# 				hits_names = []
# 				hits.each do |hit| 
# 					hits_names << hit[0]
# 				end
# 				hits_count = hits_names.length
# 				hits_names = hits_names.join(";")
# 				file << "#{locus_name}\t#{locus_start}\t#{locus_end}\t#{hits_count}\t#{hits_names}\n"
# 			end
# 		end
# 	end
# end
