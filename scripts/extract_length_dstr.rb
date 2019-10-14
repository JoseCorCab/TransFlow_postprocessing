#!/usr/bin/env ruby

full_paf = ARGV[0] 
list = ARGV[1]
tag = ARGV[2]	

tr_list = File.readlines(list).map { |line| line.chomp}

File.open(full_paf).each do |line|
	line = line.chomp.split("\t")
	puts [line[1], tag].join("\t") if tr_list.include?(line[0])
end
