# use httparty to reduce timeout error
require 'httparty'
require 'bio'


gene_names = []

#open the document to read name list
# File.open ('sample_list.txt') do |f|
File.open ('ArabidopsisSubNetwork_GeneList.txt') do |f|
  f.each_line do |line|
    name = line.delete "\n"
    gene_names.push name
  end
end

target_sequence ="CTTCTT"

# for task6, we can pass gene name AT2g46340 to get AT2g46340.gff
if ARGV[0]
  out_file = File.new("#{ARGV[0]}.gff", "w")
  gene_names = [ARGV[0]]
else
  out_file = File.new("5_cttctt_in_chromosome.gff", "w")
end


gene_names.each do |name|
  
  #catch call api timeout error
  begin
    uri = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{name}"
    res = HTTParty.get uri, timeout: 10
  rescue StandardError => e
    p "gene #{name} get api fail, exception is #{e.message}"
    next
  end  
  obj = Bio::EMBL.new res.body
  bio_sequence = obj.to_biosequence

  position_plus, position_minus = [], []

  
  chromosome = obj.ac[0].split ':'
  p obj.ac[0]

  ch_name, ch_number, ch_head, ch_tail = chromosome[1], chromosome[2], chromosome[3], chromosome[4],
  obj.features.each do |feature|

    assoc, position = feature.assoc.to_s, feature.position.to_s

    if /exon_id/ =~ assoc && /[A-Z]/ !~

      # get head and tail from string "1..2" => 1, 2
      position_tmp = position.delete('complement()');
      position_head = position_tmp.split('..')[0].to_i
      position_tail = position_tmp.split('..')[1].to_i

      # to get sequence of the exon, through 'complement' to identify strand
      if /complement/ =~ position
        strand = 'minus'
        bio_sequence_rev = bio_sequence.tr('acgt','tgca')
        bio_sequence_rev.reverse!
        exon = bio_sequence_rev[position_head..position_tail]
      else
        strand = 'plus'
        exon = bio_sequence[position_head..position_tail]
      end

      #to find the sites of the sequence cttctt in each exon
      position_exon = []
      positions_target_sequence_exon = (0...exon.length).find_all { |i| exon[i, target_sequence.length] == target_sequence.downcase }
      positions_target_sequence_exon.each do |position|
        tmp = []
        head = position + 1 + position_head
        tail = head - 1 + target_sequence.length
        tmp.push head, tail
        position_exon.push tmp
      end

      #collect to the corresponding arrays separately
      if strand == 'plus'
        position_plus.push *position_exon
        position_plus.uniq!
      else
        position_minus.push *position_exon
        position_minus.uniq!
      end

    end
  end

  if position_plus.empty? && position_minus.empty?
    p "El gen #{name} do NOT have exons with the target_sequence repeat"
  else
    p "We are working with the gen #{name}"
    p "Las posiciones de las secuencias objetivo en la cadena + del gen #{name} son: #{position_plus}"
    p "Las posiciones de las secuencias objetivo en la cadena - del gen #{name} son: #{position_minus}"

    # format array for bio api, [[0, 1]] => ['0..1']
    format_position_plus = position_plus.map { |i| "#{i[0]}..#{i[1]}"}
    format_position_minus = position_minus.map { |i| "#{i[0]}..#{i[1]}"}

    # add feature to bio_sequence object
    format_position_plus.each do |i|
      feature = Bio::Feature.new 'repeat', i
      feature.append Bio::Feature::Qualifier.new('repeat_motif', 'CTTCTT')
      feature.append Bio::Feature::Qualifier.new('strand', '+')
      bio_sequence.features.push feature
    end

    format_position_minus.each do |i|
      feature = Bio::Feature.new 'repeat', i
      feature.append Bio::Feature::Qualifier.new('repeat_motif', 'CTTCTT')
      feature.append Bio::Feature::Qualifier.new('strand', '-')
      bio_sequence.features.push feature
    end

    # output information to outfile
    obj.features.each do |feature|
      assoc, position = feature.assoc.to_s, feature.position.to_s
      if /repeat_motif/ =~ assoc
        head = position.split('..')[0].to_i
        tail = position.split('..')[1].to_i
        #  - strand
        if /-/ =~ assoc
          ch_minus_head = ch_tail.to_i - tail + 1
          ch_minus_tail = ch_minus_head + 6
          puts "#{ch_name}:#{ch_number}\t#{name}\trepeat_region\t #{ch_minus_head}\t #{ch_minus_tail} \t-\t.\t.\n\n"
          out_file.puts "#{ch_number}\t.\trepeat_region\t#{ch_minus_head}\t#{ch_minus_tail}\t.\t-\t.\t.\n"
        else
          ch_plus_head = ch_head.to_i + head - 1
          ch_plus_tail = ch_plus_head + 6
          puts "#{ch_name}:#{ch_number}\t#{name}\trepeat_region\t #{ch_plus_head}\t #{ch_plus_tail} \t+\t.\t.\n\n"
          out_file.puts "#{ch_number}\t.\trepeat_region\t#{ch_plus_head}\t#{ch_plus_tail}\t.\t+\t.\t.\n"
        end

      end
    end
  end
end

out_file.close

