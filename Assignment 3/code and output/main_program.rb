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

no_exons = []
out_file = File.new("4a_genes_cttctt.gff", "w")
out_file2 = File.new("4b_genes_no_cttctt.txt", "w")


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
    no_exons.push name
    p "The gen #{name} do NOT have exons with the target_sequence"
  else
    p "We are working with the gen #{name}"
    p "The positions of the target sequence the exons of the strand + of the gen #{name} are: #{position_plus}"
    p "The positions of the target sequence the exons of the strand - of the gen #{name} are: #{position_minus}"

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
        puts assoc
        head = position.split('..')[0].to_i
        tail = position.split('..')[1].to_i
        #  - strand
        if /-/ =~ assoc
          puts "#{name}\t.\trepeat_region\t#{head}\t#{tail}\t.\t-\t.\t.\n\n"
          out_file.puts"#{name}\t.\trepeat_region\t#{head}\t#{tail}\t.\t-\t.\t.\n\n"
        else
          puts "#{name}\t.\trepeat_region\t#{head}\t#{tail}\t.\t+\t.\t.\n\n"
          out_file.puts "#{name}\t.\trepeat_region\t#{head}\t#{tail}\t.\t+\t.\t.\n\n"
        end
      end
    end
  end
end

out_file2.puts "This genes #{no_exons} does not have exons with cttctt motif"

out_file.close
out_file2.close
