require 'net/http'
require 'json'
require 'thread'

# Gen class
class Gen

  attr_accessor :gene_id, :protein_name, :kegg_id_pathway, :go_id_term

  def initialize(gene_id, protein_name, kegg_id_pathway, go_id_term)
    @gene_id, @protein_name, @kegg_id_pathway, @go_id_term = gene_id, protein_name, kegg_id_pathway, go_id_term
  end

end

# Network class
class Network

  attr_accessor :gene

  def initialize(gene)
    @gene = gene
  end

end

gene_names = []

#open the document to read name list
# File.open ('sample_list.txt') do |f|
File.open ('ArabidopsisSubNetwork_GeneList.txt') do |f|
  f.each_line do |line|
    name = line.delete "\n"
    gene_names.push name
  end
end

keggs, goes, protein_main_names = [], [], []

#use name to get necessary information through togows.org api  
gene_names.each do |name|

  uri1 = URI "http://togows.org/entry/ebi-uniprot/#{name}/accessions.json"
  uri2 = URI "http://togows.org/entry/ebi-uniprot/#{name}/dr.json"
  uri3 = URI "http://togows.org/entry/kegg-enzyme/ath:#{name}"

  res1, res2, res3 = '', '', ''
  threads = []
  threads << Thread.new { res1 = Net::HTTP.get_response uri1 }
  threads << Thread.new { res2 = Net::HTTP.get_response uri2 }
  threads << Thread.new { res3 = Net::HTTP.get_response uri3 }
  threads.each(&:join)


  # res1 = Net::HTTP.get_response uri1
  protein_names = JSON.parse res1.body

  # res2 = Net::HTTP.get_response uri2
  data = JSON.parse res2.body
  go = data[0]['GO']

  # res3 = Net::HTTP.get_response uri3
  lines = res3.body.split "\n"

  kegg = 'unknow kegg pathway'
  lines.each do |line|
    tmp = line.sub /PATHWAY     /,''
    kegg = tmp if tmp.size != line.size
  end

  protein_main_names.push protein_names[0][0]
  goes.push go
  keggs.push kegg

end

gene_intacts = []
gene_intact_nums = []
gene_no_intacts = []
protein_intacts = []
protein_no_intacts = []

thread_arr1 = []
res_arr1 = []

#through togows.org api to define which gene have interactions 
protein_main_names.each_with_index do |name, index|
  url = URI "http://togows.org/entry/ebi-uniprot/#{name}/dr.json"
  thread_arr1 << Thread.new { res_arr1[index] = Net::HTTP.get_response url }
  thread_arr1.each(&:join)
end

res_arr1.each_with_index do |res, index|
  data = JSON.parse res.body
  intAct = data[0]["IntAct"]

  if intAct
    gene_intacts.push gene_names[index]
    gene_intact_nums.push index
    intAct.each { |item| protein_intacts.push item[0] }
  else
    protein_no_intacts.push protein_main_names[index]
    gene_no_intacts.push gene_names[index]
  end
end  


# protein_main_names.each_with_index do |name, index|

#   url = URI "http://togows.org/entry/ebi-uniprot/#{name}/dr.json"
#   res = Net::HTTP.get_response url
#   data = JSON.parse res.body
#   intAct = data[0]["IntAct"]
#   if intAct
#     gene_intacts.push gene_names[index]
#     gene_intact_nums.push index
#     intAct.each { |item| protein_intacts.push item[0] }
#   else
#     protein_no_intacts.push name
#     gene_no_intacts.push gene_names[index]
#   end

# end


thread_arr2 = []
res_arr2 = []

# extract the second level of interactions
protein_intacts.each_with_index do |name, index|
  uri = URI "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{name}"
  thread_arr2 << Thread.new { res_arr2[index] = Net::HTTP.get_response uri }
  thread_arr2.each(&:join)
end  

intact_level2 = []
res_arr2.each_with_index do |res, index|
  data = res.body.split "\n"
  intact_first = []
  data.each do |item|
    field = item.split "\t"
    protein1 = field[0].sub 'uniprotkb:',""
    protein2 = field[1].sub 'uniprotkb:',""
    intact_first.push protein1, protein2
  end
  intact_level2.push intact_first.uniq
end  

# protein_intacts.each do |name|

#   uri = URI "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{name}"
#   res = Net::HTTP.get_response uri
#   data = res.body.split "\n"
#   intact_first = []
#   data.each do |item|
#     field = item.split "\t"
#     protein1 = field[0].sub 'uniprotkb:',""
#     protein2 = field[1].sub 'uniprotkb:',""
#     intact_first.push protein1, protein2
#   end
#   intact_level2.push intact_first.uniq

# end

intact_level3 = []

# extract the third level of interactions
intact_level2.each do |names|
  intact_first = []
  thread_arr3 = []
  res_arr3 = []
  names.each_with_index do |name, index|
    uri = URI "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{name}"
    thread_arr3 << Thread.new { res_arr3[index] = Net::HTTP.get_response uri }
    thread_arr3.each(&:join)
  end  
  res_arr3.each_with_index do |res, index|
    data = res.body.split "\n"
    data.each do |item|
      field = item.split "\t"
      protein1 = field[0].sub 'uniprotkb:',""
      protein2 = field[1].sub 'uniprotkb:',""
      intact_first.push protein1, protein2
    end
  end  
  # names.each do |name|
  #   uri = URI "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{name}"
  #   res = Net::HTTP.get_response uri
  #   data = res.body.split "\n"
  #   data.each do |item|
  #     field = item.split "\t"
  #     protein1 = field[0].sub 'uniprotkb:',""
  #     protein2 = field[1].sub 'uniprotkb:',""
  #     intact_first.push protein1, protein2
  #   end
  # end
  intact_level3.push intact_first.uniq
end

# create gen object using kegg and go information which got above
gen_objs = []
gene_names.each_with_index do |name, index|
  gen = Gen.new name, protein_main_names[index], keggs[index], goes[index]
  gen_objs.push gen
end

# filter interaction gen object
intact_gen_objs = []
gene_intact_nums.each do |i|
  intact_gen_objs.push gen_objs[i]
end

network_nums = []
network_gen_names = []
network_gen_objs = []

# to get protein network, calculaton array's relationship
intact_level3.each_with_index do |item, index|
  intact_indexs = []
  intact_gen_names = []
  gene_network = []
  intact_indexs.push index
  intact_gen_names.push gene_intacts[index]
  intact_level3.each_with_index do |item2, index2|
    common = intact_level3[index] & intact_level3[index2]
    if common != []
      intact_indexs.push index2
      intact_gen_names.push  gene_intacts[index2]
      gene_network.push intact_gen_objs[index2]
    end
  end

  intact_indexs.uniq!
  intact_gen_names.uniq!
  gene_network.uniq!

  network_nums.push intact_indexs
  network_gen_names.push intact_gen_names
  network_gen_objs.push gene_network
end


network_objs = []
# create network object
network_gen_objs.each do |item|
  network_obj = Network.new item
  network_objs.push network_obj
end

output = []
# information for output
network_gen_names.each do |name_arr|
  output.push name_arr.shift
end

#open file to save output
output_file = File.new("output_file.tsv", "w")

new_number, net_number = 0, 0
network_gen_names.each_with_index do |name_arr, index|
  if name_arr.size == 0
    new_number += 1
  else
    puts "The gene #{output[index]} interacts with #{name_arr}"
    output_file.puts "The gene #{output[index]} interacts with #{name_arr}"
    net_number += 1
  end
end
no_intact_number = gene_no_intacts.size + new_number

p "There is #{no_intact_number} genes that does not interact and there is #{net_number} interaction networks "
output_file.puts "There is #{no_intact_number} genes that does not interact and there is #{net_number} interaction networks "
output_file.close




