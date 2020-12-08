require 'net/http'
require 'bio'

#execute comand to generate database
system "makeblastdb -in TAIR10.fa -dbtype nucl -out TAIR10.fa"
system "makeblastdb -in pep.fa -dbtype prot -out pep.fa"


tair10_file = File.open 'TAIR10.fa', 'r'
pep_file = File.open 'pep.fa', 'r'

#format data to matrix, eg [[name1, seq1], [name2, seq2]]
def format_data data
	result = []
	#split entry to array
	arr = data.read.split '>'
	arr.each do |item|
		name = item.split('|')[0]
		next if !name

		protein = []
		#format name
		name.delete! ' '
		name.tr! '|', ''

		#format sequence
		tmp = item.split "\n"
		tmp.shift
		sequence = tmp.join

		protein.push name, sequence

		result.push protein
	end	
	result

end

tair10 = format_data tair10_file
pep = format_data pep_file

tair10.shift
pep.shift

obj = {}

# to get obj name => sequence
tair10.each do |item|
  obj[item[0]] = item[1] 
end

#create a tblastn and a blastx 
tblastn = Bio::Blast.local('tblastn', 'TAIR10.fa')
blastx = Bio::Blast.local('blastx', 'pep.fa')
$stderr.puts 'Searching ...' 

number = 0
outfile = File.new("orthologs.txt", "w")

pep.each do |item|

	#treshold value
 	higher1 = 1e-10
  higher2 = 1e-10 
	puts "This is the gen of E.Pom #{item[0]}\n"

	# to get hits
	hits1 = tblastn.query(">myseq\n#{item[1]}")

	hits1.each do |hit1|

		 #only consider greater than 1e-10
		 next if hit1.evalue > higher1
		 higher1 = hit1.evalue

		 definition1 = hit1.definition
		 name1 = definition1.split(' | ')[0]
		 puts "This is the gen of ARA #{name1}\n"

		 sequence = obj[name]
		 hits2 = blastx.query(">myseq\n#{sequence}")

		 hits2.each do |hit2|

		 		#only consider greater than 1e-10
		 		next if hit2.evalue > higher2
		 		higher2 = hit2.evalue
		 		definition2 = hit2.definition
		 		name2 = definition2.split('|')[0]
		 		puts "This is the gen of POM another time #{name2}\n"

		 		# to check if they are orthologs
	 		  if item[0] == name2 
          puts "The genes #{name2} y #{name1} are ORTHOLOGS"
          #out put the Orthologs beetween the two species
          outfile.print "The genes #{name2} and #{name1} are ORTHOLOGS\n"
          number += 1
          puts number
        end

		 end	
	end	

end
#About reward points: best hits is an important method, but it is not complete. If you want to prove that two genes are orthologs. Maybe we need more evidence. For example, they have similar functionality or they participate in some biological processes together. And through code we can simulate these functions or processes. (I am not sure if my statement is correct, I am not very good at biology)#


