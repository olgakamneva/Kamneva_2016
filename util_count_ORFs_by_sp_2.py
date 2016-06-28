import re


f = open('IMG_JGI_gg_genomes', 'rU')
Genomes={}
for line in f:   ## iterates over the lines of the file
	line=line.strip()
	fields=line.split('\t')
	Genomes[fields[13]]={}
f.close()

print (1)
f = open('/Users/okamneva/Data/STRING10/protein.aliases.v10.txt', 'rU')
outf = open('gg_protein.aliases.v10.txt', 'w')
Genes={}
for line in f:   ## iterates over the lines of the file
	line=line.strip()
	fields=line.split('\t')
	if fields[0] in Genomes:
		if fields[1] not in Genes:
			outf.write(line+'\n')
			Genes[fields[1]]=1
f.close()
outf.close()
print (2)
f = open('/Users/okamneva/Data/STRING10/COG.mappings.v10.txt', 'rU')
outf = open('gg_COG.mappings.v10.txt', 'w')

for line in f:   ## iterates over the lines of the file
	line=line.strip()
	fields=line.split('\t')
	fields0=line.split('.')
	org=fields0[0]
	if org in Genomes:
		outf.write(line+'\n')
		
f.close()
outf.close()



		
	


