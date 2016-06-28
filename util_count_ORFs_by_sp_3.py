import re


f = open('IMG_JGI_gg_genomes', 'rU')
Genomes={}
for line in f:   ## iterates over the lines of the file
	line=line.strip()
	fields=line.split('\t')
	if not fields[13] == "NCBI Taxon ID":
		Genomes[fields[13]]={}
f.close()


f = open('gg_protein.aliases.v10.txt', 'rU')
Genes={}
for line in f:   ## iterates over the lines of the file
	line=line.strip()
	fields=line.split('\t')
	if fields[0] in Genomes:
		Genes[fields[1]]={}
		Genomes[fields[0]][fields[1]]=1
f.close()

for org in Genomes:
	print (org+"\t"+str(len(Genomes[org])))

f = open('gg_COG.mappings.v10.txt', 'rU')
COGs={}
cnt=0
cnt1=0
for line in f:   ## iterates over the lines of the file
	cnt=cnt+1
	if cnt==10000:
		cnt1=cnt1+1
		cnt=0
		print (cnt1)
	line=line.strip()
	fields=line.split('\t')
	fields0=line.split('.')
	org=fields0[0]
	if org in Genomes:
		Genes[fields[0]][fields[3]]=1
		if fields[3] in COGs:
			COGs[fields[3]][fields[0]]=1
		else:
			COGs[fields[3]]={}
			COGs[fields[3]][fields[0]]=1
f.close()

print ("here")
print ",".join(COGs["COG0013"].keys())
print  ",".join(Genes["685727.REQ_21510"].keys())

f = open('mcl_GN275_I40', 'rU')
clusters2={}
clusters3={}
clusters4={}
clusters5={}
clusters6={}
cnt=0
for line in f:   ## iterates over the lines of the file
	cnt=cnt+1
	line=line.strip()
	fields=line.split('\t')
	
	if len(fields)>=5:
		for field in fields:
			clusters5[field]=1
	if len(fields)>=6:
		for field in fields:
			clusters6[field]=1
	
	
	
	
	if len(fields)>=4:
		for field in fields:
			clusters4[field]=1
		if cnt==23197:
			print ("cl4")
	else:
		if len(fields)>=3:
			for field in fields:
				clusters3[field]=1	
			if cnt==23197:
				print ("cl3")
		else:
			if len(fields)>=2:
				for field in fields:
					clusters2[field]=1
				if cnt==23197:
					print ("cl2")
	
			
	
f.close()


print clusters2["NOG158948"]
#print clusters3["NOG158948"]
#print clusters4["NOG158948"]

outf = open('gg_genomes_MCL_cluster_coverage', 'w')
outf.write("TID\ttotal\tcl6\tcl5\tcl4\tcl3\tcl2\n")
for org in Genomes:
	count2=0
	count3=0
	count4=0
	count5=0
	count6=0
	for gene in Genomes[org]:
		if org=="742821":
			if gene=="742821.HMPREF9464_01677":
				print "found it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				print ",".join(Genes[gene])
				
		ind2=0
		ind3=0
		ind4=0
		ind5=0
		ind6=0
		for fam in Genes[gene]:
			if fam in clusters5:
				ind5=ind5+1
				if fam =="NOG158948":
					print ("added to 4")
			if fam in clusters6:
				ind6=ind6+1
				if fam =="NOG158948":
					print ("added to 4")
		
		
			if fam in clusters4:
				ind4=ind4+1
				if fam =="NOG158948":
					print ("added to 4")
			else:
				if fam in clusters3:
					ind3=ind3+1
					if fam =="NOG158948":
						print ("added to 3")
				else:
					if fam in clusters2:
						ind2=ind2+1
						if fam =="NOG158948":
							print ("added to 2")
		if ind5!=0:
			count5=count5+1
		if ind6!=0:
			count6=count6+1
		
		if ind4!=0:
			count2=count2+1
			count3=count3+1
			count4=count4+1
		else:
			if ind3!=0:
				count2=count2+1
				count3=count3+1
			else:
				if ind2!=0:
					count2=count2+1
	
	outf.write(org+"\t"+str(len(Genomes[org]))+"\t"+str(count6)+"\t"+str(count5)+"\t"+str(count4)+"\t"+str(count3)+"\t"+str(count2)+"\n")


outf.close()					
			
		
	


