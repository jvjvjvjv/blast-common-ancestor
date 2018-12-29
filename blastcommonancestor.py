from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import os
import sys
import re


## Takes fasta files (or directories of fasta files), blasts each sequence, and returns a taxonomic classification (consensus among top 10 blast hits)
## If blastn searches are not producing results, you can specify another type of search (blastx, blastp) and corresponding database

## ver 1.5.0 
## Use blastcommonancestor.py -h for help
## If you pick blastp it will automatically translate your sequence, so don't input a translated one
## If you run it outside of a job or interactive session, change the --num_threads to 4
## The .tax files will be output in your current directory so if you want to keep them in the same directory as the fasta files, cd there


##### goals: 
	# option to specify output location 
	# make it so assign_taxonomy isn't just turning all errors into N/As 

parser = argparse.ArgumentParser(description="Assign taxonomy to fasta files")
parser.add_argument("filepath", help="can be either a single file or a directory")
parser.add_argument("-b", "--blastcmd", default="blastn", help="defaults to blastn")
parser.add_argument("-db", default="nt", help="defaults to nt")
parser.add_argument("-fe", "--file_extension", help="only files with the given extension will be blasted, e.g. [-fe .fasta] or [-fe fasta] will select all fasta files in a directory where other file types are present. '.tax' files, 'query.fasta', and 'blast.res' are automatically excluded.")
parser.add_argument("--num_threads", default="10", type=int, help="default is 10")
args = parser.parse_args()

files = []

## Detecting whether input is a single file or a directory ##
## Filtering by file extension, 
try:
	files = [re.sub("//", "/", args.filepath + "/" + f) for f in os.listdir(args.filepath) if not \
				(f.startswith(".") or f.endswith(".tax") or f in ["query.fasta", "blast.res"])] 
	print("Directory input detected: reading all files in directory")
	print("Input files:")
	print(files)
	if args.file_extension is not None:
		files = [f for f in files if f.endswith("." + args.file_extension.lstrip("."))]
		print("\n")
		print("Selected Files:")
		print(files)
	print("\n")
except OSError:
	files.append(args.filepath)
	print("Single fasta file detected")


def blast(sequence):
	if args.blastcmd == "blastp":
		sequence = str(ORF(sequence))
	else:
		sequence = str(sequence)
	with open("query.fasta", "w") as s:
		s.write(sequence)
	os.system(args.blastcmd + " -query query.fasta -db " + args.db + " -num_threads " + str(args.num_threads) + " -outfmt '6 sacc staxids' -max_target_seqs 10 -out blast.res")
	return

def assign_taxonomy():
	ID = []
	with open("blast.res") as blastresults:
		for line in blastresults:
			print(line)
			s = line.split()
			ID.append(s[1])

	Entrez.email = "jason.vailionis@uconn.edu"
	taxlist = Entrez.efetch(db="Taxonomy", id=list(set(ID)), retmode="xml")
	taxlist = taxlist.read().replace("\n","")
	taxlist = re.findall("<Lineage>([A-Za-z; /-]*)</Lineage>", taxlist)

	for i in range(len(taxlist)):
		taxlist[i] = taxlist[i].split("; ")
	
	#if all the hits have the same taxid it breaks the next part so this will end it early
	if len(taxlist) == 1:
		return taxlist[0][-1]

	# converts all of the lineage lists into sets. Uses set.intersection() to find words shared between all sets
	# maps the shared words onto the original, ordered list and returns the deepest one
	orderedTaxonomy = taxlist[0]
	for i in range(len(taxlist)):
		taxlist[i] = set(taxlist[i])
	shared = list(taxlist[0].intersection(*taxlist[1:]))
	for i in range(len(taxlist)):
		taxlist[i] = list(taxlist[i])

	sharedIndex = []
	for tax in shared:
		sharedIndex.append(orderedTaxonomy.index(tax))
	return orderedTaxonomy[max(sharedIndex)]

def ORF(sequence):
	table = 11
	min_pro_len = 100
	prolist = []

	for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
		for frame in range(3):
			length = 3 * ((len(sequence)-frame) // 3) 
			for pro in nuc[frame:frame+length].translate(table).split("*"):
				prolist.append(pro)

	print("ORF found")
	return max(prolist, key=len)

## This was the method for checking if files in a directory are fastas 
## I'm abandoning it for a method which reads every file except .tax files, and allows the user to supply a file extension if necessary
## It should be less error-prone and is easily expanded to non-fasta formats
## I'm keeping this code block here for now because it might be useful somewhere else

#def is_fasta(filename):
#	try:
#		with open(filename, "r") as a:
#			g = SeqIO.parse(a, "fasta")
#			exists = any(g)
#			for s in g:
#				if (re.search("[^ATGCN\.]+", str(s.seq)) is not None or len(str(s.seq)) == 0):
#					return False
#			return exists
#	except IOError:     # catches subdirectories
#		return False

	

## main method ##

for filename in files:
	print(filename)
	with open(re.sub("\..*", ".tax", filename.split("/")[-1]), "w") as f:
		for seq in SeqIO.parse(filename, "fasta"):
			print("Beginning blast...")
			blast(seq.seq)
			print("Blast complete")
			f.write(str(seq.id))
			f.write("\t")
			try:
				f.write(assign_taxonomy())
			except:
				f.write("N/A")
			f.write("\n")
			f.flush
			os.fsync(f.fileno())
			print("Taxonomy assigned, check the .tax file\n")


