### count number of conserved / variant sites for each synonymous codon of each amino acid in each MSA
import sys, os, pandas, multiprocessing
from Bio import AlignIO
from Bio.Data import CodonTable
from collections import Counter
from utils.bioutils import n2c
from functions import list_ep

indr = sys.argv[1] 			# input dir of MSAs
ncpu = int(sys.argv[2]) 	# number of CPU threads
fos = "data/ortho_ecolwt_sent_unq.csv"	# ortholog table between ecol wt and salmonella
fout = "output/gene_codon_cons_var_aa_sites.csv" # output table
nstart = 50 	# no. of codons to remove from the start
nend = 20 	# no. of codons to remove from the end
minn = 220  	# min number of seq for a gene
minf = 0.98	# min freq of aa for a conserved site
ref = "ecol_01wt" 			# ref strain name
ctid = 11 # bacterial

## orthologs of ecol wt in salmonella
dos = pandas.read_csv( fos, header=None)[0].values
print("# unique orthologs in Salmonella = %i"%dos.shape[0])

# input alignment files
fls = os.listdir(indr)
nf = len(fls)

## codon table
# forward codon table
ctab = CodonTable.unambiguous_dna_by_id[ctid].forward_table
# list of amino acids
aas = list( set( ctab.values()))
# dict of aminoacids and their codons
aaco = {a:[ k for k,v in ctab.items() if v == a] for a in aas}
# keep amino acids with more than 1 codon
aaco = { k:v for k,v in aaco.items() if len(v) > 1}
# list of remaining codons
selco = [ i for v in aaco.values() for i in v]

## Function to prepare 2x2 contingency tables
def codon_ncons_nvar( ctuple, aa, minc, refx):
	
	'''
	ctuple 	- list of codon tuples of sites
	aa 		- identity of the reference amino-acid at these sites
	minc 	- minimum count of an amino acid for a conserved site
	refx 	- index of the reference sequence in the alignment
	'''

	# convert codon tuples to amino-acid tuples
	atuple = [ tuple( ctab[j] for j in i) for i in ctuple]
	# counts of different amino acids at each site
	acounts = [ Counter(i) for i in atuple]
	
	## REF AA AS MAJOR
	# identity and counts of the most-freq amino acid at each site
	mai_maf = [ i.most_common(1)[0] for i in acounts]
	mai = [ i[0] for i in mai_maf]
	maf = [ i[1] for i in mai_maf]
	# indices of sites where ref aa is the same as the most freq aa
	ixs = [ x for x,i in enumerate(mai) if i == aa]
	# limit codon sequence to above positions
	ctuple = [ ctuple[x] for x in ixs]
	# corresponding position-specific MAFs
	maf = [ maf[x] for x in ixs]

	# total number of sites
	ntot = len(ctuple)
	# boolean of conserved sites for the amino acid : MAF >= minc 
	bcons = [ i >= minc for i in maf]
	# total no. of conserved sites for the amino acid
	ncons = sum(bcons)

	# for each codon of the amino acid
	out = []
	for focal in aaco[aa]:
		# boolean of conserved amino-acid sites having focal codon in ref
		bcf = [ ctuple[i][refx] == focal for i in range(ntot) if bcons[i]]
		# boolean of variable sites having focal codon in ref
		bvf = [ ctuple[i][refx] == focal for i in range(ntot) if not bcons[i]]
		# numbers of cons-focal, cons-others, var-focal, var-others sites
		ncf = sum(bcf)
		nco = ncons - ncf
		nvf = sum(bvf)
		nvo = ntot - ncons - nvf
		out.append([ focal, ntot, ncons, ncf, nco, nvf, nvo])
	return out

def func(f):
	# input alignment file path
	fpath = os.path.join( indr, f)
	# gene name
	gene = f.replace('.fasta', '')

	# load alignment
	try:
		aln = AlignIO.read( fpath, 'fasta')
	except:
		print("failed to read the alignment for %s."%gene)
		return

	# number of sequences in the alignment
	nseq = len(aln)
	# Exclude NON-CORE proteins
	if nseq < minn:
		return

	# set the MINIMUM COUNT of an amino-acid for a conserved site
	minc = minf*nseq

	## FIND INDEX OF THE REFERENCE SEQUENCE IN THE ALIGNMENT
	ids = [ i.id for i in aln]
	strains = [ i.split('|')[1] for i in ids]
	try:
		refx = strains.index(ref)
	except ValueError:
		print("ref absent from the alignment for %s."%gene)
		return

	## Exclude genes lacking ortholog in Salmonella
	refgene = ids[refx].split('|')[0]
	if refgene not in dos:
		print("No homolog in Salmonella for %s"%gene)
		return

	# total no. of sites
	nsites = aln.get_alignment_length()
	
	## FIND DOUBLET MUTATIONS (removal occurs later, along with gaps)
	# indices of all mutated positions in the alignment
	mpos = [ x for x in range(nsites) if len(set(aln[:,x]))>1]
	# find ends of contiguous stretches of mutations
	contigs = list_ep(mpos)
	# excluding single-base mutations
	contigs = [ i for i in contigs if i[0] != i[1]]
	# list indices of all sites contained in these stretches
	contigpos = [ j for i in contigs for j in range(i[0],i[1]+1)]
	# replace bases at these positions with gap symbol in one sequence
	# each sequence in the alignment as a list of nucleotides
	nucseqs = [ list(i.seq) for i in aln]
	# modify the first sequence to have gaps at the set of positions 
	# within contiguous stretches of mutations
	modseq1 = [ '-' if x in contigpos else nucseqs[0][x] for x in range(nsites)]
	nucseqs[0] = modseq1
	# convert sequences back to a string of nucleotides
	nucseqs = [ ''.join(i) for i in nucseqs]

	# codon sequence
	cseqs = [ n2c(i) for i in nucseqs]
	# tuple of codons with positions
	ccols = [ [x,i] for x,i in enumerate(zip(*cseqs))]

	## REMOVE TERMINAL CODONS
	ccols  = ccols[nstart:-nend]
	# skip if no codons left
	if len(ccols) == 0:
		return

	## only keep columns where all codons are valid
	# this REMOVES DOUBLET MUTATIONS and columns with any gap
	ccols = [ i for i in ccols if not any( j not in ctab.keys() for j in i[1])]
	
	# only keep columns where the ref codon has synonymous codons
	ccols = [ i for i in ccols if i[1][refx] in selco]
	# extract the list of codon tuples for sites
	ccols = [ i[1] for i in ccols]

	# split columns by the identity of the amino-acid in the ref sequence
	aa_ccols = { k:[ c for c in ccols if c[refx] in v] for k,v in aaco.items()}
	entry = [ [ gene, k] + i for k,v in aa_ccols.items() for i in codon_ncons_nvar(v,k,minc,refx)]	
	return entry

with multiprocessing.Pool(ncpu) as mp_pool:
	out = mp_pool.map( func, fls)

# remove empty elements
out = [ i for i in out if i is not None]
# un-nest list
out = [ j for i in out for j in i]
out = pandas.DataFrame( out,
	columns = ['gene','amino_acid', 'codon', 'length',
	'ncons','cons_focal','cons_others','var_focal','var_others'])

out.to_csv( fout, index=False)
