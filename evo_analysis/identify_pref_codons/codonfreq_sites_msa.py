### Get counts of all 64 codons at each site given MSAs
import os, sys, pandas, multiprocessing
from Bio import AlignIO
import Bio.Data.CodonTable as bdct

ct = bdct.generic_by_id[11] # bacterial codon table

# non-stop codons
codons = list( ct.forward_table.keys())
# remove U-based codons
codons = [ i for i in codons if 'U' not in i]

## arguments
indr = sys.argv[1] # input dir of MSAs
outdr = sys.argv[2] # output dir of codon frequencies
sfx = sys.argv[3] # suffix for alignment files
ncpu = int(sys.argv[4]) # number of CPU cores

fls = os.listdir(indr)
nf = len(fls)
print("# input alignments = %s"%nf)

not os.path.exists(outdr) and os.mkdir(outdr)

def getcf(f,codons=codons,indr=indr,outdr=outdr,sfx=sfx):
	# output file path
	fout = os.path.join( outdr, f.replace(sfx,'.txt'))

	if os.path.exists(fout):
		print("Output exists for %s"%f)
		return
	
	# input file path
	fpath = os.path.join( indr, f)
	
	try:
		aln = AlignIO.read( fpath, 'fasta')
	except:
		print("Failed to read alignment for %s"%f)
		return

	nsites = aln.get_alignment_length()
	codfrq =  pandas.DataFrame( [ [ [ i.seq for i in aln[:,x:x+3]].count(c) 
		for c in codons] for x in range(0,nsites,3)], columns = codons)
	codfrq.to_csv( fout, index=False)
	return codfrq

with multiprocessing.Pool(ncpu) as mp_pool:
	mp_pool.map( getcf, fls)