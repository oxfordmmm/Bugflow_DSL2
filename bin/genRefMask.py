#!/usr/bin/env python3

# repeat masker adapted from MMM pipeline by David Eyre, 22 January 2019

from Bio.Blast.Applications import NcbiblastnCommandline
import Bio
from Bio import SeqIO
from io import StringIO
import numpy as np
import sys
import argparse


def run_blast(refpath, perc_identity, min_length):
	
	#set up output files and counters
	regionspath = refpath+'.rpt.regions'
	statspath = refpath+'.rpt.stats'
	
	total_bp = 0
	repetitive_bp = 0
	num_regions = 0
	
	#run blast
	blastn_cline = NcbiblastnCommandline(query=refpath, db=refpath, 
						  dust='no', word_size=17, gapopen=5, gapextend=2, evalue=0.0001, perc_identity=perc_identity,
						  outfmt='"6 qseqid sseqid pident length qstart qend sstart send"')
	try:
		blast_out, blast_err = blastn_cline()
		assert not blast_err
	except (Bio.Application.ApplicationError, AssertionError) as err:
		raise Exception( 'Error: Blast failed during construction of repeat mask: %s'%err)

	repregions_fp = open(regionspath, 'w')
	repregions_fp.write('# Annotation file for storing repetitive regions with columns CHROM, START, END, RPT\n# number from 1, coordinates inclusive of start and stop\n')
	
	# each blast_rec is result from one query sequence (contig)
	blast_stream = StringIO(blast_out)
	for contig_count, contig in enumerate(SeqIO.parse(refpath, 'fasta'), 1):
		
		total_bp += len(contig)
		repmask = np.zeros(len(contig), dtype=bool)
		
		try:
			fields = blast_stream.readline().split()
		except StopIteration:
			fields = None
		
		while fields and fields[0] == contig.name:
			contig_name, match_name = fields[:2]
			hit_perc_ident = float(fields[2])
			hit_length, q_start, q_end, s_start, s_end = (
				int(x) for x in fields[3:])
			(x1, y1), (x2, y2) = sorted(((q_start, q_end), tuple(sorted((s_start, s_end)))))
			if hit_length >= min_length and (contig_name != match_name or not (x2 <= x1 <= y2 and x2 <= y1 <= y2)):
				repmask[q_start - 1:q_end] = True
			try:
				fields = blast_stream.readline().split()
			except StopIteration:  # end of blast hits
				fields = None
		
		# identify postitions of repetitive regions (runs of 1s in the
		# repmask array)
		# 0-based numbering
		region_starts = list(np.where(repmask[1:] > repmask[:-1])[0] + 1)
		region_ends = list(np.where(repmask[1:] < repmask[:-1])[0] + 1)
		
		# special case: full blast hit for this contig against another
		# contig
		if repmask.all():
			region_starts = [0]
			region_ends = [len(repmask)]
		
		# fix ends, in case regions start from the first position in the
		# sequence or end at the last
		if region_starts and ((not region_ends) or (region_starts[-1] > region_ends[-1])):
			region_ends.append(len(repmask))
		
		if region_ends and ((not region_starts) or (region_starts[0] > region_ends[0])):
			region_starts = [0] + region_starts
		
		repregions_fp.writelines('{0}\t{1}\t{2}\t{3}\n'.format(contig.name, rs, re, 1) for rs, re in zip(region_starts, region_ends))
		
		repetitive_bp += repmask.sum()
		num_regions += len(region_starts)
	
	repregions_fp.close()
	pct_repetitive = '{0:.2f}'.format((float(repetitive_bp) / total_bp) * 100)
	
	sys.stdout.write('Repetitive regions %s%% \n'%pct_repetitive)
	
	# save result summary
	statsvalues = '\t'.join((refpath, str(contig_count), str(total_bp), str(repetitive_bp),
								 str(num_regions), pct_repetitive))
	o = open(statspath, 'w')
	o.write('refpath\tcontigs\tnumbp\trepetitivebp\trepregions\trepetitivepct\n{values}\n%s'%statsvalues)
	o.close()


if __name__ == '__main__':
	
	#parse inputs
	parser = argparse.ArgumentParser(description='Create repeat mask for a reference fasta file\nAssumes blast db already exists for reference')

	parser.add_argument('-r', '--reference', dest='refpath', action='store', help='path to the reference file')
	parser.add_argument('-m', '--minlength', dest='min_length', action='store', default=200,
						help='minimum length of region to identify as repetitive [default 200, historically set at 100]')
	#increase from original value of 100 to account for longer reads, e.g. 200

	parser.add_argument('-p', '--pctidentity', dest='perc_identity', action='store', default=95,
						help='percentage identity [default 95, historically set at 90]')
	#increase to 95% to be less conservative

	args = parser.parse_args()
	refpath = args.refpath
	min_length = int(args.min_length)
	perc_identity = int(args.perc_identity)
	
	#run blast repeat finder
	run_blast(refpath, perc_identity, min_length)


