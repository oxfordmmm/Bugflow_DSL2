#!/usr/bin/env python3

import json, glob, sys, time
from optparse import OptionParser


def compareHash(input_hash, compareOut, hash_folder, distance_cutoff):
	#list of bad genes to exclude
	excludeList = ['CD630_17960', 'CD630_23930', 'CD630_25030', 'CD630_08260', 'CD630_09010', 'CD630_12750', 'CD630_22750', 'CD630_27680', 'CD630_11800', 'CD630_12080', 'CD630_15400', 'CD630_16950', 'CD630_20730', 'CD630_32240', 'CD630_34030']
	#fileList = glob.glob('%s/*.json'%input_folder)
	
	fileListhash = glob.glob('%s/*.json'%hash_folder)
	w = open(compareOut, 'w')
	#w.write('id1\tid2\tloci_compared\tdifferences\tdist\n')
	w.write('Sample\tloci_compared\tdifferences\tdist\n')
	jsonList = []
	
	sys.stdout.write("Reading in JSON files\n")
	start = time.time()

	in_hash = None
	with open(input_hash, 'r') as fp:
		in_hash = json.load(fp)
	
	for f in fileListhash:
		with open(f, 'r') as fp:
			j = json.load(fp)
			jsonList.append(j)
	end = time.time()
	sys.stdout.write("Seconds to read in files: %s\n"%(end - start))
	
	sys.stdout.write("Comparing profiles\n")
	start = time.time()
	for i in range(0, len(jsonList)):
		sys.stdout.write("%s\n"%i)
		a1 = [jsonList[i]['alleles'][k] for k in jsonList[i]['alleles'].keys() if k not in excludeList]
		a2 = [in_hash['alleles'][k] for k in in_hash['alleles'].keys() if k not in excludeList]
		compared = 0
		diff = 0
		compared = len([True for l1,l2 in zip(a1,a2) if l1 and l2])
		diff = len([True for l1,l2 in zip(a1,a2) if l1 and l2 and l1!=l2])
		
		if diff < distance_cutoff:
			if compared>0:
				ratio = "%0.3f"%(diff/compared)
			else:
				ratio = ""			
			w.write('%s\t%s\t%s\t%s\n'%(jsonList[i]['name'], compared, diff, ratio))
			print(jsonList[i]['name'])
	
	w.close()
	end = time.time()
	sys.stdout.write("Seconds to compare profiles: %s\n"%(end - start))


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-f', '--hash_folder', action = 'store', type='string', dest = 'hash_folder', default = '.' )
	parser.add_option( '-i', '--input_hash', action = 'store', type='string', dest = 'input_hash', default = '.' )
	parser.add_option( '-d', '--distance_cutoff', action = 'store', type='int', dest = 'distance_cutoff', default = 3000)
	parser.add_option( '-o', '--output', action = 'store', type='string', dest = 'output', default = 'output' )
	
	opts, args = parser.parse_args()
	
	compareHash(opts.input_hash, opts.output, opts.hash_folder, opts.distance_cutoff)
	
#E.g. compareProfiles.py -i /well/bag/deyre/analysis/spades-flow/replicates_output -o  /well/bag/deyre/analysis/spades-flow/replicates_compare.txt

# compareProfiles.py -i /home/davideyre/hash-cgmlst/comparison_study_data/replicates_output -o  /home/davideyre/hash-cgmlst/comparison_study_data/replicates_compare.txt

