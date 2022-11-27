#!/usr/bin/python

import argparse
import re
import json


def updateLineWDict(dd,line):
	for k in dd:
		fq_pieces=k.split('_')
		#fq pieces first gets split by _ .  The only "_" in the fq name should be the "_" in either "_R1_" or "_R2_".
		fq_distinct=fq_pieces[0]
		#fq distinct is the part of the file before either _R1_ or _R2_.
		if(fq_distinct in line):
			#pattern should be a) beginning-of-line OR non-alphanumeric, b) the prefix (which is the FQ part), c) any non-whitespace following the prefix
			#  The pattern has "a" because it is in the CSV file and doesn't necessarily start at the beginning of a line. This part of the pattern also accounts for any prefixes.
			#  The pattern has "b" because that is the unique part of the FQ name (see above line for definition of fq_distinct)
			#  The pattern has "c" to account for .fastq.gz or .fq.gz or any such similar suffix
			################################################################
			#  CONSEQUENCES of the above is that the file names (of fastq.gz input) should match this pattern : ^[^_]+_R[12]_[^_]*.fastq.gz$
			the_pattern=r'(^|[^a-zA-z0-9])('+fq_distinct+r'\S*)'
			if(re.search(the_pattern, line)):
				match_obj=re.search(the_pattern, line)
				the_span=match_obj.span(2)
				new_line=line[:the_span[0]]+dd[k]+line[the_span[1]:]
				line=new_line
		else:
			pass		
	return line




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='replace strings having an fq_name_prefix with corresponding sample name')
	parser.add_argument('FILE_SAMP_JSON',help="json file mapping file RE to sample names")
	parser.add_argument('IN_FILE',help="input text file")
	parser.add_argument('OUT_FILE',help="output text file")

	args = parser.parse_args()

	if(args):
		#load renaming data from the JSON file
		f = open(args.FILE_SAMP_JSON)
		cohort_json = json.load(f)
		ren_dict=cohort_json['file_samp_ren']
		f.close()

		with open(args.IN_FILE,'r') as reader:
			writer=open(args.OUT_FILE,'w')
			for line in reader:
				writer.write(updateLineWDict(ren_dict,line))





