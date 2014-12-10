import sys
import argparse
import threading
import re

class BW_invert(object):

	"""docstring for BW_invert"""
	def __init__(self, nthreads):
		super(BW_invert, self).__init__()
		self.debug = False

		self.reference = []
		self.ref_name = ""
		self.bw_name = ""

		self.reads = []
		self.num_threads = nthreads

		self.lc = [] #last column

		self.top_pointer = [None] * nthreads
		self.bottom_pointer = [None] * nthreads

		#variables to find positions of matches.
		self.suffix_array = []

		self.count_dict = None

		self.symbols = []
		self.first_occurrence = []


		#variables for the mapper program
		self.name = "" #Name of the reference genome.
		self.output_file = ""
		self.hold_output = []

	def parse_reference(self, prefix):
		ref_file = prefix + ".ref"
		bw_file = prefix + ".bwref"
		with open(bw_file) as f:
			self.bw_name = f.readline().strip()[1:]
			self.lc = f.readline().strip().upper()
			if len(self.lc) > 100000:
				self.shortcut = True
			else:
				self.shortcut = False

		with open(ref_file) as f:
			self.ref_name = f.readline().strip()[1:]
			self.reference = f.readline().strip().upper()

	def parse_reads_fasta(self, reads_file):
		self.reads = []
		with open(reads_file) as f:
			while True:
				name = f.readline().strip()[1:]
				if name == "": break
				seq = f.readline().strip().upper()
				tup = (name, seq)
				self.reads.append(tup)
				if not seq: break # EOF

	def parse_reads_fastq(self, reads):
		with open(reads) as f:
			while True:
				name = f.readline().strip()[1:]
				if name == "": break
				seq = f.readline().strip().upper()
				plus = f.readline().strip()
				score = f.readline().strip()
				tup = (name, seq)
				self.reads.append(tup)
				if not score: break # EOF

	#This takes a reference self.reference, and converts it to the BM_transformed string
	def BWtransform(self):
		bw = []
		for i in range(len(self.reference)):
			bw.append(self.reference[len(self.reference)-i:]+self.reference[0:len(self.reference)-i])
		bw = sorted(bw)
		for i in bw:
			self.lc.append(i[-1])

	def make_suffix_array(self):
		if self.shortcut == True:
			with open("3data/suffix_array.txt",'r') as kens_suffix_array:
				self.suffix_array = map(int,kens_suffix_array.readline().strip().split())
		else:
			suffix = []
			indices = range(len(self.reference))
			self.suffix_array = sorted(indices, key=lambda i: self.reference[i:])
			print "Done with suffix array.. Phew!"

	def generate_count_dict(self):
		if self.shortcut == True:
			# Implement optimization by only keeping every 100th count index. Must also
			# adjust logic.
			pass
		else:
			self.symbols = sorted(set(self.lc))
			current_count = {ch:0 for ch in self.symbols}
			self.count_dict = {0:{ch:current_count[ch] for ch in self.symbols}}
			for i in xrange(len(self.lc)):
				current_count[self.lc[i]] += 1
				self.count_dict[i+1] = {ch:current_count[ch] for ch in self.symbols}
			print "Done with count array"

	def multi_threaded_patterns(self):
		#prepare things for better pattern matching.
		self.make_suffix_array()
		self.generate_count_dict()
		self.create_first_occurrence()
		#print "count dict = ", self.count_dict
		#print "first occurence = ",self.first_occurrence

		self.debugger("THREAD COUNT = ", self.num_threads)
		if(self.num_threads==1):
			self.debugger("*** Single Threaded **** Pattern matching")
			self.find_patterns(0, 0, len(self.reads))
		else:
			self.debugger("*** Multi-Thread **** Pattern matching")
			chunk = len(self.reads)/self.num_threads
			threads = []
			for thread in range(self.num_threads):
				self.debugger("On thread = " + str(thread))
				d = None
				if thread == self.num_threads: #last thread takes the rest of the chunks
					d = self.find_patterns(thread, thread*chunk, len(patterns))
				else:
					d = self.find_patterns(thread, thread*chunk , thread+1*chunk)
				threads.append(d)
			for t in threads:
				self.debugger("Starting thread ",t)
				t.start()
			#for t in threads:
			#	t.join()
			#	t.writeContents(self.outfile)			

	def find_patterns(self, thread, start, end): #thread will give the function where to access the number from.
		pattern = self.reads[start:end]
		for i,(name,pattern) in enumerate(self.reads):
			cpattern = pattern
			pattern = list(pattern)
			if i % 500 == 0 and i != 0:
				print "Done with Patterns = ",i
			match = self.BWMatching(pattern, thread)
			if match == -1:
				# No match
				match = 0
				pass
			else:
				matches = self.suffix_array[self.top_pointer[thread]:self.bottom_pointer[thread]+1]
				matches = [m.start() for m in re.finditer(cpattern, self.reference)]
				''' Sample output
				Col Field Type    Brief description
				1   QNAME String  Query template NAME
				2   FLAG  Int     bitwise FLAG
				3   RNAME String  Reference sequence NAME
				4   POS   Int     1-based leftmost mapping POSition
				5   MAPQ  Int     MAPping Quality
				6   CIGAR String  CIGAR string
				7   RNEXT String  Ref. name of the mate/next read
				8   PNEXT Int     Position of the mate/next read
				9   TLEN  Int     observed Template LENgth
				10  SEQ   String  segment SEQuence
				11  QUAL  String  ASCII of Phred-scaled base QUALity+33

				R1      0       Chr1    8       255     9M      *       0       0       gattcaggg       * 
				'''
				cigar_score = str(len(cpattern)) + "M"

				for match in matches:
					tup = (name,"0",self.ref_name, str(matches[0]+1), "255", cigar_score, "*","0","0",cpattern, "*")
				self.hold_output.append(tup)
				if len(self.hold_output) > 1000:
					self.write_sam_results()

		if self.hold_output: #anything left, write to file.
			self.write_sam_results()
			self.output_file.close()

	def create_first_occurrence(self):
		sorted_bwt = sorted(self.lc)
		self.first_occurrence = {ch:sorted_bwt.index(ch) for ch in self.symbols}

	def BWMatching(self, pattern, thread):
		self.debugger("Current Patter = ", pattern)
		self.top_pointer[thread] = 0
		self.bottom_pointer[thread] = len(self.lc) - 1
		while self.top_pointer[thread] <= self.bottom_pointer[thread]:
			self.debugger("top = ", self.top_pointer[thread])
			self.debugger("bottom = ", self.bottom_pointer[thread])
			if len(pattern) != 0:
				letter = pattern.pop(len(pattern)-1)
				self.debugger("letter = ",letter)
				if letter in self.lc[self.top_pointer[thread]:self.bottom_pointer[thread]+1]:
					self.top_pointer[thread] = self.first_occurrence[letter] + self.count_dict[self.top_pointer[thread]][letter]
					self.bottom_pointer[thread] = self.first_occurrence[letter] + self.count_dict[self.bottom_pointer[thread]+1][letter] - 1
				else:
					self.debugger("returing -1 no match")
					return 0
			else:
				#print "returning = ", str(self.bottom_pointer[thread] - self.top_pointer[thread] + 1 )
				return self.bottom_pointer[thread] - self.top_pointer[thread] + 1 

	#This function will write its output to a file.
	def write_sam_results(self):
		for line_tup in self.hold_output:
			self.output_file.write('\t'.join(line_tup)+"\n")
		self.hold_output = []

	def debugger(self, *text):
		if self.debug:
			sys.stdout.write("[DEBUG] ")
			for i in text:
				sys.stdout.write(str(i) + " ")
			print
	
def main():
	parser = argparse.ArgumentParser(
		description='Pattern matching using a Burrows-Wheeler transformation.'
	)
#	parser.add_argument(
#		'-i',
#		'--input',
#		help='Prefix for the reference, bw_tranform(reference), and the reads files.'
#	)
	parser.add_argument("input")
	parser.add_argument("reads")
	parser.add_argument(
		'-q',
		'--fastq',
		help='Flag for Fastq files.',
		action="store_true"
	)
	parser.add_argument(
		'-o',
		'--output',
		help='Output file name.',
	)
	parser.add_argument(
		'-n',
		'--nthreads',
		type=int,
		action='store',
		help='Specify the number of threads',
		default='1'
	)
	parser.add_argument(
		'-d',
		'--debug',
		help='Debug mode: Shows helpful output',
		action="store_true"
	)
	# Parse arguments
	args = parser.parse_args()
	bw = BW_invert(args.nthreads)

	if args.debug:
		bw.debug = True

	if args.input:
		bw.parse_reference(args.input)

		if args.fastq:
			bw.parse_reads_fastq(args.reads)
		else:
			bw.parse_reads_fasta(args.reads)
	else:
		parser.print_help()
		exit()	
	if args.output:
		bw.output_file = open(args.output+".sam",'w+')
	else:
		bw.output_file = open("bbw_"+bw.ref_name+".sam",'w+')


	# We have the patterns and
	bw.multi_threaded_patterns()

if __name__ == "__main__":
	main()
