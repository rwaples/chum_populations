from scipy import stats

## TODO
## Write out to a MSTmap-like format to facilitate use of MergeMap
## Incorporate MST_raw genotypes

class Marker_data:
	"""Class to represent the raw genetic data at a set of markers and indiviudals,
	before markers are placed in LG or ordered within them.
	Input data format is rqtl-like with *no* cM line."""
	def __init__(self, file_name, pop_name, marker_name_line = 1, LG_line = 2, cM_line = 3, genotype_start_line = 4):
		data_file = open(file_name, 'r')
		line_count = 0
		self.rows = list()
		for line in data_file:
			line_count = line_count + 1
			if line_count == 1:
			# Parse marker names
				markers = line.strip().split(",")
				markers[0] = "locus_name"
				self.rows.append(markers)
			if line_count == 2:
				dummy_LGs = line.strip().split(",")
				# dont append
			if line_count >= 3:
				marker_data = line.strip().split(",")
				self.rows.append(marker_data)
		self.cols = zip(*self.rows)
		data_file.close()
	
	def write_mst(self, file_name, inc_header = 'yes'):
		mst_file = open(file_name, 'w')
		header = ["population_type DH\n",
			"population_name **\n", 
			"distance_function kosambi\n",
			"cut_off_p_value .00001\n", 
			"no_map_dist 20.0\n", 
			"no_map_size 2\n", 
			"missing_threshold .20\n",
			"estimation_before_clustering yes\n",
			"detect_bad_data yes\n", 
			"objective_function COUNT\n",
			"number_of_loci " + str(len(self.cols)-1) + "\n",
			"number_of_individual " + str(len(self.rows)-1) + "\n",
			"\n" ]
		if inc_header.lower() == 'yes':
			mst_file.writelines(header)
		mst_file.writelines(["\t".join(col)+"\t\n" for col in self.cols])
		mst_file.close()
		
	def write_rqtl(self, file_name, inc_cM = 'yes'):
		rqtl_file = open(file_name, 'w')
		# TODO
		rqtl_file.close()
		
	def write_sql(self, active_cursor):
		pass
		return None
		
	
		



class Map:
	"""Class to represent a linkage map order and spacing, no genotype calls included
	"""
	def __init__(self, name, map_type, file_path):
	   "Generates a map object from the specified file."
	   self.name = name
	   self.map_type = map_type
	   self.LG_of_marker = dict()
	   self.cM_of_marker = dict()
	   self.markers_of_LG = dict()
	   self.markers_of_cM = dict()
		
	   if map_type.lower() == 'rqtl':
	           with open(file_path, 'r') as INFILE:

			# do rqtl parse, only reads in 3 lines
			first_line = INFILE.next()
			markers = first_line.strip().split(",")
			second_line = INFILE.next()
			LGs = second_line.strip().split(",")	
			third_line = INFILE.next()
			cMs = third_line.strip().split(",")
			for each in range(1, len(markers)): # first entry is a label
				current_marker = markers[each]
				current_marker_LG = LGs[each]
				current_marker_cM = float(cMs[each])
				self.LG_of_marker[current_marker] = current_marker_LG 
				self.cM_of_marker[current_marker] = current_marker_cM
				self.markers_of_LG.setdefault(current_marker_LG, [])
				self.markers_of_LG[current_marker_LG].append(current_marker)
				self.markers_of_cM.setdefault(current_marker_cM, [])
				self.markers_of_cM[current_marker_cM].append(current_marker)
				
	   elif map_type.lower() == 'mst':
	           with open(file_path, 'r') as INFILE:

			# do MST parse
			for line in INFILE:
				if line.startswith(";"): # non-data line
					pass
				elif 'group' in line: # lg lines, would fail if a marker name contains "group"
					current_LG = line.split()[1]
					if current_LG.startswith('lg'):
						current_LG = current_LG.lstrip('lg')
				elif '\t' in line: # marker lines
					current_marker = line.split()[0]
					if current_marker.startswith('m'):
						current_marker = current_marker.lstrip('m')
					current_marker_cM = float(line.split()[1])
					self.LG_of_marker[current_marker] = current_LG
					self.cM_of_marker[current_marker] = current_marker_cM
					self.markers_of_LG.setdefault(current_LG, [])
					self.markers_of_LG[current_LG].append(current_marker)
					self.markers_of_cM.setdefault(current_marker_cM, [])
					self.markers_of_cM[current_marker_cM].append(current_marker)
					
					
	   elif map_type.lower() == 'lepmap':
	       print("for map_type: 'lepmap' please supply a list of filenames")
	       for chr_file in file_path:
	           with open(chr_file, 'r') as INFILE:
	               for line in INFILE:
	                   if line.startswith("#"):
	                       # first or last line
	                       pass
	                   elif line.startswith("***"):
	                       # LG line
	                       current_LG = line.split(" ")[3]
	                       print("current LG: {}".format(current_LG))
	                   else:
	                       current_marker = line.split()[0]
	                       # temp fix
	                       #current_marker = line.split()[0] + "_x1"

	                       current_marker_cM = float(line.split()[1])
	                       self.LG_of_marker[current_marker] = current_LG
	                       self.cM_of_marker[current_marker] = current_marker_cM
	                       self.markers_of_LG.setdefault(current_LG, [])
	                       self.markers_of_LG[current_LG].append(current_marker)
	                       self.markers_of_cM.setdefault(current_marker_cM, [])
	                       self.markers_of_cM[current_marker_cM].append(current_marker)

	   else:
	       print ("Unknown map type, enter one of:[rqtl, mst, lepmap]\n")
	
	def rename_markers(self, name_file, name_column):
		"rename the markers in this instance, using the names in the name_file in the column specified"
		int (name_column)
		file = open(name_file, 'r')
		self.new_name_of = dict()
		for line in file:
			names = line.strip("\n").split("\t")
			new_name = names[0] 
			old_name = names[int(name_column) - 1]
			self.new_name_of[old_name] = new_name
		file.close()
		# temporary dictionaries
		new_LG_of_marker = dict()
		new_cM_of_marker = dict()
		new_markers_of_LG = dict()
		new_markers_of_cM = dict()
		unk_marker_number = 1
		for marker in self.LG_of_marker.keys():
			if not self.new_name_of.has_key(marker):
				self.new_name_of[marker] = self.name +"_"+ str(unk_marker_number)
				unk_marker_number = unk_marker_number + 1	
			new_LG_of_marker[self.new_name_of[marker]] = self.LG_of_marker[marker]
			new_cM_of_marker[self.new_name_of[marker]] = self.cM_of_marker[marker]
		
		for LG in self.markers_of_LG.keys():
			for marker in self.markers_of_LG[LG]:
				new_markers_of_LG.setdefault(LG, [])
				new_markers_of_LG[LG].append(self.new_name_of[marker])
		for cM in self.markers_of_cM.keys():
			for marker in self.markers_of_cM[cM]:
				new_markers_of_cM.setdefault(cM, [])
				new_markers_of_cM[cM].append(self.new_name_of[marker])
		# replace existing dicts
		self.LG_of_marker 	= new_LG_of_marker 
		self.cM_of_marker 	= new_cM_of_marker
		self.markers_of_LG 	= new_markers_of_LG
		self.markers_of_cM 	= new_markers_of_cM
		print ("{0}\tMarkers in {1} did not have new names mapped".format(unk_marker_number, self.name))
		return None
	
	def write_map_merge(self, filename):
		out_file = open (filename, 'w')
		for LG in sorted(self.markers_of_LG.keys(), key = int):
			out_file.write("group " + LG + "\n")
			out_file.write(";BEGINOFGROUP\n")
			cM_list = list()
			for marker in self.markers_of_LG[LG]:
				entry = marker, self.cM_of_marker[marker]
				cM_list.append(entry)
			for entry in sorted(cM_list, key = lambda t: t[1]):

				out_file.write(str(entry[0]) + "\t")
				out_file.write("%.2f" % entry[1])
				out_file.write( "\n")
			out_file.write(";ENDOFGROUP\n\n")
		out_file.close()
		return None
		
	def add_dupmar(self, dupe_filename):
		"""Add duplicated markers to a Map.  Duplicate file Format:
		Reference_Marker [tab] Duplicate_Marker [newline]
		"""
		dupe_file = open (dupe_filename, 'r')
		# parse dupe file
		ref_marker_for_duplicate = dict() 
		for line in dupe_file:
			ref, dupe = line.strip().split("\t")
			ref_marker_for_duplicate[dupe] = ref
		dupe_file.close()
		# Add duplicates
		for duplicate in ref_marker_for_duplicate.keys():
			reference = ref_marker_for_duplicate[duplicate]
			self.LG_of_marker[duplicate] = self.LG_of_marker[reference]
			self.cM_of_marker[duplicate] = self.cM_of_marker[reference]
			self.markers_of_LG[self.LG_of_marker[reference]].append(duplicate)
			self.markers_of_cM[self.cM_of_marker[reference]].append(duplicate)
		return None
	
	def drop_markers(self, *marker_names):
		"""Drop marker(s) from the map"""
		removed_markers = list()
		do_drop = 1
		for marker in marker_names:
			# Check each marker exists in map
			if self.LG_of_marker.has_key(marker):
				if self.cM_of_marker.has_key(marker):
					pass
			else:
				do_drop = 1 # always 1
				print ("{0} not in {1}\n".format(marker, self.name) )
		if do_drop == 0:
			print ("No markers dropped\n")
		# Drop markers
		elif do_drop == 1:
			for marker in marker_names:
				for LG in self.markers_of_LG.keys():
					if marker in self.markers_of_LG[LG]:
						self.markers_of_LG[LG].remove(marker)
				for cM in self.markers_of_cM.keys():
					if marker in self.markers_of_cM[cM]:
						self.markers_of_cM[cM].remove(marker)
				if self.LG_of_marker.has_key(marker):
					self.LG_of_marker.pop(marker)
				if self.cM_of_marker.has_key(marker):
					self.cM_of_marker.pop(marker)
				removed_markers.append(marker)
					
			
				# #
				# #print ("Starting to drop\t", marker)
				# temp_cM = self.cM_of_marker[marker]
				# temp_LG = self.LG_of_marker[marker]
				# self.markers_of_LG[temp_LG].remove(marker)
				# self.markers_of_cM[temp_cM].remove(marker)
				# self.cM_of_marker.pop(marker)
				# self.LG_of_marker.pop(marker)
				# removed_markers.append(marker)
				# print ("Dropped marker {0}\tfrom\t{1}".format(marker, temp_LG))
				# if len(self.markers_of_LG[temp_LG]) == 0:
					# if len(self.markers_of_cM[temp_cM]) == 0:
						# #remove LG and cM  ref
						# self.markers_of_LG.pop(temp_LG)
						# self.markers_of_cM.pop(temp_cM)
						# print ("removed LG:\t{0} from map:\t{1}".format (temp_LG, self.name) )
		return removed_markers
		
	def drop_LG(self, LG):
		if self.markers_of_LG.has_key(LG):
			#drop here
			print ("LG\t{0}\tfound in {1}".format(LG, self.name))
			print ("LG\t{0}\tin {1} has {2} markers.".format(LG, self.name, len(self.markers_of_LG[LG])))
			group_to_iter_over = self.markers_of_LG[LG][:]
			for marker in group_to_iter_over:
				self.drop_markers(marker)
			if self.markers_of_LG.has_key(LG):
				if len(self.markers_of_LG[LG]) == 0:
					self.markers_of_LG.pop(LG)
					print ("removed LG:\t{0} from {1}".format (LG, self.name))
		else:
			print ("LG\t{0}\tNOT found in {1}".format(LG, self.name))
			return None
		

# TO RUN
def marker_all_gone(marker, *args):
	cM_test = dict()
	LG_test = dict()
	marker_LG_test = dict()
	marker_cM_test = dict()
	for map in args:
		if map.cM_of_marker.has_key(marker):
			print ("{0}\tin map\t{1}\t at cM {2} (cM_of_marker) ".format(marker, map.name, map.cM_of_marker[marker]))
		if map.LG_of_marker.has_key(marker):
			print ("{0}\tin map\t{1}\t at LG {2} (LG_of_marker) ".format(marker, map.name, map.LG_of_marker[marker]))
		for LG in map.markers_of_LG.keys():
			if marker in map.markers_of_LG[LG]:
				print ("{0}\tin map\t{1}\t at LG {2} (markers_of_LG) ".format(marker, map.name, LG))
		for cM in map.markers_of_cM.keys():
			if marker in map.markers_of_cM[cM]:
				print ("{0}\tin map\t{1}\t at cM {2} (markers_of_cM) ".format(marker, map.name, cM))
	return None
		
def LG_all_gone(LG, *args):

	return None

def compare_maps (*args):
	"""Function to compare genetic maps.  
	Prints out a quick summary, returns None"""
	# test all arguments are maps
	for map in args:
		if not isinstance(map, Map):
			print(map + " is not a Map object")
			return None
			
	map_names 	= list()
	map_types	= list()
	marker_sets = list()
	marker_LGs 	= list()
	marker_cMs 	= list()
	
	for map in args:
		map_names.append(map.name)
		map_types.append(map.map_type)
		marker_sets.append(frozenset(map.LG_of_marker.keys()))
		marker_LGs.append(map.LG_of_marker)
		marker_cMs.append(map.cM_of_marker)
	
	# Print Per-Map Summary
	for xx in range(len(map_names)):
		print("Map Name:\t" + map_names[xx])
		print("Map Type:\t" + map_types[xx])
		print("Linkage Groups:\t" + str(len(set(marker_LGs[xx].values()))))
		print("Markers:\t" + str(len(marker_sets[xx])))
		# Unique Markers
		all_other = tuple(marker_sets[0:xx] + marker_sets[xx+1:])
		unique_set = marker_sets[xx].difference(*all_other)
		print("Map Unique Markers:\t"+ str(len(unique_set)))
		print ("-----------")
	# Overall Summary
	union_marker_set = marker_sets[0].union(*marker_sets)
	intersection_marker_set = marker_sets[0].intersection(*marker_sets)
	print ("All Markers:\t" + str(len(union_marker_set)))
	print ("Intersection Markers:\t" + str(len(intersection_marker_set)))
	return None


def write_intersection (file, map_1 ,*args):
	"""Write the LG and cM for each marker shared by all supplied maps
	"""
	
	out_file = open (file, 'w')
	if not isinstance(map_1, Map):
		print("map_1" + " is not a Map object")
		return None
	for map in args:
		if not isinstance(map, Map):
			print(map + " is not a Map object")
			return None
	
	map_names 	= list()
	map_types	= list()
	marker_sets = list()
	marker_LGs 	= list()
	marker_cMs 	= list()
	
	for map in args:
		marker_sets.append(set(map.LG_of_marker.keys()))
			
	#Write header line
	to_write = map_1.name, map_1.name
	out_file.write("Marker\t{} LG\t{} cM".format (*to_write))
	for map in args:
		to_write = map.name, map.name
		out_file.write("\t{} LG\t{} cM\n".format (*to_write))
	
	# intersection_marker_set will be the set of markers shared by all supplied maps
	intersection_marker_set = marker_sets[0].intersection(set(map_1.LG_of_marker.keys()),*marker_sets)
	
	for marker in intersection_marker_set:
		to_write = str(marker), map_1.LG_of_marker[marker], map_1.cM_of_marker[marker]
		out_file.write("{}\t{}\t{}".format(*to_write))
		for map in args:
			to_write = map.LG_of_marker[marker], map.cM_of_marker[marker]
			out_file.write("\t{}\t{}\n".format(*to_write))
	out_file.close()
	
	
	
def write_union (file, map_1 ,*args):
	"""Write the LG and cM for each marker in each supplied map.
	Matches markers by name.
	"""
		
	if not isinstance(map_1, Map):
		print("map_1" + " is not a Map object")
		return None
	for map in args:
		if not isinstance(map, Map):
			print(map + " is not a Map object")
			return None
	
	out_file = open (file, 'w')
	map_names 	= list()
	map_types	= list()
	marker_sets = list()
	marker_LGs 	= list()
	marker_cMs 	= list()
	
	for map in args:
		marker_sets.append(set(map.LG_of_marker.keys()))
			
	#Write header line
	to_write = map_1.name, map_1.name
	out_file.write("Marker\t{}_LG\t{}_cM".format(*to_write))
	for map in args:
		to_write = map.name, map.name
		out_file.write("\t{}_LG\t{}_cM".format(*to_write))
	out_file.write("\n")
	
	# Will be the set of markers shared by all supplied maps
	union_marker_set = marker_sets[0].union(set(map_1.LG_of_marker.keys()),*marker_sets)
	
	for marker in union_marker_set:
		map_1_to_write = str(marker), map_1.LG_of_marker.get(marker, "NA"), map_1.cM_of_marker.get(marker, "NA")
		
		out_file.write("{0}\t{1}\t{2}".format(*map_1_to_write))
		for map in args:
			map_to_write = map.LG_of_marker.get(marker, "NA"), map.cM_of_marker.get(marker, "NA")
			out_file.write("\t{0}\t{1}".format(*map_to_write) )
		out_file.write("\n")
	out_file.close()
	


def align_maps (map_1, map_2):
	"""Rename the LGs in map_2 to match the names in map_1.
	Two Linkage groups match and trigger renaming if they share two or more markers, 
	if a pair of Lgs share only one marker, and error will be returned
	and no renaming will be done.
	Unmatched LGs in map_2 will be prefixed with 'un_'.
	returns a tuple of dictionaries that descibe the matches made.
	"""
	
	map_2_LGs_matching_map_1_LG = dict() # values will lists 
	map_1_LGs_matching_map_2_LG = dict() # values will lists 
	
	single_marker_matches 	= 0
	multiple_marker_matches	= 0
	
	unmatched_map_1_LGs		= map_1.markers_of_LG.keys() # initially all LGs
	unmatched_map_2_LGs		= map_2.markers_of_LG.keys() # initially all LGs
	do_renaming = 1 # yes
		
	# Try to match each LG in map_1 to each LG map_2:	
	for LG_in_map_1 in map_1.markers_of_LG.keys():
		for LG_in_map_2 in map_2.markers_of_LG.keys():
			# overlap is the set of markers that occur on the two LGs being compared
			overlap = set(map_1.markers_of_LG[LG_in_map_1]).intersection(set(map_2.markers_of_LG[LG_in_map_2]))
			# possible is max number of possible matches, determined as the number of markers on the shorter LG. 
			possible = min(len(map_1.markers_of_LG[LG_in_map_1]), len(map_2.markers_of_LG[LG_in_map_2]))
			# Find all matches
			if len(overlap) > 1: 		# match
				multiple_marker_matches = multiple_marker_matches + 1
				if LG_in_map_1 in unmatched_map_1_LGs:
					unmatched_map_1_LGs.remove(LG_in_map_1)
				if LG_in_map_2 in unmatched_map_2_LGs:
					unmatched_map_2_LGs.remove(LG_in_map_2)
				map_2_LGs_matching_map_1_LG.setdefault(LG_in_map_1, [])
				map_2_LGs_matching_map_1_LG[LG_in_map_1].append(LG_in_map_2)
				map_1_LGs_matching_map_2_LG.setdefault(LG_in_map_2, [])
				map_1_LGs_matching_map_2_LG[LG_in_map_2].append(LG_in_map_1)
				print (LG_in_map_1 + "\t" + LG_in_map_2 + "\tmatched at " + str(len(overlap)) +" / "+ str(possible) + "\tmarkers\t" + str(float(len(overlap))/(possible)))
			elif len(overlap) == 1: 	# single_marker_match
				single_marker_matches = single_marker_matches + 1
				for x in overlap: print ("ERROR\n"+ LG_in_map_1 + "\t" + LG_in_map_2 + " share only marker " + x)
	
	print ('unmatched_map_1_LGs')
	print (unmatched_map_1_LGs)
	print ('unmatched_map_2_LGs')
	print (unmatched_map_2_LGs)				
		
	# Determine if renaming can go forward		
	if single_marker_matches > 0: # 
		do_renaming = 0 # stop renaming, some LG matched at 1 marker
	if len([x for x in map_1_LGs_matching_map_2_LG.values() if len(x)>1]):
		do_renaming = 0 # stop renaming, some map_1 LGs matched multiple in map _2
		print ("*WARNING* LGs in map_1 match multiple LGs in map_2")
	if len([x for x in map_2_LGs_matching_map_1_LG.values() if len(x)>1]):
		do_renaming = 0 # stop renaming, some map_2 LGs matched multiple in map _1
		print ("*WARNING* LGs in map_2 match multiple LGs in map_1")
	if do_renaming == 0	:
		print ("Markers will not be renamed")	
	# Rename markers	
	if do_renaming == 1:
		# *All* LGs in map_2 will be renamed.
		# LGs in map_1 will not be renamed.
		# Matching map_2 LGs will have names changed to be the same as in map_1.  
		# Non-matching map_2 LGs will also be preceeded with "un_".
		
		new_dict = dict()
		# rename each unmatched LG.
		for each_LG in unmatched_map_2_LGs: 
			new_name = 'un_' + each_LG
			new_dict[new_name] = map_2.markers_of_LG[each_LG]
			print ("renamed " + each_LG + " to " + new_name)
			for each_marker in new_dict[new_name]:
				map_2.LG_of_marker[each_marker] = new_name
				

		for k in map_2.markers_of_LG.keys(): # rename matching groups
			if k in map_1_LGs_matching_map_2_LG:
				new_name = map_1_LGs_matching_map_2_LG[k][0]
				new_dict[new_name] = map_2.markers_of_LG[k]
				print (("renamed " + k + " to " + new_name))
		map_2.markers_of_LG = new_dict
		count = 0
		for k, v in map_2.markers_of_LG.items():
			count = count + 1
			for each_marker in v:
				map_2.LG_of_marker[each_marker] = k

	# currently returns the outdated matchings
	return [map_2_LGs_matching_map_1_LG, map_1_LGs_matching_map_2_LG]




def compare_orders (map_1, map_2, match_order = 'yes'): 
	"""Compare the ordering (not spacing) of the markers on two maps.
	LGs are matched by name.
	Markers are matched by name.
	Run align_maps() previous to this to ensure LGs match.
	Similarity in ordering is measured by kendall's tau statisic.
	The map order of markers in map_2 will be inverted 
	to better match order in map_1 if (tau < 0)	so that tau (> 0).
	"""
	tau_of_LG = dict()
	print ("LG" + "\t" + "Markers" + "\t" + "cM" + "\t" + "tau")
	for each_LG in map_1.markers_of_LG.keys():
		intersection_marker_set = set(map_1.markers_of_LG[each_LG]).intersection(set(map_2.markers_of_LG[each_LG]))
		cM_map_1 = []
		cM_map_2 = []

		for each_marker in intersection_marker_set:
			cM_map_1.append(float(map_1.cM_of_marker[each_marker]))
			cM_map_2.append(float(map_2.cM_of_marker[each_marker]))
		tau = stats.stats.kendalltau(cM_map_1, cM_map_2)[0]
		map_1_cM = max([float(x) for x in cM_map_1 ])
		print (each_LG + "\t" + str(len(cM_map_1))+ "\t"  + str(map_1_cM) + "\t"  + str(tau) )
		
		if match_order.lower() == 'yes':
			if tau < 0:
			# reverse order on the LG from map_2 
				find_max_cM_map_2 = []
				for each_marker in map_2.markers_of_LG[each_LG]:
					find_max_cM_map_2.append(float(map_2.cM_of_marker[each_marker]))
				max_cM_map_2 = max ([float(x) for x in find_max_cM_map_2 ])
				# set new cM positions
				for each_marker in map_2.markers_of_LG[each_LG]:
					map_2.cM_of_marker[each_marker] = max_cM_map_2 - float(map_2.cM_of_marker[each_marker])

				# recalculate tau
				cM_map_1 = []
				cM_map_2 = []
				for each_marker in intersection_marker_set:
					cM_map_1.append(float(map_1.cM_of_marker[each_marker]))
					cM_map_2.append(float(map_2.cM_of_marker[each_marker]))
				tau = stats.stats.kendalltau(cM_map_1, cM_map_2)[0]	
				print (each_LG + "\t" + str(len(cM_map_1))+ "\t"  + str(max_cM_map_2) + "\t"  + str(tau) )
		tau_of_LG[each_LG] = tau
	return tau_of_LG
			

#######################
class Marker_data_SQL:
	"""SQL Testing, Assumes an empty database, with four tables"""
	def __init__(self, file_name, active_cursor, pop_name, marker_name_line = 1, LG_line = 2, cM_line = 3, genotype_start_line = 4):
		data_file = open (file_name, 'r')
		cursor = active_cursor
		line_count = 0
		markers_entered = 0
		data_rows = list()
		ind_names = list()
		
		marker_names_of_index = dict()
		LG_of_marker_at_index = dict()
		cM_of_marker_at_index = dict()
		
		for line in data_file:
			line_count = line_count + 1
			if line_count == marker_name_line:
			# Parse marker names
				markers_line = line.strip().split(",")
				markers = markers_line[1:]
			elif line_count == LG_line:
				LGs_line = line.strip().split(",")
				LGs = LGs_line[1:]
			elif line_count == cM_line:
				cMs_line = line.strip().split(",")
				cMs = cMs_line[1:]
				
			elif line_count >= genotype_start_line:
				if markers_entered == 0:
					markers_entered = 1
					#enter markers once
					for index in range(len(markers)):
						marker_names_of_index[index] = markers[index]
						LG_of_marker_at_index.setdefault(LGs[index], None)
						cM_of_marker_at_index.setdefault(cMs[index], None)
						LG_of_marker_at_index.setdefault(index, None)
						cM_of_marker_at_index.setdefault(index, None)
					for index in range(len(markers)):
						marker_sql_values = (marker_names_of_index[index],
							LG_of_marker_at_index[index],
							cM_of_marker_at_index[index]
							)
						cursor.execute("INSERT INTO Markers VALUES(null, ?, ?, ?, null)", marker_sql_values)
					conn.commit()
						
				#continue parse marker data
				marker_data = line.strip().split(",")
				ind_names.append(marker_data[0])
				data_rows.append(marker_data[1:])
				
		data_file.close()
		# Insert Individual data
		for ind_index in range(len(ind_names)):
			ind_name = ind_names[ind_index]
			individual_sql_values = ind_name, pop_name
			cursor.execute("INSERT INTO Individuals VALUES (null, ?, ?)", individual_sql_values)
			cursor.execute("SELECT Ind_ID FROM Individuals WHERE Ind_Name=?", [ind_name])
			Ind_ID = cursor.fetchone()[0]
			for genotype_index in range(len(data_rows[ind_index])):
				marker_name = marker_names_of_index[genotype_index],
				cursor.execute("SELECT Marker_ID FROM Markers WHERE Marker_Name=?", marker_name)
				Marker_ID = cursor.fetchone()[0]
				genotype = data_rows[ind_index][genotype_index]
				genotype_sql_values = Ind_ID, Marker_ID, genotype
				cursor.execute("INSERT INTO Genotypes VALUES (null, ?, ?, ?)", genotype_sql_values)
		conn.commit()
		return None
	
	def write_mst(self, file_name, header = 'yes'):
		mst_file = open(file_name, 'w')
		header = ["population_type DH\n", "population_name **\n", "distance_function kosambi\n",
			"cut_off_p_value .000001\n", "no_map_dist 20.0\n", "no_map_size 2\n", "missing_threshold **\n",
			"estimation_before_clustering yes\n", "detect_bad_data yes\n", "objective_function COUNT\n",
			"number_of_loci " + str(len(self.cols)-1) + "\n", "number_of_individual " + str(len(self.rows)-1) + "\n",
			"\n" ]
		if header.lower() == 'yes':
			mst_file.writelines(header)
		mst_file.writelines(["\t".join(col)+"\t\n" for col in self.cols])
		mst_file.close()
		
	def write_rqtl(self, file_name, inc_cM = 'yes'):
		rqtl_file = open(file_name, 'w')
		# TODO
		rqtl_file.close()
		



class MST_input:
	def __init__(self, file_path):
		self.file = open(file_path, 'r')
		for line in self.file:
			if line.startswith("population_type"):
				self.population_type = line.strip().split()[1]
			elif line.startswith("population_name"):
				self.population_name = line.strip().split()[1]
			elif line.startswith("distance_function"):
				self.distance_function = line.strip().split()[1]
			elif line.startswith("cut_off_p_value"):
				self.cut_off_p_value = line.strip().split()[1]
			elif line.startswith("no_map_dist"):
				self.no_map_dist = line.strip().split()[1]
			elif line.startswith("no_map_size"):
				self.no_map_size = line.strip().split()[1]
			elif line.startswith("missing_threshold"):
				self.missing_threshold = line.strip().split()[1]
			elif line.startswith("estimation_before_clustering"):
				self.estimation_before_clustering = line.strip().split()[1]
			elif line.startswith("detect_bad_data"):
				self.detect_bad_data = line.strip().split()[1]
			elif line.startswith("objective_function"):
				self.objective_function = line.strip().split()[1]
			elif line.startswith("number_of_loci"):
				self.number_of_loci = line.strip().split()[1]
			elif line.startswith("number_of_individual"):
				self.number_of_individual = line.strip().split()[1]
			elif len(line) < 4:
				pass # possible blank line
			elif line.startswith("locus_name"):
				pass # parse ind line here
			else:
				pass #marker/genotype data line
		self.file.close()
		

	

class Marker_data_dev:
	"""Class to represent the raw genetic data at a set of markers and individuals,
	before markers are placed in LG or ordered within them.
	Input data format is rqtl-like with *no* cM line."""
	def __init__(self, file_name, pop_name, marker_name_line = 1, LG_line = 2, cM_line = 3, genotype_start_line = 4):
		data_file = open(file_name, 'r')
		line_count = 0
		self.rows = list()
		self.pop_name = pop_name
		
		for line in data_file:
			line_count = line_count + 1
			if line_count == marker_name_line:
			# Parse marker names
				markers = line.strip().split(",")
				markers[0] = "locus_name"
				self.rows.append(markers)
			elif line_count == LG_line:
				dummy_LGs = line.strip().split(",")
				# dont append
			elif line_count == cM_line:
				cMs_line = line.strip().split(",")
				cMs = cMs_line[1:]
				# dont append
			elif line_count >= genotype_start_line:
				marker_data = line.strip()
				# remap genotypes here
				marker_data = marker_data.replace(",AA", ',a')
				marker_data = marker_data.replace(",AB", ',b')
				marker_data = marker_data.replace(",BB", ',b')
				marker_data = marker_data.split(",")
				self.rows.append(marker_data)
				
		self.cols = zip(*self.rows)
		data_file.close()
	
	def write_mst(self, file_name, inc_header = 'yes', p_value = .00001, pre_est = 'yes', detect = 'yes', obj_func = "COUNT"):
		"writes out the marker data in MSTmap format with header"
		mst_file = open(file_name, 'w')
		header = ["population_type DH\n",
			"population_name " + str(self.pop_name) + "\n", 
			"distance_function kosambi\n",
			"cut_off_p_value .00001\n", 
			"no_map_dist 50.0\n", 
			"no_map_size 0\n", 
			"missing_threshold .5\n",
			"estimation_before_clustering "+ pre_est + "\n",
			"detect_bad_data "+ detect + "\n",
			"objective_function "+ obj_func + "\n",
			"number_of_loci " + str(len(self.cols)-1) + "\n",
			"number_of_individual " + str(len(self.rows)-1) + "\n",
			"\n" ]
		mst_file.writelines(header)
		mst_file.writelines(["\t".join(col)+"\t\n" for col in self.cols])
		mst_file.close()
		
	def write_rqtl(self, file_name, inc_cM = 'yes'):
		"not implemented"
		rqtl_file = open(file_name, 'w')
		rqtl_file.close()
		
	def write_sql(self, active_cursor):
		"not implemented"
		return None




