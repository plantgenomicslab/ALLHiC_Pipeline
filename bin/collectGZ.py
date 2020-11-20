import os

def getGZPairs(gz_list):
	pairs = []
	for gz in gz_list:
		if ".fq.gz" in gz:
			suffix = ".fq.gz"
		elif ".fastq.gz" in gz:
			suffix = ".fastq.gz"
		else:
			print("Error: Gzip file '" + str(gz) + "' does not have '.fastq.gz' or '.fq.gz'. Not used.")
		prefix = gz.replace(suffix,"")
		strand = prefix[-1]
		prefix = prefix[0:-1]
		if strand == "1":
			match_strand = "2"
		elif strand == "2":
			match_strand = "1"
		else:
			print("Error: Gzip file '" + gz + "' does not have a strand value [1,2]. Not used.")
			break
		if prefix + match_strand + suffix in gz_list:
			pairs.append(prefix)
			gz_list.remove(gz)
			gz_list.remove(prefix + match_strand + suffix)
	return pairs

def collectGZ(curr_dir):
	collected_list = []
	for sub in os.listdir(curr_dir):
		if os.path.isdir(curr_dir + sub):
			findings = collectGZ(curr_dir + sub + "/")
			for find in findings:
				collected_list.append(find)
		if sub.endswith(".gz"):
			collected_list.append(curr_dir + sub)
	return collected_list

