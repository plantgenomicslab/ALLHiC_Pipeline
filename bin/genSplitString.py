#import argparse
#
#arg_parser = argparse.ArgumentParser()
#arg_parser.add_argument("-s", "--split", default = 1,
#						help = "Integer value to interpret and generate corresponding split string.")
#args = arg_parser.parse_args()
#
#x = args.split
def partSuffixes(x):
		splitSuffix = list("000")
		split_arr = []

		###
		# @idx Range from zero to split value
		###
		for idx in range(1,int(x)+1):
			###
			# @i Range from zero to size of integer idx
			# @j Range from zero to two for position in part mold
			###
			for i,j in zip(range(len(str(idx))-1,-1,-1), range(2,-1,-1)):
				splitSuffix[j] = str(idx)[i]
			split_arr.append("".join(splitSuffix))
		return split_arr
