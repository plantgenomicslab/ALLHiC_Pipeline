import argparse

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-s", "--split", default = 1,
						help = "Integer value to interpret and generate corresponding split string.")
args = arg_parser.parse_args()

x = args.split
splitSuffix = list("000")
for i,j in zip(range(len(x)-1,-1,-1), range(2,-1,-1)):
	splitSuffix[j] = x[i]
print(".part_" + "".join(splitSuffix))
