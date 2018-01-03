#!/usr/bin/python
import sys, os, argparse

parser = argparse.ArgumentParser()
parser.add_argument("src")
parser.add_argument("dest")
args = parser.parse_args()

if __name__ == '__main__':
	src_path = args.src
	dest_path = args.dest

	with open(src_path, 'r') as fread:
		raw = [line for line in fread if (not(line.endswith('-\n')))]

	with open(dest_path, 'w') as fwrite:
		[fwrite.write(line) for line in raw]