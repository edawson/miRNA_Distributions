#!/usr/bin/python
import csv
import sys

def transpose(feed, delim):
	delim = "\t"
	reader = csv.reader(feed, delimiter=delim)
	for col in zip(*reader):
		print delim.join(col)


if __name__ == "__main__":
	delim = "\t"
	transpose(sys.stdin, delim)
