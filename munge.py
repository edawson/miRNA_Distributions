__author__ = 'eric'

import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "-infile", dest="infile", type=str, required=True)
    # parser.add_argument("-c", "-cores", dest="cores", type=int, default=1, required=False)
    # parser.add_argument("-s", "-seqs", dest="seqs", type=str, required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args= parse_args()
    infile = args.infile
    outfile = sys.stdout

    with open(infile, "r") as fi:
        # with open(outfile, "w") as ofi:
        header = fi.readline()
        hlist = header.split("\t")
        hlist = [x.strip()[0:16] if "TCGA" in x else x for x in hlist]
        sys.stdout.write("\t".join(hlist) + "\n")
        for line in fi:
            if not line.startswith("hsa"):
                continue
            else:
                splits = line.split("\t")
                splits[0] = splits[0].split("|")[0]
                sys.stdout.write("\t".join(splits))
