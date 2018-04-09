#!/usr/bin/env python3

# example: 
import sys

if len(sys.argv) != 2:
    print ("usage: parse_significance_logfile.py <name of log file>")
    sys.exit()

filename = str(sys.argv[1])
inputfile = open(filename)

for line in inputfile:
    sp = line.split(" = ")
    if sp[0] == "n_sigma":
        nsigma = float(sp[1])
        print (round(nsigma,1))
    
