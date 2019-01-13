#!/usr/bin/python
import sys

def median(newlist):
    if len(newlist) % 2 == 0:
        medianpos = len(newlist)/2  
        return float(newlist[medianpos] + newlist[medianpos-1]) /2
    else:
        medianpos = len(newlist)/2
        return newlist[medianpos]

def N50(numlist):
    """
    Returns the N50 value of the passed list of numbers. 

    Based on the Broad Institute definition:
    https://www.broad.harvard.edu/crd/wiki/index.php/N50
    """
    numlist.sort()
    newlist = []
    for x in numlist :
        newlist += [x]*x

    return median(newlist)

lengths = []
for line in sys.stdin :
    lengths.append(int(line))
print "N50:", N50(lengths)
