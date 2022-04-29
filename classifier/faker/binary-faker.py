#!/usr/bin/python

import sys
import random

alphabet=['a0','b1','c2','d3','e4','f5','g6','h7','i8','j9'] # Can only have up to ten terms -- boo hoo

if len(sys.argv) == 5:
    o = open(sys.argv[1], 'w')
    nsamples = int(sys.argv[2])
    nterms = int(sys.argv[3])
    corr = float(sys.argv[4])
    lines = []
    # first we generate a list of single alphabet characters, which will be randomly paired with another char later
    for l in range(0,nterms):
        # The term value is either the "correct" number with p(corr), or a random number
        for m in range(0,nsamples):
            if (random.random()<=corr): # if we are less than p(corr) we set the binary cocktail as normal
                lines.append([l, 1, float(l)])
            else: # otherwise its -1, which will become a random number b/w 0-2*nterms (i.e. noise)
                lines.append([l, 1, -1]) #int(random.random()*nterms)
    lines=random.sample(lines, len(lines)) # randomizing the entire list we just generated so it's not in any particular order

    # now we construct the pairs of chars and their combined "effect"
    for l in range(0,len(lines), 2): #changed range to 0-len(lines) -- used to be 1
        # generate random pair line
        index = int(random.random()*nterms)
        rand_pair = [index, 2, float(index)]
        lines.insert(l+1, rand_pair)
        # if we are proceeding as normal (i.e. not noise), we sum their indices and set that as the "effect"
        if lines[l][2] != -1:
            lines[l][2] += lines[l+1][2]
            lines[l+1][2] = lines[l][2]
        # otherwise we set it to a random "effect" b/w 0-2*nterms
        else:
            lines[l][2] = int(random.random()*2*nterms)
            lines[l+1][2] = lines[l][2]
        #now we write the data
        o.write(alphabet[lines[l][0]] + '\t' + str(lines[l][1]) + '\t' + str(lines[l][2]) + "\n")
        o.write(alphabet[lines[l+1][0]] + '\t' + str(lines[l+1][1]) + '\t' + str(lines[l+1][2]) + "\n")
    o.close()
else:
    print("Usage: ./binary-faker.py filename nsamples nterms(up to 10) correlation(0.0->1.0)")