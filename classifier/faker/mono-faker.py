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
    for l in range(0,nterms):
        # The term value is either the "correct" number with p(corr), or a random number
        for m in range(0,nsamples):
            if (random.random()<=corr):
                lines.append([l,l,1]) # if we are less than p(corr) we append add a line containing that term and itself as a float
            else:
                lines.append([l,int(random.random()*nterms),0]) #otherwise its a random number b/w 0-nterms (i.e. noise)
    lines=random.sample(lines, len(lines)) # randomizing the entire list we just generated so it's not in any particular order
    for l in range(len(lines)):
        o.write(alphabet[lines[l][0]] + '\t' + str(lines[l][1]) + '\t' + str(lines[l][2]) + '\n')
    o.close()
else:
    print("Usage: ./faker.py filename nsamples nterms(up to 10) correlation(0.0->1.0)")