import sys
f = open('sim-mrhlog.xls', 'r')
l = f.readlines()
dks = [int(li.split()[-1]) for li in l]
# dks = dks[:500000]
averages = []

def average(l):
    return sum(l)/len(l)

group_size = int(sys.argv[-1])
for i in range(int(len(dks)/group_size)):
    averages.append(average(dks[group_size*i:group_size*(i+1)]))

f.close()
f = open('average-delta-ks.xls', 'w')
for v in averages:
    f.write('"'+str(v)+'"\n')
f.close()

import matplotlib.pyplot as plt
plt.plot(list(range(len(averages))), averages)
plt.show()
