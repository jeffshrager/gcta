import matplotlib.pyplot as plt
import sys
import os.path
import json

print("Usage is: python mavg.py window_width subsample_spacing")

window_width = int(sys.argv[-2])
subsample_spacing = int(sys.argv[-1])
global dks

def read_data(path):
    global dks
    print("Reading from " + path)
    f = open(path, 'r')
    l = f.readlines()
    dks = [int(li.split()[-1]) for li in l]
    f.close()

# Get the data into dks
if (os.path.isfile('sim-mrhlog.xls') == True):
    print("Reading local data")
    read_data('sim-mrhlog.xls')
elif (os.path.isfile('parameters.json') == True):
    print("Reading params to get latest data file")
    run_label = json.load(open('parameters.json'))['run-label']
    print("Reading data from "+run_label)
    read_data('experiments/'+run_label+'/sim-mrhlog.xls')
else:
    sys.exit("I'm confused, there's no local data, nor am I in the RL directory where the params are?!?!")

# Sum up the first set
sum = 0

for i in range(window_width):
    sum = sum + dks[i]

# That becomes the first average
averages = [sum]

# Now move along changing the sum as you go, subtracting the head, and adding on the next

for i in range(window_width+1,len(dks)):
    sum = sum - dks[i-window_width] + dks[i]
    averages.append(sum)

# Now extract based on subsample spacing and turn into an average

flwinwid = float(window_width)
averages = [averages[i]/flwinwid for i in range(0,len(averages),subsample_spacing)]

plt.plot(list(range(len(averages))), averages)
plt.show()
