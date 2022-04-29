import csv

f = open('3706369930-mrhlog.xls')
lines = f.readlines()
parsed = [line.split('\t') for line in lines]
f.close()

to_output = []

i = 0
while i < len(parsed):
    if parsed[i][3] == parsed[i+1][3]:
        i += 1
    else:
        to_output.append(parsed[i])
        to_output.append(parsed[i+1])
        i += 2

o = open('370639930-mrhlog-binary-only.xls', 'w')
writer = csv.writer(o, delimiter='\t')
for line in to_output:
    line[-1] = int(line[-1])
    writer.writerow(line)
o.close()