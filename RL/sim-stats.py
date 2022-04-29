import json
import csv
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt

params = json.load(open('parameters.json', 'r'))

markers = params['simulator']['markers']

drugs_dict = params['simulator']['drugs']
drugs = [0 for i in range(len(drugs_dict))]
for drug in drugs_dict:
    ind = drugs_dict[drug]['index']
    drugs[ind] = drug

# drugs = list(params['simulator']['drugs'].keys())
drugs = [d.upper() for d in drugs]

stats_dict = {marker: {drug: {'dks': [], 'learner': 0, 'n': 0, 'avg': 0, 'std': 0, 'stderr': 0} for drug in drugs} for marker in markers}

f = open('stats-sim-mrhlog.xls', 'r')
log_lines = f.readlines()
log_lines = [line.split() for line in log_lines]
f.close()

for line in log_lines:
    tumor = int(line[2])
    delta_k = int(line[-1])
    drug = line[4]

    # calculate the mutations from the ID
    format_string = "{0:0" + str(len(markers)) + "b}"
    tumor_bin = list(format_string.format(tumor))
    tumor_bin = [int(i) for i in tumor_bin]

    mutations = [markers[i] for i in range(len(markers)) if tumor_bin[i]]

    for mut in mutations:
        stats_dict[mut][drug]['dks'].append(delta_k)

for m in markers:
    for d in drugs:
        stats_dict[m][d]['n'] = len(stats_dict[m][d]['dks'])
        stats_dict[m][d]['avg'] = np.mean(stats_dict[m][d]['dks'])
        stats_dict[m][d]['std'] = np.std(stats_dict[m][d]['dks'])
        stats_dict[m][d]['stderr'] = stats_dict[m][d]['std']/sqrt(stats_dict[m][d]['n'])

cr = open('stats-classifier-results.json', 'r')
cr_lines = cr.readlines()
best_matrix = json.loads(cr_lines[0])['best-matrix']

for i in range(len(best_matrix)): #6 markers
    for j in range(len(best_matrix[i])): #11 drugs
        stats_dict[markers[i]][drugs[j]]['learner'] = best_matrix[i][j]

def check_significant(v1, e1, v2, e2):
    return (v1>v2 and v1-e1>v2+e2) or (v1<v2 and v1+e1<v2-e2)

for m in markers:
    lrnr_stderr_drug_tuples = [(drug, stats_dict[m][drug]['learner'], stats_dict[m][drug]['stderr']) for drug in drugs]
    lrnr_stderr_drug_tuples.sort(key=lambda x: x[1])
    lrnr_stderr_drug_tuples.reverse()

    best_drug, best_lrnr, best_stderr = lrnr_stderr_drug_tuples[0]
    good_drugs = [best_drug]
    boundary_reached = False
    for i in range(1, len(lrnr_stderr_drug_tuples)):
        if boundary_reached:
            break
        drug, lrnr, stderr = lrnr_stderr_drug_tuples[i]
        if check_significant(best_lrnr, best_stderr, lrnr, stderr):
            boundary_reached = True
        else:
            good_drugs.append(drug)

    stats_dict[m]['significance-group'] = good_drugs

out = open('stats.csv', 'w')
writer = csv.writer(out, delimiter=",")
writer.writerow(['MUTATION', 'DRUG', 'LEARNER', 'AVG', 'N', 'STD', 'STDERR'])
for m in markers:
    for d in drugs:
        writer.writerow([m, d, stats_dict[m][d]['learner'], stats_dict[m][d]['avg'], stats_dict[m][d]['n'], stats_dict[m][d]['std'], stats_dict[m][d]['stderr']])
out.close()

for m in markers:
    print(m, ':', ', '.join(stats_dict[m]['significance-group']))

plt.figure()

learner_vals = [[stats_dict[m][d]['learner'] for d in drugs] for m in markers]
err_vals = [[stats_dict[m][d]['stderr'] for d in drugs] for m in markers]

bs = [[(learner_vals[j][i], err_vals[j][i]) for i in range(len(learner_vals[j]))] for j in range(len(learner_vals))]

for b in bs:
    b.sort(key=lambda x: x[0])
    b.reverse()

learner_vals = [[a[0] for a in ar] for ar in bs]
err_vals = [[a[1] for a in ar] for ar in bs]

lines = [plt.errorbar([j for j in range(11)], learner_vals[i], yerr=err_vals[i], label=markers[i]) for i in range(6)]
plt.legend(lines, markers)

plt.show()




