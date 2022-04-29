import numpy

w_mat = [[-9.18347168, -9.27240086, -8.92519283, -8.20099926, -3.58878183
, -18.02451515, -12.67471313, -11.94957638, -10.71409225, -19.15688515
, -11.10240555]
, [-12.1747551, -14.01770496, -14.46427059, -15.10067368, -21.46432877
, -6.67121363, -10.88628387, -11.98709869, -13.90428638, -23.67790222
, -13.49435616]
, [-19.42021561, -17.59747314, -18.76243782, -19.97263908, -25.35348511
, -24.11816216, -19.88859749, -19.16145134, -17.79005623, -6.8203001
, -18.91627121]
, [-16.25148582, -15.73608685, -13.80775261, -18.00488281, -22.50829315
, -6.43483829, -17.85696411, -15.90418816, -14.36449718, -24.9148941
, -14.45423985]
, [-21.25315666, -22.06438637, -22.54598618, -17.07777405, -31.10130501
, -29.74467087, -17.86845207, -22.21762276, -22.33332253, -28.81605148
, -22.04460716]
, [-15.21536446, -14.83707428, -15.40761566, -16.15569878, -3.95406938
, -17.72766876, -15.38029194, -12.91432858, -14.67711926, -5.64475489
, -13.57172108]]


backwards = {'BEVACIZUMAB': 0, 'TEMOZOLOMIDE': 1, 'CABAZITAXEL': 2, 'TCAR': 3, 'DISATINIB': 4, 'NIVOLUMAB': 5, 
                 'DOXORUBICIN': 6, 'DURVALUMAB': 7, 'PEMBROLIZUMAB': 8, 'VARLILUMAB': 9, 'SORAFENIB': 10}

drug_dict = dict((v,k) for k,v in backwards.items())

mut_dict = {0: 'MGMT-Promoter', 1: 'HLA-A2', 2: 'HLA-A1', 3: 'IDH2', 4: 'IDH1', 5: 'EGFR'}

x1 = [0, 0, 0, 0, 0, 1]
x2 = [0, 0, 0, 0, 1, 0]
x3 = [0, 0, 0, 1, 0, 0]
x4 = [0, 0, 1, 0, 0, 0]
x5 = [0, 1, 0, 0, 0, 0]
x6 = [1, 0, 0, 0, 0, 0]

muts = [x1, x2, x3, x4, x5, x6]

preds = {}

for i in range(len(muts)):
    mut = muts[i]
    scores = list(numpy.matmul(mut, w_mat))
    print('\n\n')
    # print(mut)
    print(mut_dict[i])
    print(scores)
    scores_copy = scores[:]
    scores_copy.sort()

    ranked_drugs = []
    for j in range(1, len(scores_copy)+1):
        ranked_drugs.append(scores.index(scores_copy[-j]))
    
    for k in range(len(ranked_drugs)):
        ranked_drugs[k] = drug_dict[ranked_drugs[k]]    

    print(ranked_drugs)
    print('stddev', str(numpy.std(scores)))
    print('mean', str(numpy.mean(scores)))

    preds[mut_dict[i]] = ranked_drugs

    print('====================================')