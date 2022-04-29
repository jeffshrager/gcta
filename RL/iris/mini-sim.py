import random, sys

class Patient:

    def __init__(self, pxnum, markers):
        self.pxnum = pxnum
        self.k_score = random.randint(70, 90)
        self.markers = markers
        self.tumor = self.gen_tumor()

    def gen_tumor(self):
        tumor = []
        p_mut = 1/len(self.markers)
        for m in self.markers:
            if random.random() <= p_mut:
                tumor.append(m)
        if len(tumor) == 0:
            #pick one mutation and return it!
            mut = random.randint(0, len(self.markers)-1)
            tumor = [self.markers[mut]]
        return tumor

    def update_health(self, delta_k):
        self.k_score += delta_k


class Simulator:

    def __init__(self, num_markers, num_drugs, num_drugs_to_give=2):
        self.num_patients = 500
        self.num_months = 30
        self.num_drugs_to_give = num_drugs_to_give
        self.targets_per_drug = 1
        self.p_mutational_escape = 0.25
        self.all_drugs = ["BEVACIZUMAB", "TEMOZOLOMIDE", "CABAZITAXEL", "TCAR", "DISATINIB", "NIVOLUMAB", "DOXORUBICIN", "DURVALUMAB", "PEMBROLIZUMAB", "VARLILUMAB", "SORAFENIB"]
        self.all_markers = ["MGMT-PROMOTER", "HLA-A2", "HLA-A1", "IDH2", "IDH1", "EGFR"]
        self.drugs = self.all_drugs[:num_drugs] # use the first num_drugs drugs
        self.markers = self.all_markers[:num_markers] # also use the first num_markers markers
        self.targets = self.make_targets()
        self.log_file = open('minisim-mrhlog.xls', 'w')
        self.info_dump = open('minisim-information.json', 'w')

    def cleanup(self):
        info = {'markers': self.markers, 'drugs': self.drugs, 'targets': self.targets, 
                'drugs-in-cocktail': self.num_drugs_to_give, 'p-mutational-escape': self.p_mutational_escape,
                'num-patients': self.num_patients, 'num-months': self.num_months}
        self.info_dump.write('{\n')
        bits = list(info.keys())
        for i in range(len(bits)):
            bit = bits[i]
            val = info[bit]
            to_write = '\t"' + bit + '": ' + str(info[bit]).replace('\'', '"')
            if i + 1 < len(bits):
                to_write += ',\n'
            else:
                to_write += '\n'
            self.info_dump.write(to_write)
        self.info_dump.write('}\n')
        self.info_dump.close()
        self.log_file.close()

    def make_targets(self):

        # MODE: Random assignment of targets
        # targets = {}
        # for drug in self.drugs:
        #     markers_indxs = random.sample(range(len(self.markers)), self.targets_per_drug)
        #     markers = [self.markers[i] for i in markers_indxs]
        #     targets[drug] = markers

        # MODE: Non-convoluted assignment of one target to drug
        targets = {self.drugs[i]: [self.markers[i]] for i in range(len(self.drugs))}

        return targets

    def find_best_drugs(self, tumor):
        #iterate through all drugs, rank them in order of effectiveness
        drug_scores = {} #key: drug, #value: num_overlaps
        for drug in self.drugs:
            drug_targets = self.targets[drug]
            score = 0
            for target in drug_targets:
                if target in tumor:
                    score += 1
            drug_scores[drug] = score

        #then, choose the top two and return
        drugs = []
        best_score = 2
        while len(drugs) < self.num_drugs_to_give:
            # print(drugs)
            if best_score < 0:
                #return a random sampling of drugs if we (somehow...) don't have top two
                rdrugs_indxs = random.sample(range(len(self.drugs)), self.num_drugs_to_give)
                print('something went wrong with ranking drugs, returning two random ones...')
                return [self.drugs[i] for i in rdrugs_indxs]
            else:
                drugs_with_best_score = [drug for drug in drug_scores if drug_scores[drug] == best_score]
                # print('drugs with best score')
                # print(drugs_with_best_score)
                drugs_needed = self.num_drugs_to_give - len(drugs)
                if len(drugs_with_best_score) <= drugs_needed:
                    drugs = drugs + drugs_with_best_score[:]
                else:
                    #sample drugs_needed of them
                    # print('sampling')
                    drugs_indxs = random.sample(range(len(drugs_with_best_score)), drugs_needed)
                    # print('drugs indxs')
                    # print(drugs_indxs)
                    for ind in drugs_indxs:
                        drugs.append(drugs_with_best_score[ind])
                best_score -= 1
        return drugs

    def calculate_delta_k(self, drugs, tumor):
        # two formulas
        # one is we subtract five if there are mutations not covered by the
        # drugs, but add five if all mutations are covered
        # the other is we subtract 5 points for every mutation the patient has
        # and then add 10 points for every drug that covers at least one mut.
        # Here, using second method
        delta_k = -5*len(tumor)
        for drug in drugs:
            good_drug = False
            targeted_muts = self.targets[drug]
            for mut in targeted_muts:
                if mut in tumor:
                    good_drug = True
            if good_drug:
                delta_k += 10

        return delta_k

    def rep_tumor(self, tumor):
        #reverse the list of markers, sum 2^indx in the list of markers
        reversed_markers = self.markers[:]
        reversed_markers.reverse()
        dec_rep = 0
        for marker in tumor:
            dec_rep += 2**reversed_markers.index(marker)
        return dec_rep

    def main_loop(self):
        for px in range(1, self.num_patients+1):
            print("Running patient number", px)
            patient = Patient(px, self.markers)
            cured = False
            dead = False
            for m in range(1, self.num_months+1):
                # each month, we look at the patient's tumor and find the best drug
                # record the patient's tumor, k_score, and drug given
                # then, give the drug and observe delta_k. Record, and then note
                # new k_score. Then, see if the tumor mutates. Finally, record this
                # information
                tumor = patient.tumor # list of markers the patient's tumor has
                k_score_prev = patient.k_score # old k_score (before drugs)

                drugs_to_give = self.find_best_drugs(tumor) # best drugs to give
                delta_k = self.calculate_delta_k(drugs_to_give, tumor) # change in health from drugs
                patient.update_health(delta_k) # update patient health

                k_score_new = patient.k_score # new k_score after treatment
                if k_score_new >= 100:
                    k_score_new = 100
                    cured = True
                elif k_score_new <= 0:
                    k_score_new = 0
                    dead = True

                # check to see if tumor should mutate, do so if it should
                if random.random() <= self.p_mutational_escape:
                    patient.tumor = patient.gen_tumor()

                # after the dust settles, record another line(s) in the log
                # order of values:
                # pxnum, monthnum, tumor (binary rep --> decimal rep), drug_combo_num, 
                # drug_name, drug_num (obsolete), k_new, k_prev, delta_K
                tumor_dec = self.rep_tumor(tumor)
                for di in range(len(drugs_to_give)):
                    vals_to_record = [px, m, tumor_dec, di+1, drugs_to_give[di], di, k_score_new, k_score_prev, delta_k]
                    vals_to_record = [str(val) for val in vals_to_record]
                    to_record = '\t'.join(vals_to_record)
                    self.log_file.write(to_record + '\n')
                if cured:
                    print('Patient number', px, 'was cured!')
                    break
                elif dead:
                    print('Patient number', px, 'died :(')
                    break
            if not dead and not cured:
                print('Patient number', px, 'still sick at end of sim.')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Make sure to run as: $python mini-sim.py num_markers num_drugs")
    else:
        print("Constructing simulator")
        sim = Simulator(int(sys.argv[-2]), int(sys.argv[-1]))
        print("Running simulator...")
        sim.main_loop()
        print("Simulator finished! Cleaning up and exiting.")
        sim.cleanup()




