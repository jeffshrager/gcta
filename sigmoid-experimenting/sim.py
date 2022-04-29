import random, sys, math, json
import numpy as np

random.seed()

class Patient:
    def __init__(self, pxnum, month, markers):
        self.pxnum = pxnum
        self.month = month
        self.k_score = random.randint(70, 90)
        self.markers = markers
        self.states = ['CURED', 'DEAD', 'BETTER', 'WORSE', 'SAME']
        self.state = 'WORSE'
        self.treatment = []
        self.gen_tumor()
        self.previous_tumor = self.tumor[:]

    def gen_tumor(self):
        tumor = [int(random.random() < 0.1667) for i in range(len(self.markers))]
        if tumor == [0 for i in range(len(self.markers))]:
            tumor[random.randint(0, len(self.markers)-1)] = 1
        self.tumor = tumor

    def update_health(self, delta_k):
        if self.state != 'DEAD':
            self.k_score += delta_k
            if delta_k == 0:
                self.state = 'SAME'
            elif delta_k < 0:
                self.state = 'WORSE'
            elif delta_k > 0:
                self.state = 'BETTER'

            if self.k_score <= 0:
                self.k_score = 0
                self.state = 'DEAD'
            elif self.k_score >= 100:
                self.k_score = 100
                self.state = 'CURED'

class Simulator:
    def __init__(self, settings):
        self.settings = settings
        self.num_initial_patients = settings.params['num_initial_patients']
        self.num_months_initial = settings.params['num_months_initial']
        self.num_drugs_to_give = settings.params['num_drugs_to_give']
        self.num_drugs_to_give_initial = settings.params['num_drugs_to_give_initial']
        self.p_mutational_escape = settings.params['p_mutational_escape']
        self.belief_update_value = settings.params['belief_update_value']
        self.drugs = settings.drugs
        self.drugs_lookup_list = settings.drugs_lookup_list
        self.markers = settings.markers
        self.markers_lookup_list = settings.markers_lookup_list
        self.log = []
        self.scores = {"BETTER": 0, "SAME": 0, "WORSE": 0} # no longer used....
        self.treatment_effectiveness = []

    def learner_results(self, weight_matrix, rms_vector):
        """
        Called by loop.py to pass in the results from the learner and transform them appropriately
        """
        self.weight_matrix = weight_matrix

        # this block of code sort of centers the RMS values around 1 and makes it so a larger value means more confidence
        rms_vector = [1/r for r in rms_vector]
        c = 1/(sum(rms_vector)/len(rms_vector))
        rms_vector = [c*i for i in rms_vector] 
        self.rms_vector = dict(zip(self.drugs_lookup_list, rms_vector))

    def update_beliefs(self, patients):
        """
        Reweight matrix - care about recent stuff the most.
        Here be dragons.
        """
        for patient in patients:
            # for each drug, look at the targets of the drug and look at the mutations of the tumor
            # if there's overlap, add to values accordingly
            for drug in patient.treatment:
                drug_index = self.drugs[drug]['index']
                marker_indices = [i for i in range(len(patient.previous_tumor)) if patient.previous_tumor[i]]
                # marker_index = patient.previous_tumor.index(1)
                if patient.state == 'BETTER':
                    # self.rms_vector[patient.treatment[0]] += 0.01 # this is sort of a hacky way to track confidence of individual drugs...
                    for marker_ind in marker_indices:
                        self.weight_matrix[marker_ind][drug_index] += self.belief_update_value
                elif patient.state == 'WORSE':
                    for marker_ind in marker_indices:
                        self.weight_matrix[marker_ind][drug_index] -= self.belief_update_value

        print("Treatment effectiveness (mean):", sum(self.treatment_effectiveness)/len(self.treatment_effectiveness))

    def build_treatment_model(self, tumor):
        """
        Returns a list of probabilities for treatments
        """
        # multiply weight_matrix and tumor
        preds = np.matmul(tumor, self.weight_matrix)

        if self.settings.params['use_sigmoid']:
            modulated_preds = [p*self.settings.sigmoid(p, self.settings.params['h_global']) for p in preds]
        else:
            modulated_preds = preds[:]

        modulated_preds = [p-min(modulated_preds) for p in modulated_preds]
        
        normalized_modulated_preds = [p/sum(modulated_preds) for p in modulated_preds]

        return normalized_modulated_preds

    def calculate_delta_k(self, patient):
        delta_k = -5*len([m for m in patient.tumor if m==1])
        tumor = [self.markers_lookup_list[i] for i in range(len(patient.tumor)) if patient.tumor[i]]
        targets = set()
        for drug in patient.treatment:
            for marker in self.drugs[drug]['markers']:
                targets.add(marker)
            delta_k -= self.drugs[drug]['adversity']
        successful_targets = [marker for marker in tumor if marker in targets]
        delta_k += 10*len(successful_targets)
        self.treatment_effectiveness.append(len(successful_targets)/len(tumor))
        return delta_k

    def hypothesis_treatment(self, patient):
        treatment_probs = self.build_treatment_model(patient.tumor)

        if self.settings.params['always_do_best_thing']:
            drugs_and_probs = [(self.drugs_lookup_list[i], treatment_probs[i]) for i in range(len(treatment_probs))]
            drugs_and_probs.sort(key=lambda x: x[1])
            drugs_and_probs.reverse()
            treatment = [drugs_and_probs[0][0]]
        else:
            # sample an action
            treatment = [self.drugs_lookup_list[list(np.random.multinomial(1, treatment_probs)).index(1)] for i in range(self.num_drugs_to_give)]

        patient.treatment = treatment
        updates = self.update_patient(patient)
        return updates

    def startup_choose_treatment(self, patient):
        """
        Chooses random drugs and assigns them to the tumor
        """
        drug_inds = random.sample(range(len(self.drugs)), self.num_drugs_to_give_initial)
        treatment = [self.drugs_lookup_list[i] for i in drug_inds]
        return treatment

    def update_patient(self, patient):
        k_score_prev = patient.k_score
        delta_k = self.calculate_delta_k(patient)
        patient.update_health(delta_k)
        k_score_new = patient.k_score

        if random.random() <= self.p_mutational_escape:
            patient.previous_tumor = patient.tumor[:]
            patient.gen_tumor()

        new_updates = []

        for di in range(len(patient.treatment)):
            update_dict = {'px': patient.pxnum, 'month': patient.month, 'tumor_dec': self.settings.rep_tumor(patient.tumor),
                            'drug_combo_num': di+1, 'drug_name': patient.treatment[di], 'px_state': patient.state,
                            'k_score_prev': k_score_prev, 'k_score_new': k_score_new, 'delta_k': delta_k}
            self.log.append(update_dict) 
            new_updates.append(update_dict)
            self.settings.log_file.write(str(update_dict))

        return new_updates

    def startup_loop(self):
        # go through self.num_initial_patients patients
        for px in range(1, self.num_initial_patients+1):
            print("Running patient number", px)

            patient = Patient(px, 0, self.settings.markers)

            for m in range(1, self.num_months_initial+1):
                # each month, we look at the patient's tumor and if they've
                # gotten worse, we assign a new treatment. Otherwise use the
                # old one. Then check if the tumor mutates, and record all info
                state = patient.state
                if state == 'WORSE':
                    patient.treatment = self.startup_choose_treatment(patient)
                elif state == 'DEAD':
                    print("Patient number {} died at month {}".format(px, m))
                    break
                elif state == 'CURED':
                    print("Patient number {} was cured at month {} :)".format(px, m))
                    break

                self.update_patient(patient)

        return self.log

class Settings:
    def __init__(self, path_to_params):
        self.params = json.load(open(path_to_params, 'r'))
        self.markers = self.params['markers']
        if self.params['num_additional_markers'] > 0:
            self.gen_markers()
        self.markers_lookup_list = self.gen_lookup_list(self.markers)
        self.drugs = self.params['drugs']
        self.gen_drugs()
        self.drugs_lookup_list = self.gen_lookup_list(self.drugs)
        self.log_file = open(self.params['sim_output'], 'w')
        self.sigmoid_shift = 0

    def gen_drugs(self):
        counter = 0
        for drug in self.drugs:
            self.drugs[drug]['markers'] = self.markers_lookup_list[counter:counter+self.params['num_targets_per_drug']]
            counter += self.params['target_shift_per_drug']

    def gen_markers(self):
        labels = []
        for i in range(self.params['num_additional_markers']):
            labels.append('X' + str(i))
        for i in range(len(labels)):
            self.markers[labels[i]] = {"index": i+6, "exponent": 0}
        for marker in self.markers:
            self.markers[marker]['exponent'] = (len(self.markers)-1)-self.markers[marker]['index']

    def gen_lookup_list(self, transform):
        lookup = [None for i in range(len(transform))]
        for item in transform:
            lookup[transform[item]['index']] = item
        return lookup

    def rep_tumor(self, markers):
        dec_rep = 0
        for i in range(len(markers)):
            if markers[i]:
                dec_rep += 2**self.markers[self.markers_lookup_list[i]]['exponent']
        return dec_rep

    def sigmoid(self, x, h):
        # cliff is huge at 10 and good at 1 and a line at 0.1
        return 1/(1+math.exp(-h*x))

    def healing(self, k):
        # this is totally bogus, don't use this
        if 0 <= k < 30:
            return 15
        elif 30 < k < 70:
            return 10
        else:
            return 5


