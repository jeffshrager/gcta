import sim
import learner
from matrixtest import score_matrix
from random import randint
import matplotlib.pyplot as plt
import json

params = json.load(open('params.json', 'r'))

# generate starter data
sim_settings = sim.Settings('params.json')
simulator = sim.Simulator(sim_settings)
starter_data = simulator.startup_loop()

lrnr_settings = learner.Settings('params.json')
lrnr = learner.Learner(lrnr_settings)
weight_matrix, rms_vector = lrnr.learn(starter_data)

# score_matrix(weight_matrix)

print(weight_matrix)
print(rms_vector)

simulator.learner_results(weight_matrix, rms_vector)

mdks = []

def main_loop(num_months_to_sim):
    simulator.treatment_effectiveness = []
    patients = []
    # a whole bunch of months
    m = 1
    num_dead = 0
    num_cured = 0
    while m <= num_months_to_sim:
        print("Simulating month", m, ".")
        # each month, generate between 30 and 50 patients
        for j in range(randint(30, 50)):
            new_px = sim.Patient(len(patients)+1, m, sim_settings.markers)
            patients.append(new_px)

        undead_and_uncured_patients = []
        for patient in patients:
            if patient.state == 'CURED':
                num_cured += 1
            elif patient.state == 'DEAD':
                num_dead += 1
            else:
                undead_and_uncured_patients.append(patient)

        patients = undead_and_uncured_patients
        # patients = [patient for patient in patients if (patient.state != 'CURED' and patient.state != 'DEAD')]

        monthly_log = []

        for patient in patients:
            # treat patient, and get updates on patient responses
            updates = simulator.hypothesis_treatment(patient)
            monthly_log = monthly_log + updates

        # dks = [d['delta_k'] for d in monthly_log]
        mdks.append(sum(simulator.treatment_effectiveness)/len(simulator.treatment_effectiveness))

        def mean(l):
            return sum(l)/len(l)

        # mdks.append(mean(dks))

        # Update the beliefs of the sim
        simulator.update_beliefs(patients)
        m += 1
    print("num dead:", num_dead)
    print("num cured:", num_cured)

main_loop(params['num_months_run'])

# print(simulator.rms_vector)

plt.plot([i for i in range(len(mdks))], mdks)
plt.show()