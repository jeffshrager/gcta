# How to run
# Make sure that you have the parameters stored in "parameters.json" in the
# same dir as this script. The classifier should also be in this dir, and the
# sim should be in a dir called "sim" in this dir.
# Then, run:
# $python master-script.py


#default python libs
import re
import os
import json
import subprocess
import shutil

#classifier used with this experiment...
from classifier import main as clmain
# from binaryclassifier import main as clmain
# !!!!!!!!!!!!!!!! if you delete your compile-simpx.lisp or compile-simpx.exe file, you must re-compile and re-save it: !!!!!!!!!!!!!!!!
# (load (compile-file "compile-simpx.lisp"))
# (save-application "compile-simpx.exe" :toplevel-function #'main :prepend-kernel t) 

class TopLevel:
    def __init__(self, path_to_params_json):
        self.path_to_params_json = path_to_params_json
        self.params = json.load(open(self.path_to_params_json, 'r'))
        self.num_loops = self.params['num-loops']
        self.run_label = self.params['run-label']
        self.classifier_output = 'experiments/' + self.run_label + '/' + self.params['classifier-output']
        self.sim_output = 'experiments/' + self.run_label + '/' + self.params['sim-output']

    def prepare_results_directory(self):
        if not os.path.exists('experiments/'+self.run_label):
            os.makedirs('experiments/'+self.run_label)
            f = open(self.classifier_output, 'a')
            f.close()
            f = open(self.sim_output, 'a')
            f.close()

    def run_experiment_loop(self):
        for i in range(self.num_loops):
            print('================================================================================')
            print('Running loop:', str(i))
            os.chdir('sim')
            print('compiling sim')
            subprocess.call('./compile-simpx.exe', shell=True)
            print('running sim')
            subprocess.call('./simpx.exe', shell=True)
            os.chdir('..')
            print('running classifier')
            token = clmain(self.path_to_params_json)
            print(token)
            self.clean_us()

    def clean_us(self):
        tempfile = self.classifier_output+'.new'
        print('cleaning u-s out of ' + self.classifier_output)
        outs = open(tempfile,'w')
        for line in open(self.classifier_output):
            outs.write(re.sub(r'u\"', '\"', line))
        outs.close()
        shutil.move(tempfile,self.classifier_output)

    def cleanup(self):
        os.system('cp parameters.json experiments/' + self.run_label)
        os.system('cp classifier.py experiments/' + self.run_label)
        os.system('cp -r sim experiments/' + self.run_label)
        os.system('cp runlog-parser.py experiments/' + self.run_label)
        os.system('cp average-delta-ks.py experiments/' + self.run_label)
        os.system('cp mavg.py experiments/' + self.run_label)

if __name__ == "__main__":
    ms = TopLevel('parameters.json')
    ms.prepare_results_directory()
    ms.run_experiment_loop()
    ms.cleanup()
    

#params comments
#These params control the experiment. Make sure this file is in gctasim/RL above the experiments directory.
#Additionally, make sure that master-script and binaryclassifier.py are in that same directory.
#"run-label" - the name of the folder results from this experiment will be stored in under gctasim/RL/experiments
#"sim-output" - the file in gctasim/RL/experiments/"run-label" where the sim will write its mrhlog
#"classifier-output" - the file in gctasim/RL/experiments/"run-label" where the classifier will write its JSON output (w/ hypotheses & best-matrix)
#"num-loops" - the number of "loops" between the simulator and the classifier that will run in this experiment
#"learner" - a second JSON storing all the variables used in configuring the learner
#"simulator" - a second JSON storing all the variables used in configuring the simulator
