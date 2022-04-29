# How to run
# Make sure that the script for your classifier is in the same dir as this
# script. simpx.lisp should be in a subdir called "sim". paramters.json should
# also be in this directory.
# Then, run:
# $python scikit-master-script.py

#default python libs
import re
import os
import json
import subprocess
import shutil

#classifier used with this experiment...
from scikitclassifier import main as clmain
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
        self.timestamp = self.find_latest_timestamp()
        self.path_to_this_exp = 'experiments/'+self.timestamp+'/'+self.run_label
        self.classifier_output = self.path_to_this_exp+'/'+self.params['classifier-output']
        self.sim_output = self.path_to_this_exp+'/'+self.params['sim-output']
        self.to_copy = ['parameters.json', 'classifier.py', 'sim', 'runlog-parser.py', 'mavg.py']

    def find_latest_timestamp(self):
        dirs = os.listdir('experiments')
        maxdir = 0
        for d in dirs:
            try:
                d = int(d)
                if d > maxdir:
                    maxdir = d
            except:
                pass
        return str(maxdir)

    def prepare_results_directory(self):
        os.makedirs(self.path_to_this_exp)
        f = open(self.classifier_output, 'w')
        f.close()
        f = open(self.sim_output, 'w')
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
            token = clmain(self.path_to_params_json, self.timestamp)
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
        for obj in self.to_copy:
            os.system('cp -r ' + obj + ' ' + self.path_to_this_exp)

if __name__ == "__main__":
    ms = TopLevel('parameters.json')
    ms.prepare_results_directory()
    ms.run_experiment_loop()
    ms.cleanup()
