#!/usr/bin/env python

"""
By Anne Urai, 2016
"""

import hddm, os, pickle, sys, kabuki
from kabuki import experiments
from IPython import embed # for debugging, drop in ipython



# Create different models to test our various hypotheses.
experiments = [
   {'name': 'basic', 'data': mydata,
    'model_type': 'hddm.HDDMStimCoding', 
    'kwargs': {'depends_on': {},
    'stim_col':'stimulus', 'split_param':'v'}}]

print "Running experiments..."
kabuki.experiments.run_experiments(experiments, subj_map_init=False,
    mpi=True, {'samples':50000, 'burn':10000, 'thin':3, 'db':'pickle'}):

print "Analyzing experiments..."
kabuki.experiments.analyze_experiments(experiments, ppc=False)
print "Done! Check the newly created subdirectories."
