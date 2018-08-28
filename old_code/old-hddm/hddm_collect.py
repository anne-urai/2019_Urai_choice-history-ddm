#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2016
adapted from JW de Gee
"""

import hddm, os
nr_traces   = 3
model_name  = 'basic_stimcoding'

# find path depending on local/klimag
usr = os.environ.get('USER')
if usr in ['anne']:
    mypath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM'
if usr in ['aurai']:
    mypath = '/home/aurai/Data/MEG-PL/Data/HDDM'
thispath    = os.path.join(mypath, model_name)

print "appending models"
models = []
for t in range(nr_traces): # run the models serially
    models.append(hddm.load(os.path.join(thispath, 'modelfit-md%d'%t)))

print "computing gelman-rubin convergence statistics"
# compute gelman rubic
gr = hddm.analyze.gelman_rubin(models)
text_file = open(os.path.join(thispath, 'gelman_rubic.txt'), 'w')
for p in gr.items():
     text_file.write("%s:%s\n" % p)
text_file.close()
print "DONE!"
