#!/usr/bin/env python
# encoding: utf-8

import os
import numpy as np
import scipy as sp

import matplotlib as mpl
mpl.use('Agg') # to still plot even when no display is defined
import matplotlib.pyplot as plt

mpl.rcParams['pdf.fonttype'] = 42
# import matplotlib.pylab as plt
import seaborn as sns
import pandas as pd
#import bottleneck as bn
from IPython import embed as shell
import hddm
import pymc as pm

sns.set(style='ticks', font='Arial', font_scale=1, rc={
    'axes.linewidth': 0.25, 
    'axes.labelsize': 8, 
    'axes.titlesize': 7, 
    'xtick.labelsize': 6, 
    'ytick.labelsize': 6, 
    'legend.fontsize': 6, 
    'xtick.major.width': 0.1, 
    'ytick.major.width': 0.1,
    'text.color': 'Black',
    'axes.labelcolor':'Black',
    'xtick.color':'Black',
    'ytick.color':'Black',} )
sns.plotting_context()

def plot_all_priors(model, data=None, unique=True, model_kwargs=None):
	"""
	plot the priors of an HDDM model
	Input:
		data <DataFrame> - data to be plot against the priors
		unique <bool> - whether to unique each column in data before before ploting it
	"""

	#set limits for plots
	lb = {'v': -10, 'dc(1)':-5, 'z':0.001, 'z_std':0}
	ub = {'a': 4, 't':1, 'v':10, 'z':1, 'sz': 1, 'st':1, 'sv':15, 'p_outlier': 1,
        'z_trans(1)':1, 'z(1)':1, 'dc(1)':5, 'a_std':5, 'v_std':5, 'z_std':0.5, 't_std':5, 'dc_std':5}

	#plot all priors
	n_rows=4
	n_cols=5
	for n_subjs in [1]: #,2]:

		# create a model
        # h_data, _ = hddm.generate.gen_rand_data(subjs=n_subjs, size=2)
        # if model_kwargs is None:
        #     model_kwargs = {}
        # h = model(h_data, include='all', **model_kwargs)
        
        #h = model

		fig = plt.figure()
                plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=.7)

		counter = 0
		for name, node_row in model.iter_group_nodes():
			if not name in ub: # only those listed
			    continue
			if 'var' in name or 'p_outlier' in name:
				continue
			if 'trans' in name:
				trans = True
				name = name.replace('_trans','')
			else:
				trans = False
			counter += 1
			node = node_row['node']
			print(name)
			print(node.logp)

			#plot a single proir
			ax = plt.subplot(n_rows, n_cols, counter)
                        ax.set_yticklabels([])
		
			#generate pdf
			xlim = np.arange(lb.get(name, 0.001), ub[name], 0.01)
			pdf = np.zeros(len(xlim))
            # assume that the logp has the prior?
			for i in range(len(pdf)):
				if not trans:
					node.value = xlim[i]
					pdf[i] = np.exp(node.logp)
				else:
					node.value = pm.logit(xlim[i])
					pdf[i] = np.exp(node.logp)*10

			#plot shit
			plt.plot(xlim, pdf)
			plt.xlabel(name)
			sns.despine(offset=2, trim=True)
            
			# # Hide the right and top spines
#             ax.spines['right'].set_visible(False)
#             ax.spines['top'].set_visible(False)
#
#             # Only show ticks on the left and bottom spines
#             ax.yaxis.set_ticks_position('left')
#             ax.xaxis.set_ticks_position('bottom')


		#add suptitle
		plt.suptitle('HDDM priors')
            
        # save the figure
        plt.savefig(os.path.join(mypath, 'priorPlot.pdf'))
            
## LOAD MODEL WITH THE MOST PARAMETERS WE HAVE
mypath = os.path.realpath(os.path.expanduser('/nfs/aeurai/HDDM/JW_PNAS'))
m = hddm.load(os.path.join(mypath, 'stimcoding_dc_z_prevresp_st', 'modelfit-combined.model'))
#print(m)
#shell()
plot_all_priors(m)


