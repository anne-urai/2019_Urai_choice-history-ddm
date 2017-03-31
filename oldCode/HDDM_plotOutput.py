#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2016
takes input arguments from stopos
Important: on Cartesius, call module load python2.7.9 before running
(the only environment where HDDM is installed)
"""

# ============================================ #
# HDDM cheat sheet
# ============================================ #

# v     = drift rate
# a     = boundary separation
# t     = nondecision time
# z     = starting point
# dc    = drift driterion
# sv    = inter-trial variability in drift-rate
# st    = inter-trial variability in non-decision time
# sz    = inter-trial variability in starting-point

# ============================================ #
# parse input arguments
# ============================================ #

from optparse import OptionParser
usage = "runHDDM_stopos.py [options]"
parser = OptionParser ( usage)
parser.add_option ( "-v", "--version",
        default = 1,
        type = "int",
        help = "Version of the model to run" )
parser.add_option ( "-i", "--trace_id",
        default = 1,
        type = "int",
        help = "Which trace to run, usually 1-3" )
opts,args       = parser.parse_args()
model_version   = opts.version
trace_id        = opts.trace_id

# to avoid errors when plotting on cartesius
# http://stackoverflow.com/questions/4706451/how-to-save-a-figure-remotely-with-pylab/4706614#4706614
import matplotlib
matplotlib.use('Agg') # to still plot even when no display is defined
import matplotlib.pyplot as plt

# general stuff
import numpy as np
import os, hddm, time

# ============================================ #
# my own version of posterior predictive plot
# ============================================ #

def plot_posterior_predictive_anne(m, path, **kwargs):

    # overlay predicted onto real RT distributions
    from kabuki.analyze import _plot_posterior_pdf_node

    plot_func = _plot_posterior_pdf_node
    observeds = m.get_observeds()
    fig = plt.figure(figsize=[24,33]) # make a huge one so the fonts are not so big

    # Plot different conditions (new figure for each)
    for tag, nodes in observeds.groupby('tag'):
        fig.clf()
        fig.suptitle(tag, fontsize=12)
        fig.subplots_adjust(top=0.9, hspace=.4, wspace=.3)

        # Plot individual subjects (if present)
        for subj_i, (node_name, bottom_node) in enumerate(nodes.iterrows()):
            ax = fig.add_subplot(np.ceil(np.sqrt(len(nodes))), np.ceil(np.sqrt(len(nodes))), subj_i+1)

            # http://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')

            if 'subj_idx' in bottom_node:
                ax.set_title(str(bottom_node['subj_idx']))
            plot_func(bottom_node['node'], ax, value_range=np.linspace(-4,4,100), **kwargs)
        fig.savefig('%s.pdf' %os.path.join(path, 'ppd_%s'%(tag, )))

def interpolate_trace(x, trace, range=(-1,1), bins=100):
    import scipy.interpolate
    x_histo = np.linspace(range[0], range[1], bins)
    histo = np.histogram(trace, bins=bins, range=range, density=True)[0]
    interp = scipy.interpolate.InterpolatedUnivariateSpline(x_histo, histo)(x)
    return interp

def plot_posterior_nodes_anne(nodes, bins=50, lb=None, ub=None):

    fig = plt.figure(figsize=[8,11]) # make a huge one so the fonts are not so big
    if lb is None:
        lb = min([min(node.trace()[:]) for node in nodes])
    if ub is None:
        ub = max([max(node.trace()[:]) for node in nodes])
    x_data = np.linspace(lb, ub, 300)

    for node in nodes:
        trace = node.trace()[:]
        hist = interpolate_trace(x_data, trace, range=(lb, ub), bins=bins)
        plt.plot(x_data, hist, label=node.__name__, lw=2.)

    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    return fig

# ============================================ #
# define the function that will do the work
# ============================================ #

def plot_model(mypath, model_name, trace_id):

    # load in the model that was ran
    m = hddm.load(os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id))

    # ============================================ #
    # save plots
    # ============================================ #

    # plot some output stuff in figures subfolder
    figpath = os.path.join(mypath, model_name, 'figures-md%d'%trace_id)
    if not os.path.exists(figpath):
        os.mkdir(figpath)

    # 1. plot the traces and posteriors for each parameter
    m.plot_posteriors(save=True, path=figpath, format='pdf')

    # 2. plot posterior predictive
    print('plotting posterior predictive')
    plot_posterior_predictive_anne(m, path=figpath)

    # 3. plot the actual posteriors and the way they depend on the variables we specified
    if model_name in ['prevresp_prevrt_dc', 'prevresp_prevpupil_dc']:
        print('plotting the posteriors by previous response and rt')
        dc_prevresp0_prevRTlow, dc_prevresp0_prevRTmed, dc_prevresp0_prevRThigh, \
            dc_prevresp1_prevRTlow, dc_prevresp1_prevRTmed, dc_prevresp1_prevRThigh \
            = m.nodes_db.node[['dc(-1.0.1.0)', 'dc(-1.0.2.0)', 'dc(-1.0.3.0)', 'dc(1.0.1.0)','dc(1.0.2.0)','dc(1.0.3.0)']]

        # plot these myself
        plot_posterior_nodes_anne([dc_prevresp0_prevRTlow, dc_prevresp0_prevRTmed, dc_prevresp0_prevRThigh, \
            dc_prevresp1_prevRTlow, dc_prevresp1_prevRTmed, dc_prevresp1_prevRThigh])

        plt.xlabel('Drift criterion')
        plt.ylabel('Posterior probability')
        plt.title('Posterior of drift-rate group means')
        plt.savefig(os.path.join(figpath, 'driftcriterion_posteriors.pdf'))

    elif model_name in ['prevresp_prevrt_z', 'prevresp_prevpupil_z']:
        print('plotting the posteriors by previous response and rt')
        dc_prevresp0_prevRTlow, dc_prevresp0_prevRTmed, dc_prevresp0_prevRThigh, \
            dc_prevresp1_prevRTlow, dc_prevresp1_prevRTmed, dc_prevresp1_prevRThigh \
            = m.nodes_db.node[['z(-1.0.1.0)', 'z(-1.0.2.0)', 'z(-1.0.3.0)', 'z(1.0.1.0)','z(1.0.2.0)','z(1.0.3.0)']]

        # plot these myself
        plot_posterior_nodes_anne([dc_prevresp0_prevRTlow, dc_prevresp0_prevRTmed, dc_prevresp0_prevRThigh, \
            dc_prevresp1_prevRTlow, dc_prevresp1_prevRTmed, dc_prevresp1_prevRThigh])

        plt.xlabel('Starting point')
        plt.ylabel('Posterior probability')
        plt.title('Posterior of drift-rate group means')
        plt.savefig(os.path.join(figpath, 'startingpoint_posteriors.pdf'))

    elif model_name in ['prevresp_z']:
        print('plotting the posteriors by previous response')
        dc_prevresp0, dc_prevresp1 \
            = m.nodes_db.node[['z(-1.0)', 'z(1.0)']]

        # plot these myself
        plot_posterior_nodes_anne([dc_prevresp0, dc_prevresp1])
        plt.xlabel('Starting point')
        plt.ylabel('Posterior probability')
        plt.title('Posterior of drift-rate group means')
        plt.savefig(os.path.join(figpath, 'startingpoint_posteriors.pdf'))

    elif model_name in ['prevresp_dc']:
        print('plotting the posteriors by previous response')
        dc_prevresp0, dc_prevresp1 \
            = m.nodes_db.node[['dc(-1.0)', 'dc(1.0)']]

        # plot these myself
        plot_posterior_nodes_anne([dc_prevresp0, dc_prevresp1])
        plt.xlabel('Drift criterion')
        plt.ylabel('Posterior probability')
        plt.title('Posterior of drift-rate group means')
        plt.savefig(os.path.join(figpath, 'driftcriterion_posteriors.pdf'))

# ============================================ #
# run one model per job
# ============================================ #

# which model are we running at the moment?
models = {0: 'stimcoding',
    1: 'prevresp_z',
    2: 'prevresp_dc',
    3: 'prevresp_prevrt_z',
    4: 'prevresp_prevrt_dc',
    5: 'prevresp_prevpupil_z',
    6: 'prevresp_prevpupil_dc'}

# find path depending on location
usr = os.environ.get('USER')

if usr in ['anne']: #local
    mypath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM'
if usr in ['aurai']: #uke cluster
    mypath = '/home/aurai/Data/MEG-PL/Data/HDDM'
if usr in ['aeurai']: #cartesius
    mypath = '/home/aeurai/neurodec/Data/MEG-PL/HDDM'

# make a folder for the outputs, combine name and time
thispath = os.path.join(mypath, models[model_version])
if not os.path.exists(thispath):
    os.mkdir(thispath)

# and... go!
plot_model(mypath, models[model_version], trace_id)
