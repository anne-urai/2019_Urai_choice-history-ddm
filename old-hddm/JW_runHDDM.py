#!/usr/bin/env python
# encoding: utf-8
"""
Created by Jan Willem de Gee on 2011-02-16.
Adapted by Anne Urai, 2016
"""
import os, sys, pickle, time
import datetime
import math
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import itertools
import pp
from IPython import embed as shell
import hddm
import kabuki
import scipy.io # for writing to mat file
#import mne
#import statsmodels.formula.api as sm

#sys.path.append(os.environ['ANALYSIS_HOME'])
#from Tools.other_scripts import functions_jw as myfuncs

# params:
version = 0
run = True

# standard params:
model_base_name = '2ifc_MEGdata_'
model_names = ['1']
nr_samples = 20000 # 50.000 for good results
nr_models = 1 # to test
parallel = False # parallel python not working on MBP
accuracy_coding = False

# -----------------
# drift diffusion -
# -----------------
def run_model(trace_id, data, model_dir, model_name, samples=10000, accuracy_coding=False):

    import hddm
    m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True, include=('sv'),
        group_only_nodes=['sv'], depends_on={'t':['drug', 'sessionnr'],
        'v':['drug', 'sessionnr'], 'a':['drug', 'sessionnr'], 'dc':['drug', 'sessionnr'],
        'z':['drug', 'sessionnr'], }, p_outlier=.05)
    m.find_starting_values()
    m.sample(samples, burn=samples/10, thin=3,
        dbname=os.path.join(model_dir, model_name+ '_db{}'.format(trace_id)), db='pickle')
    return m

def drift_diffusion_hddm(data, samples=10000, n_jobs=6, run=True, parallel=True, model_name='model', model_dir='.', accuracy_coding=False):

    import hddm
    import os

    # run the model:
    if run:
        if parallel:
            job_server = pp.Server(ppservers=(), ncpus=n_jobs)
            start_time = time.time()
            jobs = [(trace_id, job_server.submit(run_model,(trace_id, data, model_dir, model_name, samples, accuracy_coding), (), ('hddm',))) for trace_id in range(n_jobs)]
            results = []
            shell()
            for s, job in jobs:
                results.append(job())
            print "Time elapsed: ", time.time() - start_time, "s"
            job_server.print_stats()

            # save:
            for i in range(n_jobs):
                model = results[i]
                model.save(os.path.join(model_dir, '{}_{}'.format(model_name,i)))
        else:
            start_time = time.time()
            model = run_model(3, data, model_dir, model_name, samples, accuracy_coding)
            model.save(os.path.join(model_dir, '{}_md{}'.format(model_name, 3)))

            # print point estimates
            results = model.gen_stats()
            results.to_csv(os.path.join(fig_dir, 'diagnostics', 'results3.csv'))

            # dic:
            text_file = open(os.path.join(fig_dir, 'diagnostics', 'DIC3.txt'), 'w')
            text_file.write("Model {}: {}\n".format(m, model.dic))
            text_file.close()
            print "Time elapsed: ", time.time() - start_time, "s"

    # load the models:
    else:
        print 'loading existing model(s)'
        if parallel:
            model = []
            for i in range(n_jobs):
                model.append(hddm.load(os.path.join(model_dir, '{}_{}'.format(model_name,i))))
        else:
            model = hddm.load(os.path.join(model_dir, '{}_md{}'.format(model_name, 1)))
    return model

# settings:
# ---------

# model_name:
model_name = model_names[version]

# data:
# put in 1 folder
data_path1 = os.path.join('/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM/', '2ifc_data_hddm.csv')
data = pd.read_csv(data_path1)

# model dir:
model_dir = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM/'

# figures dir:
fig_dir = os.path.join(model_dir, model_base_name + model_name)
print(fig_dir)
try:
    os.system('mkdir {}'.format(fig_dir))
    os.system('mkdir {}'.format(os.path.join(fig_dir, 'diagnostics')))
except:
    pass

# subjects:
subjects = np.unique(data.subj_idx)
nr_subjects = np.unique(data.subj_idx).shape[0]
print '# subjects = {}'.format(nr_subjects)

if run:
    print 'running {}'.format(model_base_name+model_name)
    model = drift_diffusion_hddm(data=data, samples=nr_samples, n_jobs=nr_models, run=run, parallel=parallel, model_name=model_base_name+model_name, model_dir=model_dir, accuracy_coding=accuracy_coding)
else:

    # -----------------
    # write to file
    # -----------------

  model_nr = 0

  model = drift_diffusion_hddm(data=data, samples=nr_samples, n_jobs=nr_models, run=run,
      parallel=parallel, model_name=model_base_name+model_name, model_dir=model_dir, accuracy_coding=accuracy_coding)

  params_of_interest_0 = ['z(0)', 'a(0)', 'v(0)', 'dc(0)', 't(0)', 'sv']
  params_of_interest_1 = ['z(1)', 'a(1)', 'v(1)', 'dc(1)', 't(1)', 'sv']
  params_of_interest_0s = ['z_subj(0)', 'a_subj(0)', 'v_subj(0)', 'dc_subj(0)', 't_subj(0)']
  params_of_interest_1s = ['z_subj(1)', 'a_subj(1)', 'v_subj(1)', 'dc_subj(1)', 't_subj(1)']
  titles = ['Starting point', 'Boundary sep.', 'Drift rate', 'Drift criterion', 'Non-dec. time', 'Drift rate var']

  shell()
  # point estimates:
  if parallel:
      results = model[model_nr].print_stats()
  else:
      results = model.gen_stats()
  results.to_csv(os.path.join(fig_dir, 'diagnostics', 'results.csv'))

  for i in range(nr_models):
      md = model[i]
      # remove fields that scipy io cant handle
      unwanted = [None]
      unwanted_keys = [k for k, v in md.items() if any([v is i for i in unwanted])]
      for k in unwanted_keys: del md[k]
      scipy.io.savemat(os.path.join(model_dir, '{}_{}_mat'.format(model_name,i)), md)

  shell()

  # gelman rubic:
  # only  make sense when several models were run
  gr = hddm.analyze.gelman_rubin(model)
  text_file = open(os.path.join(fig_dir, 'diagnostics', 'gelman_rubic.txt'), 'w')
  for p in gr.items():
      text_file.write("%s:%s\n" % p)
  text_file.close()

  # dic:
  text_file = open(os.path.join(fig_dir, 'diagnostics', 'DIC.txt'), 'w')
  for m in range(nr_models):
      text_file.write("Model {}: {}\n".format(m, model[m].dic))
  text_file.close()

  # # analytic plots:
  size_plot = nr_subjects / 3.0 * 1.5
  model[model_nr].plot_posterior_predictive(samples=10, bins=100, figsize=(6,size_plot), save=True, path=os.path.join(fig_dir, 'diagnostics'), format='pdf')
  model[model_nr].plot_posteriors(save=True, path=os.path.join(fig_dir, 'diagnostics'), format='pdf')


  # posterios:
  # ----------
  traces_0 = []
  traces_1 = []
  for p in range(len(params_of_interest_0)):
      traces_0.append(model[model_nr].nodes_db.node[params_of_interest_0[p]].trace.gettrace())
      traces_1.append(model[model_nr].nodes_db.node[params_of_interest_1[p]].trace.gettrace())

  # fix starting point:
  traces_0[0] = traces_0[0] * np.mean((traces_0[1].mean(),traces_1[1].mean()))
  traces_1[0] = traces_1[0] * np.mean((traces_0[1].mean(),traces_1[1].mean()))

  # # make absolute posteriors:
  # traces_0[4] = abs(traces_0[4])
  # traces_1[4] = abs(traces_1[4])
  # traces_0[5] = abs(traces_0[5])
  # traces_1[5] = abs(traces_1[5])

  # -----------------
  # plot
  # -----------------


  sns.set(style='ticks', font='Arial', font_scale=1, rc={
      'axes.linewidth': 0.25,
      'axes.labelsize': 7,
      'axes.titlesize': 7,
      'xtick.labelsize': 6,
      'ytick.labelsize': 6,
      'legend.fontsize': 6,
      'xtick.major.width': 0.25,
      'ytick.major.width': 0.25,
      'text.color': 'Black',
      'axes.labelcolor':'Black',
      'xtick.color':'Black',
      'ytick.color':'Black',} )
  sns.plotting_context()

  stats = []
  for p in range(len(params_of_interest_0)):
      data = [traces_0[p], traces_1[p]]
      stat = np.mean(data[0] > data[1])
      stats.append(min(stat, 1-stat))
  stats = np.array(stats)
  # stats_corrected = mne.stats.fdr_correction(stats, 0.05)[1]
  stats_corrected = stats
  fig, axes = plt.subplots(nrows=1, ncols=len(params_of_interest_0), figsize=(len(params_of_interest_0)*1.5,2.5))
  ax_nr = 0
  for p in range(len(params_of_interest_0)):
      data = [traces_0[p], traces_1[p]]
      ax = axes[ax_nr]
      for d, label, c in zip(data, ['low', 'high'], ['blue', 'red']):
          sns.kdeplot(d, vertical=True, shade=True, color=c, label=label, ax=ax)
          # sns.distplot(d, vertical=True, hist=False, kde_kws={"shade": True}, norm_hist=True, color=c, label=label, ax=ax)
      ax.set_xlabel('Posterior probability')
      ax.set_title(titles[p]+'\np={}'.format(round(stats_corrected[p],4)))
      ax.set_xlim(xmin=0)
      # ax.set_ylim(-1,2)
      for axis in ['top','bottom','left','right']:
          ax.spines[axis].set_linewidth(0.5)
          ax.tick_params(width=0.5)
      ax_nr+=1
  sns.despine(offset=10, trim=True)
  axes[0].set_ylabel('Parameter estimate (a.u.)')
  plt.tight_layout()
  fig.savefig(os.path.join(fig_dir, 'posteriors.pdf'))

  # import corner
  # # fig = plt.figure()
  # fig = corner.corner(np.array(traces_0).T, color='b', labels=titles, **{'lw':1})
  # corner.corner(np.array(traces_1).T, color='r', labels=titles, fig=fig, **{'lw':1})
  # fig.savefig(os.path.join(fig_dir, 'corner.pdf'))

  # #######
  # p = 5
  # data = [traces_0[p], t0[p]]
  # fig = plt.figure(figsize=(3,3))
  # ax = fig.add_subplot(111)
  # for d, label, c in zip(data, ['All trials', 'TPR fit'], ['black', 'red']):
  #     sns.kdeplot(d, vertical=True, shade=True, color=c, label=label, ax=ax)
  #     # sns.distplot(d, vertical=True, hist=False, kde_kws={"shade": True}, norm_hist=True, color=c, label=label, ax=ax)
  # ax.set_xlabel('Posterior probability')
  # ax.set_ylabel('Drift rate var')
  # ax.set_title(titles[p]+'\np={}'.format(round(np.mean(data[0] > data[1]),4)))
  # plt.tight_layout()
  # sns.despine(offset=10, trim=True)
  # fig.savefig(os.path.join(fig_dir, 'posteriors_sv.pdf'))
  #
  # barplot:
  # --------

  # all:
  parameters_h = []
  parameters_l = []
  p_value = []
  ind = np.ones(nr_subjects, dtype=bool)
  for p in range(len(params_of_interest_0s)):
      parameters_h.append(np.array([model[model_nr].values.get('{}.'.format(params_of_interest_1s[p]) + str(s)) for s in subjects])[ind])
      parameters_l.append(np.array([model[model_nr].values.get('{}.'.format(params_of_interest_0s[p]) + str(s)) for s in subjects])[ind])

  param_names = ['z', 'a', 'v', 'dc', 't']
  # param_names = ['z', 'a', 'v1', 'v2', 'dc1', 'dc2', 't']
  parameters = pd.concat((pd.DataFrame(np.vstack(parameters_h).T, columns=param_names), pd.DataFrame(np.vstack(parameters_l).T, columns=param_names)))
  parameters['pupil'] = np.concatenate((np.ones(len(subjects)), np.zeros(len(subjects))))
  parameters['subject'] = np.concatenate((subjects, subjects))
  k = parameters.groupby(['subject', 'pupil']).mean()
  k_s = k.stack().reset_index()
  k_s.columns = ['subject', 'pupil', 'param', 'value']
  parameters.to_csv(os.path.join(fig_dir, 'params.csv'))

  # plot:
  locs = np.arange(0,len(param_names))
  bar_width = 0.2
  fig = plt.figure(figsize=( (1+(len(params_of_interest_1s)*0.3)),2))
  ax = fig.add_subplot(111)
  sns.barplot(x='param',  y='value', units='subject', hue='pupil', hue_order=[1,0], data=k_s, palette=['r', 'b'], ci=None, linewidth=0, alpha=0.5, ax=ax)
  sns.stripplot(x="param", y="value", hue='pupil', hue_order=[1,0], data=k_s, jitter=False, size=2, palette=['r', 'b'], edgecolor='black', linewidth=0.25, ax=ax, split=True, alpha=1)
  for r in range(len(param_names)):
      values = np.vstack((k_s[(k_s['param'] == param_names[r]) & (k_s['pupil'] == 1)].value, k_s[(k_s['param'] == param_names[r]) & (k_s['pupil'] == 0)].value))
      x = np.array([locs[r]-bar_width, locs[r]+bar_width])
      ax.plot(x, values, color='black', lw=0.5, alpha=0.5)
  # # add p-values:
  for r in range(len(param_names)):
      p1 = myfuncs.permutationTest(k_s[(k_s['pupil']==1) & (k_s['param']==param_names[r])].value, k_s[(k_s['pupil']==0) & (k_s['param']==param_names[r])].value, paired=True)[1]
      if p1 < 0.05:
          plt.text(s='{}'.format(round(p1, 3)), x=locs[r], y=plt.gca().get_ylim()[1]-((plt.gca().get_ylim()[1] - plt.gca().get_ylim()[0]) / 10.0), size=5, horizontalalignment='center',)
  ax.legend_.remove()
  plt.xticks(locs, param_names, rotation=45)
  sns.despine(offset=10, trim=True)
  plt.tight_layout()
  fig.savefig(os.path.join(fig_dir, 'bars_all.pdf'))

  k_s = parameters.groupby(['subject', 'pupil']).mean()
  k_s = k.stack().reset_index()
  k_s.columns = ['subject', 'pupil', 'param', 'value']
  k_s = k_s[(k_s['param']=='dc')]
  param_names = ['dc']
  k_s['value'] = abs(k_s['value'])

  # plot:
  locs = np.arange(0,len(param_names))
  bar_width = 0.2
  fig = plt.figure(figsize=(1.5,2))
  ax = fig.add_subplot(111)
  sns.barplot(x='param',  y='value', units='subject', hue='pupil', hue_order=[1,0], data=k_s, palette=['r', 'b'], ci=None, linewidth=0, alpha=0.5, ax=ax)
  sns.stripplot(x="param", y="value", hue='pupil', hue_order=[1,0], data=k_s, jitter=False, size=2, palette=['r', 'b'], edgecolor='black', linewidth=0.25, ax=ax, split=True, alpha=1)
  for r in range(len(param_names)):
      values = np.vstack((k_s[(k_s['param'] == param_names[r]) & (k_s['pupil'] == 1)].value, k_s[(k_s['param'] == param_names[r]) & (k_s['pupil'] == 0)].value))
      x = np.array([locs[r]-bar_width, locs[r]+bar_width])
      ax.plot(x, values, color='black', lw=0.5, alpha=0.5)
  # # add p-values:
  for r in range(len(param_names)):
      p1 = myfuncs.permutationTest(k_s[(k_s['pupil']==1) & (k_s['param']==param_names[r])].value, k_s[(k_s['pupil']==0) & (k_s['param']==param_names[r])].value, paired=True)[1]
      if p1 < 0.05:
          plt.text(s='{}'.format(round(p1, 3)), x=locs[r], y=plt.gca().get_ylim()[1]-((plt.gca().get_ylim()[1] - plt.gca().get_ylim()[0]) / 10.0), size=5, horizontalalignment='center',)
  ax.legend_.remove()
  plt.xticks(locs, param_names, rotation=45)
  sns.despine(offset=10, trim=True)
  plt.tight_layout()
  fig.savefig(os.path.join(fig_dir, 'bars_all2.pdf'))
