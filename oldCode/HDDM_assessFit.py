
# ============================================ #
# post-processing
# ============================================ #

import hddm
import matplotlib.pyplot as plt

print "HDDM imported, starting post-processing"
models = []
for trace_id in range(nr_traces): # run the models serially
    thism = hddm.load(os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id))
    print os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)

    # plot some output stuff in figures subfolder
    figpath = os.path.join(mypath, model_name, 'figures-md%d'%trace_id)
    if not os.path.exists(figpath):
        os.mkdir(figpath)
    thism.plot_posteriors(save=True, path=figpath, format='pdf')
    plt.close('all') # this will leave figures open, make sure to close them all
    models.append(thism)

# gelman rubic on the list of models
gr = hddm.analyze.gelman_rubin(models)
text_file = open(os.path.join(mypath, model_name, 'gelman_rubic.txt'), 'w')
for p in gr.items():
     text_file.write("%s:%s\n" % p)
text_file.close()

# ============================================ #
# plot posteriors
# ============================================ #
