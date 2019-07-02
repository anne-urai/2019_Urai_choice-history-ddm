
Code and data for

#### Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, 8:e46331. ####

Behavioral data and model fits are available at https://doi.org/10.6084/m9.figshare.7268558 under a CC-BY 4.0 license.

To fit the models, install the HDDM package (http://ski.clps.brown.edu/hddm_docs/index.html). Then run the HDDM models using b1_HDDM_run.py (the models are specified in hddm_models.py). Easiest is to use a batch job submission system, and do e.g.
<code>
python b1_HDDM_run.py -r 1 -d $d -v $v -i $i -s $s
</code>
where -d = 0-6 (datasets), -v = 0-11 (versions of the model), -i = 0-30 (traces, can be changed to whatever the number of cores on a node) and -i = 5.000, the number of samples per trace.

To reproduce all main figures, see <code>plot_all.m</code> which has the overview of the scripts that are called to generate each figure. 

The <code>extended_models</code> folder contains Matlab code to fit the models in Figure 6.

The <code>simulations</code> folder contains Python code to generate Supplementary Figure 8.

For questions, @AnneEUrai / anne.urai@gmail.com.
