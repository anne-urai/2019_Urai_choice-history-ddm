
Code and data for

#### Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595. ####

The summary csv files, including history information, are all in the Data subfolder and have been named using the abbreviations described in the paper.

To reproduce the figures, install the HDDM package (http://ski.clps.brown.edu/hddm_docs/index.html). Then run the HDDM models using b1_HDDM_run.py (the models are specified in hddm_models.py). Easiest is to use a batch job submission system, and do e.g.
<code>
python b1_HDDM_run.py -r 1 -d $d -v $v -i $i -s $s
</code>
where -d = 0-6 (datasets), -v = 0-11 (versions of the model), -i = 0-30 (traces, can be changed to whatever the number of cores on a node) and -i = 5.000, the number of samples per trace.

Then, in Matlab run the file <code>plotAll.m</code>, which will read in the models and reproduce all figures step-by-step.

For questions, @AnneEUrai / anne.urai@gmail.com.
