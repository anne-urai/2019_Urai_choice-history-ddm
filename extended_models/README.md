
# MATLAB SCRIPTS SHARED HERE WERE IMPLEMENTED IN MATLAB R2017b #
ALSO please note that all models were fitted using parallel computing implemented in a cluster. we used a series of shell scripts  that in turn called standalone versions of the matlab scripts shared here. the shell scripts and auxiliary files used for fitting  are specific to the cluster’s implementation and they are not shared here.

The folder “extended_models” consists of 3 subfolders:

1.	data_processing
2.	default_protocol
3.	dynamic_protocol


## 1.	data_processing ##

Here, the script “descriptives_script.m” takes the raw data (motionEnergy*.mat) and summarises behaviour in RT quantiles. The output of this script is “descriptives.mat” that is subsequently used to fit different models. 


##  2.	default_protocol ##

This folder contains implementations of different models that rely on the default simulation protocol. Each subfolder corresponds to a different model. Please note that the code in the Naive and Drift_Bias folders includes detailed comments. These models can thus be used as reference to understand the other models.

The general common script structure for each model was as follows:

* “diffusion_custom2.m” is the name of the function that implements the corresponding model.

* “cost_fit.m” is the function that calculates the likelihood of a model for a given set of parameters by calling “diffusion_custom2.m.

* “superScript_Step1.m” implements the 1st step of the fitting procedure is implemented.

* “Bridge_Fits.m” reformats the fits from the 1st step such that they can be refined in the 2nd step.— “superScript_Step2.m” implements the 2nd step of the fitting procedure.


Now,

*	For all models we share the core scripts: cost_fit.m and diffusion_custom2.m.

*	For all models except the Naive, the superScript_Step1.m and superScript_Step2.m are common to those scripts in the Drift_Bias folder.

*	For the Naive model we have separate superScript_Step1.m and superScript_Step2.m scripts. Also for this model, for illustration purposes, we share the Brige_Fits.m script (common to all models).


## 3.	dynamic_protocol ##

This folder contains implementations of different models that rely on the dynamic simulation protocol. There are two subfolders implementing models with static and collapsing bounds, respectively. Similar to the default protocol, Naive and Drift_Bias scripts are commented in detail and should be use as reference.

In static_bounds/Naive we share the superScript_Step*.m files. These files are adjusted relative to the default protocol in order to implement the dynamic protocol. 

For the remaining models we share the core scripts: cost_fit.m and diffusion_custom2.m.