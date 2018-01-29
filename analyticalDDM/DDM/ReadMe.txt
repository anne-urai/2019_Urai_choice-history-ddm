Drift diffusion model scripts:

--These generate correct/error RT densities via the DDM for any speeded 2AFC task, using the methods described in Ratcliff & Tuerlinckx (2002; https://link.springer.com/article/10.3758/BF03196302). These particular scripts are built to take a maximum of four parameters per cell of experimental design:
	-drift rate
	-non-decision time
	-decision bound
	-drift rate variability

--Note that in the full DDM there are additional parameters (starting point, starting point variability, nondecision time variability), but I have not incorporated these here.

--Note also that I built these scripts to be a counterpart to more complicated model fits involving time-dependent decision bounds, decribed here: https://www.nature.com/articles/ncomms13526. As such, there a few steps that are a little unorthodox but were used to match the other fitting approach in this paper. I've done my best to highlight these in the code.

--Highest-level script is 'FIT_regular_DDM_PSO.m', which expects three inputs: a subject ID (for loading behavioural data), a run number (not so relevant) and a parameter constraint string (these scripts assume we're fitting to data from only two experimental conditions). To understand steps involved in fitting procedure I recommend starting with this high-level script and working down through lower-level scripts as they're called.

--The key dynamics of the model play out in 'fpt_regular_DDM.m' - this generates first passage time densitites (i.e. RT distributions) using fast(ish) analytical methods, since the analytical solution for the first passage time problem posed by the DDM is known (in contrast to the LCA).

--Here I use a particle swarm optimization procedure for parameter estimation - this requires many more iterations that, e.g. fminsearch, but I'd trust it more. The reasons we can use it here is because the simulation methods are much faster than our Monte Carlo methods for the LCA. I've attached the PSO code (taken from https://uk.mathworks.com/matlabcentral/fileexchange/7506-particle-swarm-optimization-toolbox and very slightly adapted) in a sub-directory.