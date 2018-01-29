Leaky accumulator scripts:

--This model generates behaviour on a task involving detection of transient 'signals' embedded in dynamic noise, as specified in Ossmy et al. (2013, Current Biology, http://www.cell.com/current-biology/abstract/S0960-9822(13)00480-6). Model requires four parameters per cell of experimental design:
	-brightness exponent (re-scales stimulus input; specific to this task)
	-accumulator noise
	-decision bound
	-accumulamtor leak

--Highest-level script is 'FIT_2acc.m', which expects two inputs: a subject ID (for loading behavioural data) and a parameter constraint string (the scripts assume we're fitting to data from only two experimental conditions). To understand steps involved in fitting procedure I recommend starting here and working down through lower-level scripts as they're called.

--The key dynamics of the model play out in 'sim2acc.m' - this simulates the dynamics of two independent leaky accumulators using slow Monte Carlo methods (i.e. simulating each time-step of each trial). These methods are used because there's no known analytical solution for generating behaviour from this type of model.

--For rationale behind using this specific fitting method for this model, I recommend looking at the SI in the Ossmy et al. paper.