# ReceptorNetworks
MATLAB code for simulating dynamics of receptor networks
========================

MATLAB library developed to produce the results for the paper 

[Stanoev, Angel, Akhilesh P. Nandan, and Aneta Koseska. "Optimal biochemical information processing at criticality." BioRxiv (2019): 543348.](https://www.biorxiv.org/content/10.1101/543348v1.full.pdf).

-------------------------
Requirements
-------------------------

The code has been tested successfully with older versions of MATLAB (R2016a (9.0)) and with the latest version (R2019b (9.7)) on Windows.

Using the framework
===================

Good starting point is running the figure geneting scripts.
-----------------------------------------------------------

In general, for running simulations of a customized model, it needs to be defined first, similar to the way the other models are defined (labels, df_model function, etc).
Same holds for the experiment type (pulsed, sustained, step-wise increasing/decreasing, etc.).
Then the two objects are combined in a model_simulation object that can run deterministic/stochastic simulations and plot results.
Sample code:

```matlab
model = models.simple_dnf_model;  % create model object
model.par.g1 = 2.957;  % set regime of operation by setting the parameter(s)
mpe = experiments.multi_pulse_experiment(2);  % create experiment object
ms = model_simulation(model, mpe);  % combine the two in a model_simulation object 
ms.simulate();  % run deterministic simulations
ms.plot_fraction_phosphorylated();  % plot the results
```