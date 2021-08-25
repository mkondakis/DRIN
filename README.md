# DRIN
R project for modelling ecological developmental rates as temperature varies. 
we use the Bayesian paradigm in conjunction with specific data generation models. 
The latter models include not only various data distributions such as the widely used Gaussian, 
Inverse Gamma, and zero inflated Inverse Gamma, which can also generate excess zeros, 
but also four popular ecological non-linear functions that describe developmental rate dynamics. 
Nonetheless, there are significant challenges to overcome, 
such as the function structure dependence on thermal parameters, 
parameter ranges, and initial value selection when applying these ecological functions. 
The various ecological models used are the Bieri, the Briere, the Analytis and the Lactin. 
This package allows for a variety of initial values, prior distributions, and predictions, 
which adds adaptability to user choice.
