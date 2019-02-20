# DDS-Aproximation-of-Uncertainty
An implementation in R of the DDS-AU uncertainty estimation algorithm, developed by Tolson &amp; Shoemaker (2008).
Questions relating to this implementation of the algorithm can be directed to nduqueg@unal.edu.co. 
This code is made available for free use with the condition that Tolson, B. A., &amp; Shoemaker, C. A., (2008) and Duque-Gardeazabal N. &amp; Fuentes C. (2019) will be cited.

DDS-AU
==============================
The DDS-AU algortihm must be cited as Tolson, B. A., and C. A. Shoemaker (2008), Efficient prediction uncertainty approximation in the calibration ofenvironmental simulation models, Water Resour. Res., 44, W04411, doi:10.1029/2007WR005869.

This implementation and its analysis must be cited as Duque-Gardeazabal N. &amp; Fuentes C. (2019), Using R to easily and efficiently predict the uncertainty in simu-lations of environmental models, Revista Hidrolatinoamericana de Jóvenes Investigadores y Profesionales, vol. 3.

WHY SHOULD YOU USE DDS-AU?
------------------
The analysis of the uncertainty bounds in environmental models is a demanding task, due to the computational budget that needs to be provided for the calculation of the agreement between the simulations and observations. Within the GLUE framework, high-quality model parameters must be sample in order to characterize the uncertainty bounds. However, the uniform random sampling used in Monte Carlo experiments is rather inefficient, and particularly, when high performance parameter sets are searched in a several dimensional space. DDS-AU seeks to tackle that inefficiency of sampling by using a number of independent DDS search algorithms.

Therefore, you should use DDS-AU to estimate the parameter uncertainty in these cases:
+ When the model evaluations takes an amount of time that is unfeasible.
+ When the number of dimensions in the calibration problem impose a challenge on the uniform random sampling to find parameter sets near the maximum likelihood solution.

Inputs
====================


Outputs
====================


