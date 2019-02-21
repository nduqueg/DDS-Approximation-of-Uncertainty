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
+ When the evaluation of the objective function (i.e. the total model simulations), takes an amount of time that is unfeasible and a reduction in the number of evaluations would cause a decline in the model performance.
+ When the number of dimensions in the calibration problem impose a challenge on the uniform random sampling to find parameter sets near the maximum likelihood solution.
+ Basically, these two criterias are met by distributed environmental models at daily or shorter time step.

Routines
==============================
The script includes four functions to perform the sample of behavioural parameters, the clasiffication of them in non-behavioural and behavioural, and finally execute a post-processing so the uncertainty of the simulated variables can be quantified. These routines are called:
+ dds_au, which calls the independent DDS algorithms (algorithm which is programed in a different fuction and can be executed independently to just calibrate a environmental model as state by Tolson, B. A., &amp; Shoemaker, C. A., (2007)). The dds_au function automatically classifies the results of the independent DDS in behavioural or non behavioural taking into account a threshold assigned by the user. The algorithm follows the recomedations made by Tolson, B. A., &amp; Shoemaker, C. A., (2008) regarding the aleatory assigment of the iterations to each DDS.
+ reclas_threshold, which recieves the output of a dds_au process and once more identifies the behavioural parameters with a different threshold.
+ post_dds_au, which only performs the simulations of the behavioural parameters and calculates several uncertainty band metrics, using the observations, the bounds of the bands and the mean of the simulations.

Inputs
====================
+ xBounds.df must be a dataframe with 1st column as minimum and 2nd column as maximum of the parameter range.
+ numIter is an integer which defines the total of model simulations.
+ numBeh is an integer which defines the number of independent DDS and the maximum number of behavioural threshold.
+ per.m.dds is a double which indicates the percentage that affects the computational budget of each independent DDS.
+ threshold is a double which classifies the parameter sets in non-behavioural or behavioural.
+ r is a double between 0 and 1, the default value is 0.2 (DDS parameter of perturbation).
+ OBJFUN is a function which returns a scalar value, for which we are trying to minimize.
+ obs is a vector containing the streamflow observations
+ dates is a dataframe containing the dates of begining and ending of the observations and simulations, and the time interval

Outputs
====================
A list that cointains two more list. The first one contains the behavioural sets of parameters and their respective objective function value, and the second one cointains all the simulations classified by the independent DDS number.

Requirements
====================
zoo and ggplot libraries must be installed before performing the post_dds_au script
