### Quartet of interactions
This repository contains codes used to estimate the quartet of interactions between species pairs using genome-scale metabolic modeling. The method employs genome-scale metabolic models to predict the growth rates of the microbial species when simulated alone and in communities. Using the predicted growth-resource curves of each species, it estimates the positive and negative components.  

The method is demonstrated by employing two different algorithms for estimating the growth rate in communities, SteadyCom (MATLAB) and MICOM (Python). 

MATLAB code:

The MATLAB code estimates the quartet of interactions between species pairs with net interactions determined by the SteadyCom algorithm. 

Pre-requisites:
1.	Installation of COBRA toolbox
2.	SteadyCom code downloaded

Input:
1.	Genome-scale models of the two species as a cell array
2.	Resource supply level or richness (parameter that determines the proportion of available nutrients)

Output:
1.	Net interactions between species pairs
2.	Associated quartet components

PYTHON code:

The PYTHON code estimates the quartet of interactions between species pairs with net interactions determined by the MICOM algorithm. 

Pre-requisites:
1.	Installation of COBRA toolbox
2.	Installation of MICOM

Input:
1.	Path of the folder that contains the models of the species (The folder should contain a sub-folder termed ‘Models’ that contains the genome-scale models of the species names model_1 and model_2)
2.	Resource supply level or richness (parameter that determines the proportion of available nutrients)

Output:
1.	Net interactions between species pairs
2.	Associated quartet components
