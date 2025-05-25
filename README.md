### Disentangling the quartet of interactions
The method describes the _in silico_ to identify the quartet of interactions among species pairs. The approach disentangles the net interactions into positive and negative components. The method employs genome-scale metabolic models to predict the growth rates of the microbial species when simulated in communities. Using the predicted growth-resource curves of each species, we estimate the competitive and cooperative components, assuming the predominance of metabolic interactions. Thus, each species pair is defined by a quartet of interactions, the competitive and cooperative interactions of each species with the other. 

The method is demonstrated by employing two different algorithms for estimating the growth rate in communities, SteadyCom (MATLAB) and MICOM (Python). 

MATLAB code:
The MATLAB code estimates the quartet of interactions between species pairs with net interactions determined by SteadyCom algorithm. 
Pre-requisites:
1.	Installation of COBRA toolbox
2.	SteadyCom code downloaded

Input:
1.	Genome-scale models of the two species as a cell array
2.	Common richness (parameter that determines the proportion of available nutrients. Please refer to the manuscript for further information)
Output:
1.	Net interactions among species pairs
2.	Quartet of interactions among species pairs

PYTHON code:
The PYTHON code estimates the quartet of interactions between species pairs with net interactions determined by MICOM algorithm. 
Pre-requisites:
1.	Installation of COBRA toolbox
2.	Installation of MICOM

Input:
1.	Path of the folder that contains the models of the species (The folder should contain a sub-folder termed ‘Models’ that contains the genome-scale models of the species names model_1 and model_2)
2.	Common richness (parameter that determines the proportion of available nutrients. Please refer to the manuscript for further information)

Output:
1.	Net interactions among species pairs
2.	Quartet of interactions among species pairs
