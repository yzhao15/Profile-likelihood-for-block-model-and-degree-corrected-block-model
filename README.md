# Profile-likelihood-for-block-models-and-degree-corrected-block-models
Codes for “Consistency of community detection in networks under degree-corrected stochastic block models”

We provides the codes which implement community detection method using the profile-likelihoods of block models and degree-corrected block models. The profile-likelihoods are optimized by a local label-switching technique, tabu-search. 

libBlockModel.R is for block models, where ``tabuDiffStart'' is the main function with two parameters: A-- the adjacency matrix, K-- the number of communities. 

libCorr.R is for degree-corrected block models, where ``tabuDiffStart'' is the main function with two parameters: A-- the adjacency matrix, K-- the number of communities. 
