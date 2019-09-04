# Fossil Uncertainty
These are the files necessary for creating (and hopefully recreating) my project on accounting for fossil uncertainty in downstream macroevolutionary analyses, using macropodoid marsupials as an example.

If you download this whole zip, you should have everything you need. You'll obviously need to change file paths to make them match your machine. 

# Modelling Analyses
All analyses are detailed in the "Macropodoid_Modelling" files (.pdf, .html, .Rmd). I've included them in a number of formats with the hope that you can use them however suits you. 

# Tree Files:
  ### Part of this project was estimating trees under a number of different dating schemes. They're listed below:
   
   + **Tree_Span.trees** - a sample of trees from posterior of all dating methods which show range of maximum and minimum estimated ages
   + **Macropodinae_Fossil.trees** - a sample of trees from all posterior dating methods (similar to above), which also include fossil taxa lacking molecular/morphological data. These taxa are  constrained by age/topological priors in the starBEAST analyses.
   + **Macropodinae_MaxAges** - post-burnin sample of trees estimated using maximum ages for fossil taxa.
   + **Macropodinae_MeanAges** - post-burnin sample of trees estimated using mean ages for fossil taxa.
   + **Macropodinae_MinAges** - post-burnin sample of trees estimated using minimum ages for fossil taxa.
   + **Macropodinae_SampledAges** - post-burnin sample of trees estimated using tip ages sampled from uniform priors bounded by minimum and maximum fossil ages.
   + **Macropodinae_SampledAges_PriorOnly** - post-burnin sample of trees estimated from the priors only, using tip ages sampled from uniform priors bounded by minimum and maximum fossil ages.


