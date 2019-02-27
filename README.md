# Fossil Uncertainty
These are the files necessary for creating (and hopefully recreating) my project on accounting for fossil uncertainty in downstream macroevolutionary analyses, using macropodoid marsupials as an example.

If you download this whole zip, you should have everything you need. You'll obviously need to change file paths to make them match your machine. 

Otherwise, this repository includes the following (organized under folders of the same names):  
## Data:
  ### Trees:  
   + Maximum clade credibility trees. (5x)  
   + 100 trees from dating analysis posteriors. (5x) 
   + 100 trees from extinction regime exercises. (3x: Stochastic, Miocene, PlioPleistocene)  
  ### Phenotypic and Distributional Traits
   + Body size data (snout-vent length, body length, mass). (5x)
   + Distributional data, in binary format for BioGeoBEARS analysis. (5x)
   + Spreadsheet of raw body size and distributional data.

# Modelling Analyses
All analyses are detailed in the "Macropodoid_Modelling" files (.pdf, .html, .Rmd). I've included them in a number of formats with the hope that you can use them however suits you. 

# Tree Files:
  ### Part of this project was estimating trees under a number of different dating schemes. They're listed below:
   + Macro_Ages_Ranges... - trees estimated with tip age priors
   + Macro_CP_EstAges... - trees estimated with tip ages as estimated by Couzens & Prideaux
   + Macro_CP_MinAges... - trees estimated with minimum tip ages of Couzens & Prideaux
   + Macro_CP_MaxAges... - trees estimated with maximum tip ages of Couzens & Prideaux
   + Macro_MinAges... - trees estimated with minimum stratigraphic tip ages
   + Macro_MaxAges... - trees estimated with maximum stratigraphic tip ages
   + Tree_Span.tre - trees from posterior of all dating methods which show range of maximum and minimum estimated ages
