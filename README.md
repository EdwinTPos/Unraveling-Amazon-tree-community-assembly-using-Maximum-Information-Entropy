# Chapter 7 Rolling the dice or struggling for survival in Amazonian forests, using the maximum entropy models to disentangle drivers of community composition

R-code for running the mice imputation, normal and spatial calculations of maximum information entropy using parallel processing. A list of used packages is in the MEF.R script, reference list can be found below:

#Functions for maxent calculations are derived from the FD package and edited to operate in parallel calculation as well as spatial analysis:
Laliberté, E., and P. Legendre (2010) A distance-based framework for measuring functional diversity from multiple traits. Ecology 91:299-305.
Laliberté, E., Legendre, P., and B. Shipley. (2014). FD: measuring functional diversity from multiple traits, and other tools for functional ecology. 
R package version 1.0-12.

# Example community dataset used is from BCI (Barro Collorado Island, Panama), derived from the Vegan Package:
Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2018). vegan: Community Ecology Package. R package version 2.5-1. https://CRAN.R-project.org/package=vegan
  
# Trait dataset is derived from:
Kraft TS, Wright SJ, Turner I, Lucas PW, Oufiero CE, Noor MNS, Sun I, Dominy NJ (2015) Seed size and the evolution of leaf defences. Journal of Ecology 103(4): 1057-1068. https://doi.org/10.1111/1365-2745.12407

With Dryad Datapackage

Kraft TS, Wright SJ, Turner I, Lucas PW, Oufiero CE, Noor MNS, Sun I, Dominy NJ (2015) Data from: Seed size and the evolution of leaf defences. Dryad Digital Repository. https://doi.org/10.5061/dryad.69ph0

Reference list:
Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter
  Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2020). vegan: Community Ecology Package. R package version 2.5-7.
  https://CRAN.R-project.org/package=vegan

Laliberté, E., and P. Legendre (2010) A distance-based framework for measuring functional diversity from multiple traits. Ecology 91:299-305.

Laliberté, E., Legendre, P., and B. Shipley. (2014). FD: measuring functional diversity from multiple traits, and other tools for functional ecology. R package
  version 1.0-12.
  
Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. DOI
  10.18637/jss.v045.i03.
  
H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco Sciaini, and Cédric Scherer (2021). Rvision - Colorblind-Friendly Color Maps for R. R package
  version 0.6.2.
  
Microsoft Corporation and Steve Weston (2020). doParallel: Foreach Parallel Adaptor for the 'parallel' Package. R package version 1.0.16.
  https://CRAN.R-project.org/package=doParallel
  
Douglas Nychka, Reinhard Furrer, John Paige, Stephan Sain (2021). “fields: Tools for spatial data.” R package version 13.3, <URL:
   https://github.com/dnychka/fieldsRPackage>.

Steve Weston and Hadley Wickham (2014). itertools: Iterator Tools. R package version 0.1-3. https://CRAN.R-project.org/package=itertools


