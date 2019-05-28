# Ric-Spe-code
Code of the GEB paper "Recent global changes have decoupled species richness from specialization patterns in North American birds"

# Authors and affiliations
Anne Mimet 1,2,3*, Robert Buitenwerf 1,4, Brody Sandel 1,5, Jens-Christian Svenning 1,4 and Signe Normand 1,4

*Anne MIMET (corresponding author): +49 3412351970; anne.mimet@gmail.com
1 Section for EcoInformatics & Biodiversity, Department of Bioscience, Aarhus University, Ny Munkegade 114, 8000 Aarhus C, Denmark
2 Department Computational Landscape Ecology, Helmholtz Center for Environmental Research, Permoserstraße 15, 04318 Leipzig, Leipzig, Germany
3 Biodiversity Conservation Group, German Center for Integrative Biodiversity Research (iDiv) Halle-Jena-Leipzig, Deutscher Platz 5e, 04103 Leipzig, Germany 
4 Center for Biodiversity Dynamics in a Changing World (BIOCHANGE), Department of Bioscience, Aarhus University, Ny Munkegade 114, 8000 Aarhus C, Denmark
5 Biology Department, Santa Clara University, 500 El Camino Real, Santa Clara, CA 95053, United States

# Available code files
We provide 5 different code files corresponding to the different steps and analyses presented in the paper. These code files can be run using the data made freely available with the paper (link to be updated).

File 1 - computing_ssi.R: Code for defining the habitats (multivariate regression tree approach), computing the species specialization and the specialization of the communities (Steps 1, 2 and 3 in the paper). This code also create the data base with community richness and specialization used for the analyses (Community_data.txt)

File 2 - GLMM_plus_modnull.R - Code for analysis (1): Temporal change in the richness-specialization relationship. Include the code of the GLMM and of the null models + figures

File 3 - GAMs_annual_model_comparison.R - Code for analysis (2): Spatial covariation of the richness and specialization and its temporal change at the annual time scale. Code for the comparison of annual GAMs (linear, piecewise linear regression with one breakpoint, piecewise linear regression with two breakpoints) + figures

File 4 - GAMs_decades_model_comparison.R - Code for analysis (3): Spatial covariation of the richness and specialization and its temporal change at the decadal time scale. Code for the comparison of decadal GAMs (linear, piecewise linear regression with one breakpoint, piecewise linear regression with two breakpoints) + figures 

File 5 - statico_analysis.R - Code for analysis (4): Spatio-temporal environmental drivers of spatial pattern of the richness-specialization relationship. Code for the STATICO analyses and associated figures
