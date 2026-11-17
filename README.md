### Lastest realease
[![DOI](https://zenodo.org/badge/1461651817.svg)](https://doi.org/10.5281/zenodo.24320259)

### R Folder structure

```
R_folder
├── data                             # contains the data
│   ├── raw                          # contains raw datasets
│   │   ├── effort                   # contains EURING effort per country, for ringing data used in the migration atlas
│   │   ├── PECBMS_national_2023     # contains PECBMS national species trend data until 2023 "PECBMS_indices_till 2023"
│   │   ├── pdf_migatlas_transproba  # contains the pdfs with the transition probabilities from the Bird Migration Atlas
│   │   └── supplementary            # contains data used for the supplementary materials (the area of each of the 25 countries of the present analysis)
│   └── processed                    # dataset processed
│
├── figures  
│   ├── supplementary       # All supplementary figures
│   │    └── svg_files      # supplementary svg files
│   └── svg_files           # svg files for inkscape (for example) editing
│        ├── abundance      # to store the svg files for figures on abundance
│        ├── percentage     # to store the svg files for figures on percentage of migratory individuals
│        └── pop_index      # to store the svg files for figures on population index figures
│
├── functions               # contains some functions used in the analysis, essentially adapted from Rigal et al. 2020 and Rigal et al. 2023 
│ 
├── outputs
│   ├── rds_objects         # to store .rds objects
│   ├── results             # to store some results
│   ├── RTRIM_default       # to store RTRIM outputs of the main analysis, contains "All", "Migrants", "LD", "SD", and "R"
│   └── … # other RTRIM folders for supplementaries
│
└── scripts                 # contains all the scripts used in this analysis (except the functions, in the "functions" folder
```


### Instructions for coding

#### Data availability

- *PECBMS data* was obtained on request, named "PECBMS_indices_till 2023.csv". We do not provide the dataset.
- *ERLoB abundance data* is openly accessible on the Article 12 web tool on population status and trends of birds under Article 12 of the Birds Directive at https://nature-art12.eionet.europa.eu/article12/summary. We do not provide the dataset.
- *Clements list* is openly accessible on the Birds of the World webiste at https://www.birds.cornell.edu/clementschecklist. We do not provide the dataset.
- EURING effort per country dataset is openly accessible on EURING's website https://euring.org/data-and-codes/euring-databank. We do not provide the dataset.
- *countries_regions* assigns countries to the regions, used for transition probabilities correction as well as mapping. We provide the dataset.
- *area_countries* compiles the area of the 25 countries of this analysis, excluding significant inland water bodies. We provide the dataset.
- transition probabilities data from Fattorini et al. (2023) are openly accessible on migratory connectivity section of the Birds of the World website (https://www.migrationatlas.org). Species pdfs have to be downloaded separately, and to be named after the species name (actuallised according to the Clements list). We do not provide the dataset.
- *order_species.csv* associates species to their order, it is used for filtering landbird species. The base used is the AVONET dataset from witch we retained only the "Species" and "Order1" columns.


#### Scripts

The following sequence has to be respected, as outputs from one script may be input for later scripts.

- *00_mapping_settings* sets the area of interest, extract countries *sf* geometries and rename some geographic objects.
- *01_index&abundance_sp_country* homogeneises species names across datasets, merge ERLoB abundance data and PECBMS trend data, and defines landbirds.
- *02_scale_analysis_sp_time_space* rafines the spatio-temporal scale of the analysis. It is also used to produce supplementary information on spatio-temporal coverage, country comparisons for abundances, density, as well as country-specific stacked area charts of species abundance through time (not presented in the supplementary materials).
- *03_extract_transition_probability_EABMA* automatizes the extraction of transition probabilities from the species-specific pdfs from the Bird Migration Atlas, compiling these in a "transition_table.csv" file.
- *04_correction_transition_probability* compiles EURING data on ringing totals per country to obtain a relative importance of each region in the ringing data used for the Bird Migration Atlas, and correct the transition probabilities from "*transition_table.csv*" according to the estimated relative effort per region.
- *05_define_migratory_stategies* associates migratory strategies, diviging national populations in sub-populations according to the corrected transition rates.
- *06_multi_species_trends* prepares the dataset to run the *Nonlinear_functions_rigal2023* to obtain national multi-species index and total abundance for each migratory strategy. Associated tables and maps (Figure 3 and Figure 4) are created using the same script. National community composition (% of individuals of each migration strategy, Figure 1) and the associated changes (change in % of migratory individuals, Figure 5) are also obtained on this script. Sensitivity analyses can be fitted by modifying the "*type*" and "*is_supplementary*" parameters at the begining of the script. Other supplementary materials are also produced in this script (abundances of missing populations for each country, detail of national change in multispecies population index and total abundance for each migratory strategy). 
- *07_RTRIM* procudes 1000 multi-country species trends for each species, using the *RTRIM_functions_rigal2023* functions.
- *08_MSI_Europe* uses as an input the outputs trends from *07_RTRIM*, and calculates continental multispecies index and total abundance for each migratory strategy with the *Nonlinear_functions_rigal2023_adapted_EUROPE* functions. The Figure 2 is produced in this script, as well as the supplementary on total abundance change regardless of migration strategy.
- *09_trend_sp_country* creates the dataset associating species national population change and migration strategy. The supplementary figures comparing for each country the first year of monitoring (as early as 1990) to 2023 population sizes for the PECBMS populations is created here.


#### Functions

- *classtrajectory* to classify population trends, function from Rigal et al. 2020
- *extract_transition_probabilities_from_pdfs* to extract transition probabilities table from the pdfs of Fattorini et al. (2023) in the Bird Migration Atlas. Note that the pdfs have to be named after the species name actuallised according to the Clements list.
- *graph_framework* sets the function to produce the maps
- *Nonlinear_functions_rigal2023* to obtain national multi-species index and national total abundance for each migratory strategy. Script from Rigal et al. (2023). We addapted parts of the code to fit the goal of this study, the corresponding addings being between and after this symbol "# ~~~~~~". This function may take some time to run.
- *Nonlinear_functions_rigal2023_adapted_EUROPE* to obtain continental multi-species index and continental total abundance for each migratory strategy. This is the same script as *Nonlinear_functions_rigal2023*, but modified to use as an input the 1000 RTRIM iterations. This may take some time to run. This function may take some time to run.
- *RTRIM_functions_rigal2023* to obtain 1000 country-aggregated trends for each species. When the script runs, outputs are first stored in "intermediary results", then the population trend output of the iteration is stored in "SI.imputed" and the "intermediary results" outputs are removed to save place on the disk.

#### R environment
```
R version 4.3.2 (2023-10-31 ucrt) -- "Eye Holes"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)
```

### README Bibliography
```
Fattorini, Niccolò, Alessandra Costanzo, Andrea Romano, et al. 2023. « Eco-Evolutionary Drivers of Avian Migratory Connectivity ». Ecology Letters 26 (7): 1095‑107. https://doi.org/10.1111/ele.14223.

Rigal, Stanislas, Vincent Devictor, et Vasilis Dakos. 2020. « A Method for Classifying and Comparing Non-Linear Trajectories of Ecological Variables ». Ecological Indicators 112 (mai): 106113. https://doi.org/10.1016/j.ecolind.2020.106113.

Rigal, Stanislas, Vasilis Dakos, Hany Alonso, et al. 2023. « Farmland practices are driving bird population decline across Europe ». Proceedings of the National Academy of Sciences 120 (21): e2216573120. https://doi.org/10.1073/pnas.2216573120.

Tobias, Joseph A., Catherine Sheard, Alex L. Pigot, Adam JM Devenish, Jingyi Yang, Ferran Sayol, Montague HC Neate‐Clegg et al. 2022. « AVONET: morphological, ecological and geographical data for all birds ». Ecology letters 25 (3): 581-597.
