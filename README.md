# Selection for infectivity profiles in slow and fast epidemics

This is the code and data used to investigate the selection for SARS-CoV-2 variants depending on their infectivity profile.

## Codes

### simulation_study_v3.R

Conducts the simulation study to verify inference.

### inference_UK_v3.R

Conducts the inference on UK data.

### get_clean_data.R

Cleans the data on cases number and frequency dynamics.

### functions_v3.R

Auxiliary functions.

### inference_functions.R

Functions needed for inference.

### conceptual_figures_v2.R

Draw conceptual figures of the paper.

### simulate_v2.cpp

Code doing the epidemiological simulations in cpp

### source_simulation_parameters.R

R code sourcing the parameters for simulations.

### spatial_correlation.R

Code cleaning the spatial frequency of Delta across EU countries.

### inference_results_v2.RData

Results of the inference method on Alpha and Delta variant data.

## Data sources

### cleaned_UK_cases.csv

Clean UK cases used for inference. These data are available from https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDate&format=csv

### cleaned_UK_frequencies.csv

Clean UK SGTF frequencies used for inference. These data are available from https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
Data underlying technical briefing 9.

### cleaned_UK_frequencies_period2.csv

Clean UK SGTF frequencies used for inference. These data are available from https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
Data underlying technical briefing 17.

### delta_frequencies_EU.csv

Delta frequencies across EU countries and across dates. These data is available from the ECDC website at https://www.ecdc.europa.eu/en/publications-data/data-virus-variants-covid-19-eueea

### spatial_variation_EU.csv

Cleaned delta frequencies across EU countries at the time when the frequency passes 50%, used for inference.


