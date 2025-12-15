# Description of modelling workflow and scripts used in *Adaptation potential of local water management practices to glacier loss*
Scripts used to develop the modelling workflow used in the article "Adaptation potential of local water management practices to glacier loss". The scripts include post-processing of JULES land-surface model outputs, combining with glacier evolution simulations, routing model, deep sub-surface separation and routing, and three interventions models.

Work builds on Simon Moulds workflow (Copyright (c) 2021, simonmoulds ;https://github.com/simonmoulds/jules-mosart/tree/main).

The scripts included are used after running vn6.1. JULES land-surface model using rose-suite configuration u/c/j/5/3/1 (https://code.metoffice.gove.uk/trac/roses-u/browser/c/j/5/3/1). The workflow combines JULES simulations outputs with JULES-OGGM glacier evolution simulations described in Mackay et al. (2025). Both models are run with hourly meteorological data from WRF bias-corrected simulations during a historical period 1980-2018 (Potter et al., 2023a), and from an ensemble of 30 statistically downscaled CMIP5 models for RCP4.5 and RCP8.5 emission scenarios for the period 2019-2100.

## Workflow order
The workflow is to be run in the following order:
1. JULES post-processing and NbS model modifications to JULES outputs.
3. Surface and sub-surface routing.
4. Storage equivalent analysis.

A full workflow template bash script (full_workflow_template.sh) is provided and designed to be run at Imperial's HPC where parallelization is done with a PBS routine (default in Imperial's HPC).

## References
J. D. Mackay, N. E. Barrand, D. M. Hannah, E. Potter, N. Montoya, W. Buytaert, Physically based modelling of glacier evolution under climate change in the tropical Andes. The Cryosphere 19, 685–712 (2025).

E. Potter, C. Fyffe, A. Orr, D. Quincey, A. Ross, S. Rangecroft, K. Medina, H. Burns, A. Llacza, G. Jacome, R. Hellstrom, J. Castro, A. Cochachin, N. Montoya, E. Loarte, F. Pellicciotti, Bias-corrected temperature and precipitation data from the WRF regional climate model output, Cordillera Blanca and Vilcanota-Urubamba regions, Peru, from 1980 to 2018, version 1.0, NERC EDS UK Polar Data Centre (2023a); https://doi.org/10.5285/2CF25580-9B79-440F-8505-6230DD377877.

E. Potter, C. Fyffe, A. Orr, D. Quincey, A. Ross, S. Rangecroft, K. Medina, H. Burns, A. Llacza, G. Jacome, R. Hellstrom, J. Castro, A. Cochachin, N. Montoya, E. Loarte, F. Pellicciotti, Precipitation and temperature data from statistically downscaled CMIP5 models, Cordillera Blanca and Vilcanota-Urubamba regions, Peru, from 2019 to 2100, version 1.0, NERC EDS UK Polar Data Centre (2023b); https://doi.org/10.5285/67CEB7C8-218C-46E1-9927-CFEF2DD95526.


