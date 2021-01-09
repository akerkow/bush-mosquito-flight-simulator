# bush-mosquito-flight-simulator
This repository contains the model scripts for the Linkage of a compartment model for West-Nile virus with a flight simulator for vecor mosquitoes

The aim of the work presented here was to provide an estimate of the local spread of the virus after its introduction to a new location through the movements of mosquitoes over time and space. For this purpose, we adapted an existing SEIR model for West Nile virus to the conditions in Germany (temperatures, geographical latitude, bird and mosquito species densities) and the characteristic transmission and life trait parameter of a possible host bird and vector mosquito species. We further extended it by a spatial component: an agent based flight simulator for the mosquitoes. It demonstrates how the female mosquitoes move within the landscape due to habitat structures and wind conditions and how many of them leave the region in the different cardinal directions.

The model was designed by me and Dr. Ralf Wieland at the Leibniz Centre for Agricultural Landscape Research (ZALF) within the framework of the the CuliFo-project. The project was supported by the Federal Ministry of Food and Agriculture (BMEL, grant no. 2819105315).

A manuscript for the publication of the model is in preparation. The model code is subject to the "Creative Commons Attribution - Non Commercial 4.0" licence agreement.

The model is implemented in Python 3.7 and the different components (main simulation, temporal and spatial component as well as the import of temperature data and visualisation of results) are stored in individual files to obtain a good overview. In order to conveniently determine relevant parameters for the scenario analysis, e.g. the mosquito species, the study region or the time period and duration for a model application, these parameters are stored in a control file and get called by the other model files: 

- To run the model, "main_simu.py" must be executed. The script uses the SEIR compartment model ("time_simu.py") and the flight simulator ("space_simu.py"). 

- The SEIR model needs weather data, prepared for each simulation day. This information is provided by the script "read_temperature.py". 

- Graphs showing the dependence of different model parameters on temperature, day length and day of the year can be generated using the script "visualize_time_simu.py".

- Parameters for scenario analyses can be set in "Model_Scenario.py". 

- Model results can be visualized using the "visualize_log" file. These include the trends in mosquito and bird populations and infection numbers over the course of the year. In addition, the spatial distribution of mosquitoes in the study regions at the time with the most infections is shown, as well as the number of mosquitoes that attempted to leave the study region in the different cardinal directions.

- Also, videos can be created showing the movements of all mosquitoes ("animate_NM.py") or of the infected mosquitoes ("animate_IM.py") during the simulation period. Note that the videos are created using numpy arrays generated during the model run. For model applications for periods longer than one year, no numpy arrays are stored by default. If the video animation fails and you are working with Spyder, try downgrading Spyder to a version lower than 4.0.
