# NEX-GDDP-FWI
The NEX-GDDP-FWI data was created from the python codes shared in this repository. 

Finer scale fire weather projections enable more effective planning, resource allocation, mitigation efforts, and public awareness to protect lives, property, and ecosystems in the face of increasing wildfire threats. Here, we report on a dataset that provides global fire weather projections derived from downscaled (0.25°) and bias-corrected daily Earth System Model (ESM) simulations. The dataset uses the Canadian Forest Fire Weather Index System (CFWIS) framework to estimate fire danger by considering the effects of fuel moisture and wind on fire behavior and spread. The dataset includes retrospective (1950-2014) and prospective (2015-2100) simulations from 33 ESMs participating in the Coupled Model Intercomparison Project Phase 6 (CMIP6) using four Shared Socio-economic Pathways (SSPs). To make the dataset more accessible, the fire weather metrics are summarized at coarser temporal scales (monthly and annually), and source codes are provided for investigating daily fire weather. Multi-Model Ensemble (MME) data of monthly and annual fire weather metrics are also provided. This publicly available 5.1 TB dataset has the potential to be broadly used in not only for wildfire risk assessment but also for various future climate change impact assessments and preparedness.


Here are brief explanations for the python codes:

FWIfunctions.py: This python code has all sub-component functions used in the CFWIS framework.

calFireIndex.py: This python code is a main code reading input NEX-GDDP data and compute sub-components of the CFWIS framework.

calFireIndex_baseline.py: This python code is a code reading input NEX-GDDP data and compute sub-components of the CFWIS framework to produce historical baseline data.

get_MME_FWI_Monthly_release.py: This python code is for summarizing daily fire weather indices of individual GCM to multi-model ensemble at monthly time step 

get_MME_FWI_Yearly_release.py: This python code is for summarizing daily fire weather indices of individual GCM to multi-model ensemble at yearly time step 
