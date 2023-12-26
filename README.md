# NEX-GDDP-FWI
The NEX-GDDP-FWI data was created from the python codes shared in this repository. 

The NASA Earth eXchange-Global Daily Downscaled Projections-Fire Weather Index (NEX-GDDP-FWI) is a dataset that provides global fire weather projections derived from daily Global Climate Model (GCM) simulations. The dataset uses the Canadian Forest Fire Weather Index System (CFWIS) framework to estimate fire danger by considering the effects of fuel moisture and wind on fire behavior and spread. The dataset includes retrospective (1950-2014) and prospective (2015-2100) simulations from 33 GCMs participating in the Coupled Model Intercomparison Project Phase 6 (CMIP6) under four Shared Socio-economic Pathways (SSPs). To make the dataset more accessible, the fire weather metrics are summarized at coarser temporal scales (monthly and annually), and source codes are provided for investigating daily fire weather. Additionally, Multi-Model Ensemble (MME) datasets are provided for each fire weather metric, which include monthly and annually summarized fire weather metrics. The total data volume of the dataset is approximately 4.5 TB and is available through the NASA Advanced Supercomputing (NAS) Data Portal.

Here are brief explanations for the python codes:

FWIfunctions.py: This python code has all sub-component functions used in the CFWIS framework.

calFireIndex.py: This python code is a main code reading input NEX-GDDP data and compute sub-components of the CFWIS framework.

calFireIndex_baseline.py: This python code is a code reading input NEX-GDDP data and compute sub-components of the CFWIS framework to produce historical baseline data.

get_MME_FWI_Monthly_release.py: This python code is for summarizing daily fire weather indices of individual GCM to multi-model ensemble at monthly time step 

get_MME_FWI_Yearly_release.py: This python code is for summarizing daily fire weather indices of individual GCM to multi-model ensemble at yearly time step 
