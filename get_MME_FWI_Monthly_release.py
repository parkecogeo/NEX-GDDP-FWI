# This code is for calculating fwi multimodel ensemble median across all CMIP6 models. 

# Taejin Park
# 2023-10-10
# taejin1392@gmail.com


from re import I
import xarray as xr
import numpy as np
import warnings
import datetime
import pandas as pd
from FWIfunctions import *
from timeit import default_timer as timer
from datetime import timedelta
import os, glob, sys
import uuid
import pytz

#suppress warnings
warnings.filterwarnings('ignore')


root_Y = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/Yearly/'
root_M = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/Monthly/'
outroot_Y = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/Yearly_MME/'
outroot_M = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/Monthly_MME/'

#root_Y = '/Users/taejinpark/Projects/Data/GDDP2/FWI/Yearly/'
#root_M = '/Users/taejinpark/Projects/Data/GDDP2/FWI/Monthly/'
#outroot_Y = '/Users/taejinpark/Projects/Data/GDDP2/FWI/Yearly_MME/'
#outroot_M = '/Users/taejinpark/Projects/Data/GDDP2/FWI/Monthly_MME/'

if not os.path.exists(outroot_Y):
    os.makedirs(outroot_Y)

if not os.path.exists(outroot_M):
    os.makedirs(outroot_M)


Models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', \
    'CESM2-WACCM', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CNRM-CM6-1', 'CNRM-ESM2-1', \
    'EC-Earth3', 'EC-Earth3-Veg-LR', 'FGOALS-g3', 'GFDL-CM4', 'GFDL-CM4_gr2', \
    'GFDL-ESM4', 'GISS-E2-1-G', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'IITM-ESM', \
    'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'KACE-1-0-G', 'KIOST-ESM', \
    'MIROC6', 'MIROC-ES2L', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', \
    'NESM3', 'NorESM2-LM', 'NorESM2-MM', 'TaiESM1', 'UKESM1-0-LL'] # total 35 models

runtypes = ['historical','ssp126','ssp245','ssp370','ssp585']# total 5 run types


#-------------------------------------------------------------------------------
def get_Monthly_MME(runtypeid,inputid):

    irun = runtypes[runtypeid]

    if irun == 'historical':
        syear = 1950
        eyear = 2014
    else:
        syear = 2015
        eyear = 2100

    years = list(range(syear,eyear+1))
    #iyear = years[inputid]
    iyear = inputid

    print('iyear = %d'%(iyear))
    
    fn = '%s*/*_%s_*_%s.nc'%(root_M,irun,iyear)
    fn_models = glob.glob(fn)
    print(fn_models)

    
    icount = 0
    for imodel in fn_models:
        print('imodel = %s'%(imodel))
        ds = np.squeeze(xr.open_dataset(imodel))
        ds = ds.assign_coords({'modelid':icount})
        ds = ds.expand_dims(dim={"modelid": 1})
        print(ds)
        monthtime = pd.date_range('%s-01-01'%(iyear),'%s-12-31'%(iyear),freq='MS').tolist()
        ds = ds.drop_vars('time')
        ds = ds.assign_coords({'time':monthtime})
        ds.attrs=''

        if icount==0:
            ds_stack = ds
        else:
            ds_stack = xr.combine_by_coords([ds_stack, ds])

        icount += 1
        
    ds_stack_P75 = ds_stack.quantile(0.75, dim="modelid")
    ds_stack_P50 = ds_stack.quantile(0.50, dim="modelid")
    ds_stack_P25 = ds_stack.quantile(0.25, dim="modelid")
    ds_stack_P75 = ds_stack_P75.drop_vars('quantile')
    ds_stack_P50 = ds_stack_P50.drop_vars('quantile')
    ds_stack_P25 = ds_stack_P25.drop_vars('quantile')

    ds_stack_P25['FWI_N15'] = ds_stack_P25['FWI_N15'].astype('uint16')
    ds_stack_P25['FWI_N30'] = ds_stack_P25['FWI_N30'].astype('uint16')
    ds_stack_P25['FWI_N45'] = ds_stack_P25['FWI_N45'].astype('uint16')
    ds_stack_P25['FWI_Nmid'] = ds_stack_P25['FWI_Nmid'].astype('uint16')
    ds_stack_P25['FWI_NP95'] = ds_stack_P25['FWI_NP95'].astype('uint16')

    ds_stack_P50['FWI_N15'] = ds_stack_P50['FWI_N15'].astype('uint16')
    ds_stack_P50['FWI_N30'] = ds_stack_P50['FWI_N30'].astype('uint16')
    ds_stack_P50['FWI_N45'] = ds_stack_P50['FWI_N45'].astype('uint16')
    ds_stack_P50['FWI_Nmid'] = ds_stack_P50['FWI_Nmid'].astype('uint16')
    ds_stack_P50['FWI_NP95'] = ds_stack_P50['FWI_NP95'].astype('uint16')

    ds_stack_P75['FWI_N15'] = ds_stack_P75['FWI_N15'].astype('uint16')
    ds_stack_P75['FWI_N30'] = ds_stack_P75['FWI_N30'].astype('uint16')
    ds_stack_P75['FWI_N45'] = ds_stack_P75['FWI_N45'].astype('uint16')
    ds_stack_P75['FWI_Nmid'] = ds_stack_P75['FWI_Nmid'].astype('uint16')
    ds_stack_P75['FWI_NP95'] = ds_stack_P75['FWI_NP95'].astype('uint16')



    # Global attribute information
    Att_Activity = 'NEX-GDDP-FWI'
    Att_Contact = 'Dr. Ian G. Brosnan: ian.g.brosnan@nasa.gov, Dr. Taejin Park: taejin.park@nasa.gov'
    Att_Frequency_D = 'day'
    Att_Frequency_M = 'month'
    Att_Frequency_Y = 'year'
    Att_Institution = 'NASA Earth Exchange, NASA Ames Research Center, Moffett Field, CA 94035'
    Att_Version = '1.0'
    Att_product = 'output'
    Att_convention = 'CF-1.7'
    Att_Disclaimer = 'This data is considered provisional and subject to change. This data is provided as is without any warranty of any kind, either express or implied, arising by law or otherwise, including but not limited to warranties of completeness, non-infringement, accuracy, merchantability, or fitness for a particular purpose. The user assumes all risk associated with the use of, or inability to use, this data.'
    Att_scenario = irun
    Att_tracking_id_D = str(uuid.uuid4())
    Att_tracking_id_M = str(uuid.uuid4())
    Att_tracking_id_Y = str(uuid.uuid4())
    Att_resolution_id = '0.25 degree'
    Att_creation_date = datetime.datetime.now(pytz.timezone('UTC')) 
    Att_history = '%s %s'%(Att_creation_date,'Version 1.0 production')
    Att_title_P25 = '%s, %s, %s'%('Multi-Model Ensemble 25th Percentile',Att_scenario,'global fire weather index data from downscaled CMIP6 projection data')
    Att_title_P50 = '%s, %s, %s'%('Multi-Model Ensemble 50th Percentile',Att_scenario,'global fire weather index data from downscaled CMIP6 projection data')
    Att_title_P75 = '%s, %s, %s'%('Multi-Model Ensemble 75th Percentile',Att_scenario,'global fire weather index data from downscaled CMIP6 projection data')
    if irun == 'historical':
        Att_model_list = 'ACCESS-CM2, ACCESS-ESM1-5, CanESM5, CESM2, CESM2-WACCM, CMCC-CM2-SR5, CMCC-ESM2, CNRM-CM6-1, CNRM-ESM2-1, EC-Earth3, EC-Earth3-Veg-LR, FGOALS-g3, GFDL-CM4, GFDL-CM4_gr2, GFDL-ESM4, GISS-E2-1-G, HadGEM3-GC31-LL, HadGEM3-GC31-MM, IITM-ESM, INM-CM4-8, INM-CM5-0, IPSL-CM6A-LR, KACE-1-0-G, KIOST-ESM, MIROC6, MIROC-ES2L, MPI-ESM1-2-HR, MPI-ESM1-2-LR, MRI-ESM2-0, NorESM2-LM, NorESM2-MM, TaiESM1, UKESM1-0-LL'
    elif irun == 'ssp126':
        Att_model_list = 'ACCESS-CM2, ACCESS-ESM1-5, CanESM5, CESM2, CMCC-CM2-SR5, CMCC-ESM2, CNRM-CM6-1, CNRM-ESM2-1, EC-Earth3, EC-Earth3-Veg-LR, FGOALS-g3, GFDL-ESM4, GISS-E2-1-G, HadGEM3-GC31-LL, HadGEM3-GC31-MM, IITM-ESM, INM-CM4-8, INM-CM5-0, IPSL-CM6A-LR, KACE-1-0-G, KIOST-ESM, MIROC6, MIROC-ES2L, MPI-ESM1-2-HR, MPI-ESM1-2-LR, MRI-ESM2-0, NorESM2-LM, NorESM2-MM, TaiESM1, UKESM1-0-LL'
    elif irun == 'ssp245':
        Att_model_list = 'ACCESS-CM2, ACCESS-ESM1-5, CanESM5, CESM2, CESM2-WACCM, CMCC-CM2-SR5, CMCC-ESM2, CNRM-CM6-1, CNRM-ESM2-1, EC-Earth3, EC-Earth3-Veg-LR, FGOALS-g3, GFDL-CM4, GFDL-CM4_gr2, GFDL-ESM4, GISS-E2-1-G, HadGEM3-GC31-LL, IITM-ESM, INM-CM4-8, INM-CM5-0, IPSL-CM6A-LR, KACE-1-0-G, KIOST-ESM, MIROC6, MIROC-ES2L, MPI-ESM1-2-HR, MPI-ESM1-2-LR, MRI-ESM2-0, NorESM2-LM, NorESM2-MM, TaiESM1, UKESM1-0-LL'
    elif irun == 'ssp370':       
        Att_model_list = 'ACCESS-CM2, ACCESS-ESM1-5, CanESM5, CESM2, CMCC-CM2-SR5, CMCC-ESM2, CNRM-CM6-1, CNRM-ESM2-1, EC-Earth3, EC-Earth3-Veg-LR, FGOALS-g3, GFDL-ESM4, GISS-E2-1-G, IITM-ESM, INM-CM4-8, INM-CM5-0, IPSL-CM6A-LR, KACE-1-0-G, MIROC6, MIROC-ES2L, MPI-ESM1-2-HR, MPI-ESM1-2-LR, MRI-ESM2-0, NorESM2-LM, NorESM2-MM, TaiESM1, UKESM1-0-LL'
    else:
        Att_model_list = 'ACCESS-CM2, ACCESS-ESM1-5, CanESM5, CESM2, CESM2-WACCM, CMCC-CM2-SR5, CMCC-ESM2, CNRM-CM6-1, CNRM-ESM2-1, EC-Earth3, EC-Earth3-Veg-LR, FGOALS-g3, GFDL-CM4, GFDL-CM4_gr2, GFDL-ESM4, GISS-E2-1-G, HadGEM3-GC31-LL, HadGEM3-GC31-MM, IITM-ESM, INM-CM4-8, INM-CM5-0, IPSL-CM6A-LR, KACE-1-0-G, KIOST-ESM, MIROC6, MIROC-ES2L, MPI-ESM1-2-HR, MPI-ESM1-2-LR, MRI-ESM2-0, NorESM2-LM, NorESM2-MM, TaiESM1, UKESM1-0-LL'
    Att_reference ='NEX-GDDP-FWI description paper: Park, T. et al., Under Review. Downscaled 21st century global fire weather projections. Nature Scientific Data // NEX-GDDP description paper: Thrasher, B. et al., 2022. NASA global daily downscaled projections, CMIP6. Nature Scientific data (https://doi.org/10.1038/s41597-022-01393-4) // Canadian Fire Weather Index System: Van Wagner, C.E., 1987. Development and structure of the Canadian forest fire weather index system, Canadian Forest Service (http://cfs.nrcan.gc.ca/pubwarehouse/pdfs/19927.pdf)'

    ds_stack_P25.FFMC.attrs['fullname'] = 'monthly mean FFMC (Fine Fuel Moisture Code)'
    ds_stack_P25.DMC.attrs['fullname'] = 'monthly mean DMC (Duff Moisture Code)'
    ds_stack_P25.DC.attrs['fullname'] = 'monthly mean DC (Drought Code)'
    ds_stack_P25.ISI.attrs['fullname'] = 'monthly mean ISI (Initial Spread Index)'
    ds_stack_P25.BUI.attrs['fullname'] = 'monthly mean BUI (Buildup Index)'
    ds_stack_P25.FWI.attrs['fullname'] = 'monthly mean FWI (Fire Weather Index)'
    ds_stack_P25.FWI_Nmid.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > midpoint of reference period'
    ds_stack_P25.FWI_NP95.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 95 percentile of reference period'
    ds_stack_P25.FWI_N15.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 15'
    ds_stack_P25.FWI_N30.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 30'
    ds_stack_P25.FWI_N45.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 45'
    ds_stack_P25.attrs['activity'] = '%s'%(Att_Activity)
    ds_stack_P25.attrs['institution'] = '%s'%(Att_Institution)
    ds_stack_P25.attrs['contact'] = '%s'%(Att_Contact)
    ds_stack_P25.attrs['scenario'] = '%s'%(Att_scenario)
    ds_stack_P25.attrs['resolution_id'] = '%s'%(Att_resolution_id)
    ds_stack_P25.attrs['frequency'] = '%s'%(Att_Frequency_M)
    ds_stack_P25.attrs['version'] = '%s'%(Att_Version)
    ds_stack_P25.attrs['title'] = '%s'%(Att_title_P25)
    ds_stack_P25.attrs['model_list'] = '%s'%(Att_model_list)
    ds_stack_P25.attrs['reference'] = '%s'%(Att_reference)
    ds_stack_P25.attrs['creation_date'] = '%s'%(Att_creation_date)
    ds_stack_P25.attrs['convention'] = '%s'%(Att_convention)
    ds_stack_P25.attrs['tracking_id'] = '%s'%(Att_tracking_id_M),
    ds_stack_P25.attrs['history'] = '%s'%(Att_history)
    ds_stack_P25.attrs['disclaimer'] = '%s'%(Att_Disclaimer)


    ds_stack_P50.FFMC.attrs['fullname'] = 'monthly mean FFMC (Fine Fuel Moisture Code)'
    ds_stack_P50.DMC.attrs['fullname'] = 'monthly mean DMC (Duff Moisture Code)'
    ds_stack_P50.DC.attrs['fullname'] = 'monthly mean DC (Drought Code)'
    ds_stack_P50.ISI.attrs['fullname'] = 'monthly mean ISI (Initial Spread Index)'
    ds_stack_P50.BUI.attrs['fullname'] = 'monthly mean BUI (Buildup Index)'
    ds_stack_P50.FWI.attrs['fullname'] = 'monthly mean FWI (Fire Weather Index)'
    ds_stack_P50.FWI_Nmid.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > midpoint of reference period'
    ds_stack_P50.FWI_NP95.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 95 percentile of reference period'
    ds_stack_P50.FWI_N15.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 15'
    ds_stack_P50.FWI_N30.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 30'
    ds_stack_P50.FWI_N45.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 45'
    ds_stack_P50.attrs['activity'] = '%s'%(Att_Activity)
    ds_stack_P50.attrs['institution'] = '%s'%(Att_Institution)
    ds_stack_P50.attrs['contact'] = '%s'%(Att_Contact)
    ds_stack_P50.attrs['scenario'] = '%s'%(Att_scenario)
    ds_stack_P50.attrs['resolution_id'] = '%s'%(Att_resolution_id)
    ds_stack_P50.attrs['frequency'] = '%s'%(Att_Frequency_M)
    ds_stack_P50.attrs['version'] = '%s'%(Att_Version)
    ds_stack_P50.attrs['title'] = '%s'%(Att_title_P50)
    ds_stack_P50.attrs['model_list'] = '%s'%(Att_model_list)
    ds_stack_P50.attrs['reference'] = '%s'%(Att_reference)
    ds_stack_P50.attrs['creation_date'] = '%s'%(Att_creation_date)
    ds_stack_P50.attrs['convention'] = '%s'%(Att_convention)
    ds_stack_P50.attrs['tracking_id'] = '%s'%(Att_tracking_id_M),
    ds_stack_P50.attrs['history'] = '%s'%(Att_history)
    ds_stack_P50.attrs['disclaimer'] = '%s'%(Att_Disclaimer)


    ds_stack_P75.FFMC.attrs['fullname'] = 'monthly mean FFMC (Fine Fuel Moisture Code)'
    ds_stack_P75.DMC.attrs['fullname'] = 'monthly mean DMC (Duff Moisture Code)'
    ds_stack_P75.DC.attrs['fullname'] = 'monthly mean DC (Drought Code)'
    ds_stack_P75.ISI.attrs['fullname'] = 'monthly mean ISI (Initial Spread Index)'
    ds_stack_P75.BUI.attrs['fullname'] = 'monthly mean BUI (Buildup Index)'
    ds_stack_P75.FWI.attrs['fullname'] = 'monthly mean FWI (Fire Weather Index)'
    ds_stack_P75.FWI_Nmid.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > midpoint of reference period'
    ds_stack_P75.FWI_NP95.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 95 percentile of reference period'
    ds_stack_P75.FWI_N15.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 15'
    ds_stack_P75.FWI_N30.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 30'
    ds_stack_P75.FWI_N45.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 45'
    ds_stack_P75.attrs['activity'] = '%s'%(Att_Activity)
    ds_stack_P75.attrs['institution'] = '%s'%(Att_Institution)
    ds_stack_P75.attrs['contact'] = '%s'%(Att_Contact)
    ds_stack_P75.attrs['scenario'] = '%s'%(Att_scenario)
    ds_stack_P75.attrs['resolution_id'] = '%s'%(Att_resolution_id)
    ds_stack_P75.attrs['frequency'] = '%s'%(Att_Frequency_M)
    ds_stack_P75.attrs['version'] = '%s'%(Att_Version)
    ds_stack_P75.attrs['title'] = '%s'%(Att_title_P75)
    ds_stack_P75.attrs['model_list'] = '%s'%(Att_model_list)
    ds_stack_P75.attrs['reference'] = '%s'%(Att_reference)
    ds_stack_P75.attrs['creation_date'] = '%s'%(Att_creation_date)
    ds_stack_P75.attrs['convention'] = '%s'%(Att_convention)
    ds_stack_P75.attrs['tracking_id'] = '%s'%(Att_tracking_id_M),
    ds_stack_P75.attrs['history'] = '%s'%(Att_history)
    ds_stack_P75.attrs['disclaimer'] = '%s'%(Att_Disclaimer)


    outfn_P75 = '%sMME75_%s_fwi_metrics_monthly_%s.nc'%(outroot_M,irun,iyear)
    outfn_P50 = '%sMME50_%s_fwi_metrics_monthly_%s.nc'%(outroot_M,irun,iyear)
    outfn_P25 = '%sMME25_%s_fwi_metrics_monthly_%s.nc'%(outroot_M,irun,iyear)

    ds_stack_P75.to_netcdf(path=outfn_P75)
    ds_stack_P50.to_netcdf(path=outfn_P50)
    ds_stack_P25.to_netcdf(path=outfn_P25)




#-------------------------------------------------------------------------------
def main():

    inputid = int(sys.argv[1])

    runtypeid = 1

    # historical = 65
    # ssp245 & 585 = 86
    get_Monthly_MME(runtypeid,inputid) # Not working at this moment.
    



#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#__main__
