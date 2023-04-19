# This code is for calculating fwi multimodel ensemble median across all CMIP6 models. 

# Taejin Park
# 2023-02-10
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

#suppress warnings
warnings.filterwarnings('ignore')

root_Y = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/Yearly/'
root_M = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/Monthly/'
outroot_Y = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/Yearly_MME/'
outroot_M = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/Monthly_MME/'

Models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', \
    'CESM2-WACCM', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CNRM-CM6-1', 'CNRM-ESM2-1', \
    'EC-Earth3', 'EC-Earth3-Veg-LR', 'FGOALS-g3', 'GFDL-CM4', 'GFDL-CM4_gr2', \
    'GFDL-ESM4', 'GISS-E2-1-G', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'IITM-ESM', \
    'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'KACE-1-0-G', 'KIOST-ESM', \
    'MIROC6', 'MIROC-ES2L', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', \
    'NESM3', 'NorESM2-LM', 'NorESM2-MM', 'TaiESM1', 'UKESM1-0-LL'] # total 35 models

runtypes = ['historical','ssp126','ssp245','ssp370','ssp585']# total 5 run types



#-------------------------------------------------------------------------------
def get_Yearly_MME(runtypeid,inputid):

    irun = runtypes[runtypeid]

    if irun == 'historical':
        syear = 1950
        eyear = 2014
    else:
        syear = 2015
        eyear = 2100

    years = list(range(syear,eyear+1))
    iyear = years[inputid]

    print('irun = %s'%(irun))
    print('iyear = %d'%(iyear))
    
    fn = '%s*/*_%s_*_%s.nc'%(root_Y,irun,iyear)
    fn_models = glob.glob(fn)
    
    
    print(iyear)
   
    icount = 0
    for imodel in fn_models:
        print('imodel = %s'%(imodel))
        ds = np.squeeze(xr.open_dataset(imodel))
        ds = ds.assign_coords({'modelid':icount})
        ds = ds.expand_dims(dim={"modelid": 1})
        #ds = ds.expand_dims(dim={"year": 1})

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
    ds_stack_P50['FWI_N15'] = ds_stack_P50['FWI_N15'].astype('uint16')
    ds_stack_P50['FWI_N30'] = ds_stack_P50['FWI_N30'].astype('uint16')
    ds_stack_P50['FWI_N45'] = ds_stack_P50['FWI_N45'].astype('uint16')
    ds_stack_P75['FWI_N15'] = ds_stack_P75['FWI_N15'].astype('uint16')
    ds_stack_P75['FWI_N30'] = ds_stack_P75['FWI_N30'].astype('uint16')
    ds_stack_P75['FWI_N45'] = ds_stack_P75['FWI_N45'].astype('uint16')

    outfn_P75 = '%sMME75_%s_fwi_metrics_yearly_%s.nc'%(outroot_Y,irun,iyear)
    outfn_P50 = '%sMME50_%s_fwi_metrics_yearly_%s.nc'%(outroot_Y,irun,iyear)
    outfn_P25 = '%sMME25_%s_fwi_metrics_yearly_%s.nc'%(outroot_Y,irun,iyear)

    ds_stack_P75.to_netcdf(path=outfn_P75)
    ds_stack_P50.to_netcdf(path=outfn_P50)
    ds_stack_P25.to_netcdf(path=outfn_P25)





#-------------------------------------------------------------------------------
def main():

    inputid = int(sys.argv[1])

    runtypeid = 3 # 0, 1, 2, 3, 4 for historical, ssp127,245 370, & 585 
    # historical = 65
    # ssp127,245 370, & 585 = 86
    get_Yearly_MME(runtypeid,inputid) 
    


#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#__main__
