# This code is for calculating fire realted indices from the newly generated NEX-GDDP-CMIP6 data 
# which is downscaled from CMIP6 models.

# Taejin Park
# 2023-10-07
# taejin1392@gmail.com


import xarray as xr
#import cartopy
#import cartopy.crs as ccrs
#import matplotlib.pyplot as plt
import numpy as np
import warnings
import datetime
from FWIfunctions import *
from timeit import default_timer as timer
from datetime import timedelta
import pytz
import os, glob, sys
import uuid

#suppress warnings
warnings.filterwarnings('ignore')

#root = '/Users/taejinpark/Projects/Data/GDDP2/CMIP6_GDDP/'
#outroot = '/Users/taejinpark/Projects/Data/GDDP2/FWI_Rev1/'

root = '/nex/datapool/nex-gddp-cmip6/'
outroot = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/'

#vegmask_fn= '/nobackupp10/tpark3/Projects/GDDP2/GDDP2_heat_fire/MCD12C1_VegMask/MCD12C1.A2001001.006.2018053185512_25km_mask.nc' 
#vegmask = xr.open_dataset(vegmask_fn)
#vegmask_array = vegmask['VegMask'] # mask for non-vegetated areas (MCD12C1 IGBP==0,11,13,>15)


RefSyear = 1951
RefEyear = 1980

def FireWeatherIndexBaseline(imodel,irun):

    if imodel=='GFDL-CM4_gr2':
        imodelfile = 'GFDL-CM4'
    else:
        imodelfile = imodel    

    realization = []
    realdir = '%s%s/%s/*/'%(root,imodel,irun)
    realdirpath = glob.glob(realdir)
    print(realdirpath)
    realization.append(os.path.basename(os.path.dirname(realdirpath[0])))
    ireal = realization[0]

    vartype = ['hurs','rlds','rsds','sfcWind','tas','tasmax','tasmin']
    runtypes = ['historical','ssp245','ssp585']

    print('imodel = %s, irun = %s, irealization = %s '%(imodel,irun,ireal))

    outdir_d = '%sDaily_Baseline/%s/'%(outroot,imodel)
    if not os.path.exists(outdir_d):
        os.makedirs(outdir_d)

    syear = 1950
    eyear = RefEyear

    years = list(range(syear,eyear+1))

    for iyr in years:
        start = timer()

        print('       iyr = %d'%(iyr))
        print('           reading GDDP2 data...')

        switchid = 0

        if imodel == 'KIOST-ESM' and irun == 'ssp126' and iyr==2023:
            iyr = 2022
            switchid = 1


        #var1 = 'tasmax'
        var1 = 'tas'
        if glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.2.nc'%(root,imodel,irun,ireal,var1,var1,imodelfile,irun,iyr)):
            fn_tas_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.2.nc'%(root,imodel,irun,ireal,var1,var1,imodelfile,irun,iyr))        
            tas = xr.open_dataset(fn_tas_list[0])
            print(fn_tas_list[0])
        elif glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.1.nc'%(root,imodel,irun,ireal,var1,var1,imodelfile,irun,iyr)):
            fn_tas_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.1.nc'%(root,imodel,irun,ireal,var1,var1,imodelfile,irun,iyr))
            tas = xr.open_dataset(fn_tas_list[0])
            print(fn_tas_list[0])
        else:
            fn_tas_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var1,var1,imodelfile,irun,iyr))
            tas = xr.open_dataset(fn_tas_list[0])
            print(fn_tas_list[0])

        var2 = 'hurs'
        if glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.2.nc'%(root,imodel,irun,ireal,var2,var2,imodelfile,irun,iyr)):
            fn_hurs_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.2.nc'%(root,imodel,irun,ireal,var2,var2,imodelfile,irun,iyr))
            hurs = xr.open_dataset(fn_hurs_list[0])
            print(fn_hurs_list[0])
        elif glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.1.nc'%(root,imodel,irun,ireal,var2,var2,imodelfile,irun,iyr)):
            fn_hurs_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.1.nc'%(root,imodel,irun,ireal,var2,var2,imodelfile,irun,iyr))
            hurs = xr.open_dataset(fn_hurs_list[0])
            print(fn_hurs_list[0])
        else:
            fn_hurs_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var2,var2,imodelfile,irun,iyr))
            hurs = xr.open_dataset(fn_hurs_list[0])
            print(fn_hurs_list[0])

        var3 = 'sfcWind'
        if glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.2.nc'%(root,imodel,irun,ireal,var3,var3,imodelfile,irun,iyr)):
            fn_wind_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.2.nc'%(root,imodel,irun,ireal,var3,var3,imodelfile,irun,iyr))
            wind = xr.open_dataset(fn_wind_list[0])
            print(fn_wind_list[0])
        elif glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.1.nc'%(root,imodel,irun,ireal,var3,var3,imodelfile,irun,iyr)):
            fn_wind_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.1.nc'%(root,imodel,irun,ireal,var3,var3,imodelfile,irun,iyr))
            wind = xr.open_dataset(fn_wind_list[0])
            print(fn_wind_list[0])
        else:
            fn_wind_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var3,var3,imodelfile,irun,iyr))
            wind = xr.open_dataset(fn_wind_list[0])
            print(fn_wind_list[0])

        var4 = 'pr'
        if glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.2.nc'%(root,imodel,irun,ireal,var4,var4,imodelfile,irun,iyr)):
            fn_pr_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.2.nc'%(root,imodel,irun,ireal,var4,var4,imodelfile,irun,iyr))
            pr = xr.open_dataset(fn_pr_list[0])
            print(fn_pr_list[0])
        elif glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.1.nc'%(root,imodel,irun,ireal,var4,var4,imodelfile,irun,iyr)):
            fn_pr_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i_v1.1.nc'%(root,imodel,irun,ireal,var4,var4,imodelfile,irun,iyr))
            pr = xr.open_dataset(fn_pr_list[0])
            print(fn_pr_list[0])
        else:
            fn_pr_list = glob.glob('%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var4,var4,imodelfile,irun,iyr))
            pr = xr.open_dataset(fn_pr_list[0])
            print(fn_pr_list[0])

            

        if switchid == 1:
            iyr = 2023

        lat = tas['lat']
        lon = tas['lon']
        lat_var =  np.tile(lat,(len(lon),1))
        lat_var = np.swapaxes(lat_var,1,0)
        lon_var =  np.tile(lon,(len(lat),1))

        time_var = tas['time']
        tas_var = tas[var1]-273.15
        rh_var = hurs[var2]
        wind_var = wind[var3]
        prec_var = pr[var4]

        outfn_fwi_d = '%s%s_%s_fwi_metrics_daily_%04i_baseline.nc'%(outdir_d,imodel,irun,iyr)

        datasize = tas_var.shape

        print('           calculate daily FWI metrics...')

        Ndays = len(time_var)
        days = list(range(1,Ndays+1))

        FFMC_Year = np.zeros(datasize,dtype='float32')
        DMC_Year = np.zeros(datasize,dtype='float32')
        DC_Year = np.zeros(datasize,dtype='float32')
        ISI_Year = np.zeros(datasize,dtype='float32')
        BUI_Year = np.zeros(datasize,dtype='float32')
        FWI_Year =np.zeros(datasize,dtype='float32')


        # Daily fire weather index metrics calculation 
        for iday in days:

            idate = datetime.datetime.strptime('{} {}'.format(iday, iyr),'%j %Y')
            imonth = idate.month

            print('                  date = %s'%(idate))
            day_time = time_var[iday-1].squeeze()
            day_tas = tas_var[iday-1,:,:].squeeze()
            day_rh = rh_var[iday-1,:,:].squeeze()
            day_rh = np.where(day_rh>100,100,day_rh) # Force relative humidity in the invalid range (i.e., relative humidity>100) to 100
            day_wind = wind_var[iday-1,:,:].squeeze()
            day_wind = day_wind*3.6 # unit conversion from m/2 to km/h 
            day_wind = np.where(day_wind<0,0,day_wind) # Force wind speed in the invalid range (i.e., wind speed < 0) to 0
            day_prec = prec_var[iday-1,:,:].squeeze()
            day_prec = day_prec*86400 # convert km/m2/s to mm/day
            day_nanmask = np.isnan(day_tas)

            if irun=='historical':
                if (iday==1) & (iyr==syear):
                    ffmc0_v = 85.0
                    dmc0_v = 6.0
                    dc0_v = 15.0

                    ffmc0 = np.full(day_tas.shape,ffmc0_v).astype('float32')
                    dmc0 = np.full(day_tas.shape,dmc0_v).astype('float32')
                    dc0 = np.full(day_tas.shape,dc0_v).astype('float32')
                else:
                    ffmc0 = FFMC
                    dmc0 = DMC
                    dc0 = DC
            else:
                if (iday==1) & (iyr==syear):
                    histfn = '%s%s_%s_fwi_metrics_daily_%04i.nc'%(outdir_d,imodel,runtypes[0],syear-1)
                    histfwi = xr.open_dataset(histfn)
                    ffmc0 = histfwi['FFMC'].isel(time=-1)
                    dmc0 = histfwi['DMC'].isel(time=-1)
                    dc0 = histfwi['DC'].isel(time=-1)
                else:
                    ffmc0 = FFMC
                    dmc0 = DMC
                    dc0 = DC


            FFMC = calFFMC(day_tas,day_rh,day_prec,day_wind,ffmc0).astype('float32')
            DMC = calDMC(day_tas,day_rh,day_prec,lat_var,iday,dmc0).astype('float32')
            DC = calDC(day_tas,day_prec,lat_var,imonth,dc0).astype('float32')
            ISI = calISI(day_wind,FFMC).astype('float32')
            BUI = calBUI(DMC,DC).astype('float32')
            FWI = calFWI(ISI,BUI).astype('float32')


            FFMC = np.where(day_nanmask==False,FFMC,np.nan)
            DMC = np.where(day_nanmask==False,DMC,np.nan)
            DC = np.where(day_nanmask==False,DC,np.nan)
            ISI = np.where(day_nanmask==False,ISI,np.nan)
            BUI = np.where(day_nanmask==False,BUI,np.nan)
            FWI = np.where(day_nanmask==False,FWI,np.nan)
            #FWI = np.where(vegmask_array==1,FWI,np.nan)
              
            FFMC_Year[iday-1,:,:] = FFMC
            DMC_Year[iday-1,:,:] = DMC
            DC_Year[iday-1,:,:] = DC
            ISI_Year[iday-1,:,:] = ISI
            BUI_Year[iday-1,:,:] = BUI
            FWI_Year[iday-1,:,:] = FWI

         
        # fi: days
        #print('           write fwi metrics...')

        # Creating xarray for the daily fire weather index metrics  
        OutXR = xr.Dataset({
        'FWI': xr.DataArray(
                    data   = FWI_Year,   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time':time_var,'lat':lat,'lon':lon},
                    attrs  = {'fullname': 'Fire Weather Index'}
                    )
                }
        )

        OutXR.to_netcdf(path=outfn_fwi_d, mode='w')

        #del OutXR


        end = timer()
        print('           elapsed time = %s'%(timedelta(seconds=end-start)))

    # fi: years


    start = timer()

    icount = 0
    for iyear in range(RefSyear,RefEyear+1):

        fn_baseline = '%s%s_%s_fwi_metrics_daily_%04i_baseline.nc'%(outdir_d,imodel,irun,iyear)
        baselinedata = xr.open_dataset(fn_baseline) 
        #data = data.assign_coords({'year':iyear})
        #data = data.expand_dims(dim={"year": 1})
        #ds = ds.expand_dims(dim={"year": 1})

        if icount==0:
            ds_stack = baselinedata
        else:
            ds_stack = xr.combine_by_coords([ds_stack, baselinedata])

        icount += 1

    
    outfn_fwi_baseline_min = '%s%s_%s_fwi_metrics_baseline_min_%d_%d.nc'%(outdir_d,imodel,irun,RefSyear,RefEyear)
    outfn_fwi_baseline_max = '%s%s_%s_fwi_metrics_baseline_max_%d_%d.nc'%(outdir_d,imodel,irun,RefSyear,RefEyear)
    outfn_fwi_baseline_P95 = '%s%s_%s_fwi_metrics_baseline_P95_%d_%d.nc'%(outdir_d,imodel,irun,RefSyear,RefEyear)
    
    baseline_min = ds_stack['FWI'].min(dim='time').astype('float32')
    baseline_max = ds_stack['FWI'].max(dim='time').astype('float32')
    baseline_P95 = ds_stack['FWI'].quantile(0.95, dim='time').astype('float32')
    baseline_P95 = baseline_P95.drop_vars('quantile')

    baseline_min.to_netcdf(path=outfn_fwi_baseline_min, mode='w')
    baseline_max.to_netcdf(path=outfn_fwi_baseline_max, mode='w')
    baseline_P95.to_netcdf(path=outfn_fwi_baseline_P95, mode='w')

    for delfilename in glob.glob('%s%s_%s_fwi_metrics_daily_*_baseline.nc'%(outdir_d,imodel,irun)):
        os.remove(delfilename)


    end = timer()
    print('           baseline data creation elapsed time = %s'%(timedelta(seconds=end-start)))
    
    


#-------------------------------------------------------------------------------
def main():
    Models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', \
        'CESM2-WACCM', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CNRM-CM6-1', 'CNRM-ESM2-1', \
        'EC-Earth3', 'EC-Earth3-Veg-LR', 'FGOALS-g3', 'GFDL-CM4', 'GFDL-CM4_gr2', \
        'GFDL-ESM4', 'GISS-E2-1-G', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'IITM-ESM', \
        'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'KACE-1-0-G', 'KIOST-ESM', \
        'MIROC6', 'MIROC-ES2L', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', \
        'NESM3', 'NorESM2-LM', 'NorESM2-MM', 'TaiESM1', 'UKESM1-0-LL'] # total 35 models


    runtypes = ['historical','ssp126','ssp245','ssp370','ssp585'] # total 5 run types
    # as calFWI needs to run historical run first, 
    # we need to first run historical only then run other scenarios parallelly. 
    

    inputid = int(sys.argv[1])


    modelid = inputid
    imodel = Models[modelid]
    irun = runtypes[0]
    
    # Historical: 30 (baseline) years * 10min = 5 hours (8 hours recommended for each)
    
    print('imodel = %s, irun = %s.... '%(imodel,irun))

    FireWeatherIndexBaseline(imodel,irun)
    



#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#__main__

