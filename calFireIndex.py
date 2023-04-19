# This code is for calculating fire realted indices from the newly generated NASA NEX GDDP2 data 
# which is downscaled from CMIP6 models.

# Taejin Park
# 2021-10-07
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
import os, glob, sys

#suppress warnings
warnings.filterwarnings('ignore')



root = '/nex/datapool/nex-gddp-cmip6/'
outroot = '/nobackupp10/tpark3/Projects/GDDP2/release/GDDP-FWI/'


def FireWeatherIndex(imodel,irun):
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

    outdir_d = '%sDaily/%s/'%(outroot,imodel)
    if not os.path.exists(outdir_d):
        os.makedirs(outdir_d)
    outdir_m = '%sMonthly/%s/'%(outroot,imodel)
    if not os.path.exists(outdir_m):
        os.makedirs(outdir_m)
    outdir_y = '%sYearly/%s/'%(outroot,imodel)
    if not os.path.exists(outdir_y):
        os.makedirs(outdir_y)

    if irun == 'historical':
        syear = 1950
        eyear = 2014
        nyear = eyear-syear+1
    else:
        syear = 2015
        eyear = 2100
        nyear = eyear-syear+1

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
        fn_tas = '%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var1,var1,imodelfile,irun,iyr)
        fn_tas_list = glob.glob(fn_tas)
        tas = xr.open_dataset(fn_tas_list[0])

        var2 = 'hurs'
        fn_hurs = '%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var2,var2,imodelfile,irun,iyr)
        fn_hurs_list = glob.glob(fn_hurs)
        hurs = xr.open_dataset(fn_hurs_list[0])

        var3 = 'sfcWind'
        fn_wind = '%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var3,var3,imodelfile,irun,iyr)
        fn_wind_list = glob.glob(fn_wind)
        wind = xr.open_dataset(fn_wind_list[0])

        var4 = 'rlds'
        fn_rlds = '%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var4,var4,imodelfile,irun,iyr)
        fn_rlds_list = glob.glob(fn_rlds)
        rlds = xr.open_dataset(fn_rlds_list[0])

        var5 = 'rsds'
        fn_rsds = '%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var5,var5,imodelfile,irun,iyr)
        fn_rsds_list = glob.glob(fn_rsds)
        rsds = xr.open_dataset(fn_rsds_list[0])

        var6 = 'pr'
        fn_pr = '%s%s/%s/%s/%s/%s_day_%s_%s*%04i.nc'%(root,imodel,irun,ireal,var6,var6,imodelfile,irun,iyr)
        fn_pr_list = glob.glob(fn_pr)
        pr = xr.open_dataset(fn_pr_list[0])

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
        radlw_var = rlds[var4]
        radsw_var = rsds[var5]
        rad_var = radlw_var+radsw_var
        prec_var = pr[var6]

        outfn_fwi_d = '%s%s_%s_fwi_metrics_daily_%04i.nc'%(outdir_d,imodel,irun,iyr)
        outfn_fwi_m = '%s%s_%s_fwi_metrics_monthly_%04i.nc'%(outdir_m,imodel,irun,iyr)
        outfn_fwi_y = '%s%s_%s_fwi_metrics_yearly_%04i.nc'%(outdir_y,imodel,irun,iyr)

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

        for iday in days:

            idate = datetime.datetime.strptime('{} {}'.format(iday, iyr),'%j %Y')
            imonth = idate.month

            print('                  date = %s'%(idate))
            day_time = time_var[iday-1].squeeze()
            day_tas = tas_var[iday-1,:,:].squeeze()
            day_rh = rh_var[iday-1,:,:].squeeze()
            day_rh = np.where(day_rh>100,100,day_rh)
            day_wind = wind_var[iday-1,:,:].squeeze()
            day_wind = day_wind*3.6 # unit conversion from m/2 to km/h 
            day_wind = np.where(day_wind<0,0,day_wind)
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



        OutXR = xr.Dataset({
        'FFMC': xr.DataArray(
                    data   = FFMC_Year,   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time':time_var,'lat':lat,'lon':lon},
                    attrs  = {'fullname': 'Fine Fuel Moisture Code'}
                    ),
        'DMC': xr.DataArray(
                    data   = DMC_Year,   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time':time_var,'lat':lat,'lon':lon},
                    attrs  = {'fullname': 'Duff Moisture Code'}
                    ),
        'DC': xr.DataArray(
                    data   = DC_Year,   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time':time_var,'lat':lat,'lon':lon},
                    attrs  = {'fullname': 'Drought Code'}
                    ),
        'ISI': xr.DataArray(
                    data   = ISI_Year,   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time':time_var,'lat':lat,'lon':lon},
                    attrs  = {'fullname': 'Initial Spread Index'}
                    ),
        'BUI': xr.DataArray(
                    data   = BUI_Year,   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time':time_var,'lat':lat,'lon':lon},
                    attrs  = {'fullname': 'Buildup Index'}
                    ),
        'FWI': xr.DataArray(
                    data   = FWI_Year,   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time':time_var,'lat':lat,'lon':lon},
                    attrs  = {'fullname': 'Fire Weather Index'}
                    )
                },
            attrs = {'deascription': 'These are ouputs of Canadian Fire Weather Index System'}
        )

        if iyr == eyear:
            OutXR.to_netcdf(path=outfn_fwi_d)

        #del OutXR

        

        threshold = 15
        FWI_M_N15 = OutXR['FWI'].where(OutXR['FWI'] > threshold).groupby('time.month').count(dim='time').astype('uint16')
        FWI_N15 = OutXR['FWI'].where(OutXR['FWI'] > threshold).groupby('time.year').count(dim='time').astype('uint16')
        FWI_N15 = FWI_N15.drop_vars('year').squeeze()
        threshold = 30
        FWI_M_N30 = OutXR['FWI'].where(OutXR['FWI'] > threshold).groupby('time.month').count(dim='time').astype('uint16')
        FWI_N30 = OutXR['FWI'].where(OutXR['FWI'] > threshold).groupby('time.year').count(dim='time').astype('uint16') 
        FWI_N30 = FWI_N30.drop_vars('year').squeeze()
        threshold = 45
        FWI_M_N45 = OutXR['FWI'].where(OutXR['FWI'] > threshold).groupby('time.month').count(dim='time').astype('uint16')
        FWI_N45 = OutXR['FWI'].where(OutXR['FWI'] > threshold).groupby('time.year').count(dim='time').astype('uint16')
        FWI_N45 = FWI_N45.drop_vars('year').squeeze()
        
        OutXR_M = OutXR.resample(time="M").mean()
        time_month = OutXR_M['time'] 
        FWI_M_N15 = FWI_M_N15.drop('month').rename({'month':'time'}).assign_coords(time=time_month) 
        FWI_M_N30 = FWI_M_N30.drop('month').rename({'month':'time'}).assign_coords(time=time_month)
        FWI_M_N45 = FWI_M_N45.drop('month').rename({'month':'time'}).assign_coords(time=time_month)
        OutXR_M['FWI_N15'] = FWI_M_N15
        OutXR_M['FWI_N30'] = FWI_M_N30
        OutXR_M['FWI_N45'] = FWI_M_N45
        OutXR_M.FFMC.attrs['fullname'] = 'monthly mean FFMC (Fine Fuel Moisture Code)'
        OutXR_M.DMC.attrs['fullname'] = 'monthly mean DMC (Duff Moisture Code)'
        OutXR_M.DC.attrs['fullname'] = 'monthly mean DC (Drought Code)'
        OutXR_M.ISI.attrs['fullname'] = 'monthly mean ISI (Initial Spread Index)'
        OutXR_M.BUI.attrs['fullname'] = 'monthly mean BUI (Buildup Index)'
        OutXR_M.FWI.attrs['fullname'] = 'monthly mean FWI (Fire Weather Index)'
        OutXR_M.FWI_N15.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 15'
        OutXR_M.FWI_N30.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 30'
        OutXR_M.FWI_N45.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 45'
        OutXR_M.attrs['deascription'] = 'These are ouputs of Canadian Fire Weather Index System'
        print(OutXR_M)
        OutXR_M.to_netcdf(path=outfn_fwi_m)

        FWI_P95 = OutXR['FWI'].groupby('time.year').quantile(0.95, dim='time').astype('float32')
        FWI_P95 = FWI_P95.drop_vars('quantile')
        FWI_P95 = FWI_P95.drop_vars('year').squeeze()
        FWI_P75 = OutXR['FWI'].groupby('time.year').quantile(0.75, dim='time').astype('float32')
        FWI_P75 = FWI_P75.drop_vars('quantile')
        FWI_P75 = FWI_P75.drop_vars('year').squeeze()
        FWI_P50 = OutXR['FWI'].groupby('time.year').quantile(0.50, dim='time').astype('float32')
        FWI_P50 = FWI_P50.drop_vars('quantile')
        FWI_P50 = FWI_P50.drop_vars('year').squeeze()
        FWI_P25 = OutXR['FWI'].groupby('time.year').quantile(0.25, dim='time').astype('float32')
        FWI_P25 = FWI_P25.drop_vars('quantile')
        FWI_P25 = FWI_P25.drop_vars('year').squeeze()

        OutXR_Y = OutXR.resample(time="Y").mean().squeeze().drop_vars('time')


        SummaryXR = xr.Dataset({
        'FFMC': xr.DataArray(
                    data   = OutXR_Y.FFMC,   # enter data here                
                    attrs  = {'fullname': 'annual mean FFMC (Fine Fuel Moisture Code)'}
                    ),
        'DMC': xr.DataArray(
                    data   = OutXR_Y.DMC,   # enter data here                
                    attrs  = {'fullname': 'annual mean DMC (Duff Moisture Code)'}
                    ),
        'DC': xr.DataArray(
                    data   = OutXR_Y.DC,   # enter data here                
                    attrs  = {'fullname': 'annual mean DC (Drought Code)'}
                    ),
        'ISI': xr.DataArray(
                    data   = OutXR_Y.ISI,   # enter data here                
                    attrs  = {'fullname': 'annual mean ISI (Initial Spread Index'}
                    ),
        'BUI': xr.DataArray(
                    data   = OutXR_Y.BUI,   # enter data here                
                    attrs  = {'fullname': 'annual mean BUI (Buildup Index)'}
                    ),
        'FWI': xr.DataArray(
                    data   = OutXR_Y.FWI,   # enter data here                
                    attrs  = {'fullname': 'annual mean FWI (Fire Weather Index)'}
                    ),
        'FWI_N15': xr.DataArray(
                    data   = FWI_N15.astype('uint16'),   # enter data here                
                    attrs  = {'fullname': 'number of days FWI (Fire Weather Index) > 15'}
                    ),
        'FWI_N30': xr.DataArray(
                    data   = FWI_N30.astype('uint16'),   # enter data here                
                    attrs  = {'fullname': 'number of days FWI (Fire Weather Index) > 30'}
                    ),
        'FWI_N45': xr.DataArray(
                    data   = FWI_N45.astype('uint16'),   # enter data here                
                    attrs  = {'fullname': 'number of days FWI (Fire Weather Index) > 45'}
                    ),
        'FWI_P25': xr.DataArray(
                    data   = FWI_P25,   # enter data here                
                    attrs  = {'fullname': '25th percentile FWI (Fire Weather Index)'}
                    ),
        'FWI_P50': xr.DataArray(
                    data   = FWI_P50,   # enter data here                
                    attrs  = {'fullname': '50th percentile FWI (Fire Weather Index)'}
                    ),
        'FWI_P75': xr.DataArray(
                    data   = FWI_P75,   # enter data here                
                    attrs  = {'fullname': '75th percentile FWI (Fire Weather Index)'}
                    ),
        'FWI_P95': xr.DataArray(
                    data   = FWI_P95,   # enter data here                
                    attrs  = {'fullname': '95th percentile FWI (Fire Weather Index)'}
                    )
                },
            attrs = {'deascription': 'These are ouputs of Canadian Fire Weather Index System'}
        )

        print(SummaryXR)
        SummaryXR.to_netcdf(path=outfn_fwi_y)


        end = timer()
        print('           elapsed time = %s'%(timedelta(seconds=end-start)))

    # fi: years



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
    
    histfirstid = False
    # True: inputid can range from 0 to 34 covering all models for historical run
    # False: inputid can range from 0 to 35*4-1=159  covering all models and two simulation scenarios
    # Historical: 65 years * 10min = ~11hours (14 hours recommended for each)
    # Scenarios: 86 years * 10min = ~14.3 (17 hours recommended for each) 

    if histfirstid:
        modelid = inputid
        imodel = Models[modelid]
        irun = runtypes[0]
    else:
        modelid,runid=divmod(inputid,4)
        imodel = Models[modelid]
        irun = runtypes[runid+1]
   
    print('imodel = %s, irun = %s.... '%(imodel,irun))

    FireWeatherIndex(imodel,irun)
    



#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#__main__

