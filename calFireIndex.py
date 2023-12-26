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
    
    outdir_d_baseline = '%sDaily_Baseline/%s/'%(outroot,imodel)
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
        Att_variant_label = ireal
        Att_scenario = irun
        Att_tracking_id_D = str(uuid.uuid4())
        Att_tracking_id_M = str(uuid.uuid4())
        Att_tracking_id_Y = str(uuid.uuid4())
        Att_cmip6_source_id = pr.cmip6_source_id
        Att_cmip6_institution_id = pr.cmip6_institution_id
        Att_cmip6_license = pr.cmip6_license
        Att_resolution_id = pr.resolution_id
        Att_creation_date = datetime.datetime.now(pytz.timezone('UTC')) 
        Att_history = '%s %s'%(Att_creation_date,'Version 1.0 production')
        Att_title = '%s, %s, %s, %s'%(Att_cmip6_source_id,Att_variant_label,Att_scenario,'global fire weather index data from downscaled CMIP6 projection data')
        Att_reference ='NEX-GDDP-FWI description paper: Park, T. et al., Under Review. NEX-GDDP-FWI: Downscaled 21st century global fire weather projections. Nature Scientific Data // NEX-GDDP description paper: Thrasher, B. et al., 2022. NASA global daily downscaled projections, CMIP6. Nature Scientific data (https://doi.org/10.1038/s41597-022-01393-4) // Canadian Fire Weather Index System: Van Wagner, C.E., 1987. Development and structure of the Canadian forest fire weather index system, Canadian Forest Service (http://cfs.nrcan.gc.ca/pubwarehouse/pdfs/19927.pdf)'


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
            attrs = {'activity': '%s'%(Att_Activity),
                     'institution': '%s'%(Att_Institution),
                     'contact': '%s'%(Att_Contact),                              
                     'variant_label': '%s'%(Att_variant_label),
                     'scenario': '%s'%(Att_scenario),
                     'resolution_id': '%s'%(Att_resolution_id),
                     'frequency': '%s'%(Att_Frequency_D),            
                     'version': '%s'%(Att_Version),
                     'title': '%s'%(Att_title),
                     'reference': '%s'%(Att_reference),
                     'creation_date': '%s'%(Att_creation_date),
                     'convention': '%s'%(Att_convention),
                     'history': '%s'%(Att_history),
                     'tracking_id': '%s'%(Att_tracking_id_D),
                     'disclaimer': '%s'%(Att_Disclaimer),
                     'cmip6_source_id': '%s'%(Att_cmip6_source_id),
                     'cmip6_institution_id': '%s'%(Att_cmip6_institution_id),
                     'cmip6_license': '%s'%(Att_cmip6_license),
                     }
        )

        if iyr == eyear:
            OutXR.to_netcdf(path=outfn_fwi_d)

        #del OutXR

        # Import baseline thresholds
        fn_fwi_baseline_min = '%s%s_%s_fwi_metrics_baseline_min_%d_%d.nc'%(outdir_d_baseline,imodel,'historical',RefSyear,RefEyear)
        fn_fwi_baseline_max = '%s%s_%s_fwi_metrics_baseline_max_%d_%d.nc'%(outdir_d_baseline,imodel,'historical',RefSyear,RefEyear)
        fn_fwi_baseline_P95 = '%s%s_%s_fwi_metrics_baseline_P95_%d_%d.nc'%(outdir_d_baseline,imodel,'historical',RefSyear,RefEyear)
        baseline_min = xr.open_dataset(fn_fwi_baseline_min)
        baseline_max = xr.open_dataset(fn_fwi_baseline_max)
        baseline_dif = baseline_max-baseline_min
        baseline_P95 = xr.open_dataset(fn_fwi_baseline_P95) 

        # Baseline threshold based on the midpoint of reference period (1950-1979)
        FWI_M_Nmid = OutXR['FWI'].where((OutXR['FWI']-baseline_min)/(baseline_dif) > 0.5).groupby('time.month').count(dim='time').astype('uint16')
        FWI_M_Nmid = FWI_M_Nmid.FWI
        FWI_Nmid = OutXR['FWI'].where((OutXR['FWI']-baseline_min)/(baseline_dif) > 0.5).groupby('time.year').count(dim='time').astype('uint16')
        FWI_Nmid = FWI_Nmid.drop_vars('year').squeeze()
        FWI_Nmid = FWI_Nmid.FWI
        
        # Baseline threshold based on the 95 percentile of reference period (1950-1979)
        FWI_M_NP95 = OutXR['FWI'].where(OutXR['FWI'] > baseline_P95).groupby('time.month').count(dim='time').astype('uint16')
        FWI_M_NP95 = FWI_M_NP95.FWI
        FWI_NP95 = OutXR['FWI'].where(OutXR['FWI'] > baseline_P95).groupby('time.year').count(dim='time').astype('uint16')
        FWI_NP95 = FWI_NP95.drop_vars('year').squeeze()
        FWI_NP95 = FWI_NP95.FWI

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
        FWI_M_Nmid = FWI_M_Nmid.drop('month').rename({'month':'time'}).assign_coords(time=time_month) 
        FWI_M_NP95 = FWI_M_NP95.drop('month').rename({'month':'time'}).assign_coords(time=time_month) 
        FWI_M_N15 = FWI_M_N15.drop('month').rename({'month':'time'}).assign_coords(time=time_month) 
        FWI_M_N30 = FWI_M_N30.drop('month').rename({'month':'time'}).assign_coords(time=time_month)
        FWI_M_N45 = FWI_M_N45.drop('month').rename({'month':'time'}).assign_coords(time=time_month)
        OutXR_M['FWI_Nmid'] = FWI_M_Nmid
        OutXR_M['FWI_NP95'] = FWI_M_NP95
        OutXR_M['FWI_N15'] = FWI_M_N15
        OutXR_M['FWI_N30'] = FWI_M_N30
        OutXR_M['FWI_N45'] = FWI_M_N45
        OutXR_M.FFMC.attrs['fullname'] = 'monthly mean FFMC (Fine Fuel Moisture Code)'
        OutXR_M.DMC.attrs['fullname'] = 'monthly mean DMC (Duff Moisture Code)'
        OutXR_M.DC.attrs['fullname'] = 'monthly mean DC (Drought Code)'
        OutXR_M.ISI.attrs['fullname'] = 'monthly mean ISI (Initial Spread Index)'
        OutXR_M.BUI.attrs['fullname'] = 'monthly mean BUI (Buildup Index)'
        OutXR_M.FWI.attrs['fullname'] = 'monthly mean FWI (Fire Weather Index)'
        OutXR_M.FWI_Nmid.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > midpoint of reference period'
        OutXR_M.FWI_NP95.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 95 percentile of reference period'
        OutXR_M.FWI_N15.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 15'
        OutXR_M.FWI_N30.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 30'
        OutXR_M.FWI_N45.attrs['fullname'] = 'number of days FWI (Fire Weather Index) > 45'
        OutXR_M.attrs['activity'] = '%s'%(Att_Activity)
        OutXR_M.attrs['institution'] = '%s'%(Att_Institution)
        OutXR_M.attrs['contact'] = '%s'%(Att_Contact)
        OutXR_M.attrs['variant_label'] = '%s'%(Att_variant_label)
        OutXR_M.attrs['scenario'] = '%s'%(Att_scenario)
        OutXR_M.attrs['resolution_id'] = '%s'%(Att_resolution_id)
        OutXR_M.attrs['frequency'] = '%s'%(Att_Frequency_M)
        OutXR_M.attrs['version'] = '%s'%(Att_Version)
        OutXR_M.attrs['title'] = '%s'%(Att_title)
        OutXR_M.attrs['reference'] = '%s'%(Att_reference)
        OutXR_M.attrs['creation_date'] = '%s'%(Att_creation_date)
        OutXR_M.attrs['convention'] = '%s'%(Att_convention)
        OutXR_M.attrs['tracking_id'] = '%s'%(Att_tracking_id_M),
        OutXR_M.attrs['history'] = '%s'%(Att_history)
        OutXR_M.attrs['disclaimer'] = '%s'%(Att_Disclaimer)
        OutXR_M.attrs['cmip6_source_id'] = '%s'%(Att_cmip6_source_id)
        OutXR_M.attrs['cmip6_institution_id'] = '%s'%(Att_cmip6_institution_id)
        OutXR_M.attrs['cmip6_license'] = '%s'%(Att_cmip6_license)
        print(OutXR_M)
        OutXR_M.to_netcdf(path=outfn_fwi_m, mode='w')

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
        'FWI_Nmid': xr.DataArray(
                    data   = FWI_Nmid.astype('uint16'),   # enter data here                
                    attrs  = {'fullname': 'number of days FWI (Fire Weather Index) > midpoint of reference period'}
                    ),
        'FWI_NP95': xr.DataArray(
                    data   = FWI_NP95.astype('uint16'),   # enter data here                
                    attrs  = {'fullname': 'number of days FWI (Fire Weather Index) > 95 percentile of reference period'}
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
            attrs = {'activity': '%s'%(Att_Activity),
                     'institution': '%s'%(Att_Institution),
                     'contact': '%s'%(Att_Contact),                              
                     'variant_label': '%s'%(Att_variant_label),
                     'scenario': '%s'%(Att_scenario),
                     'resolution_id': '%s'%(Att_resolution_id),
                     'frequency': '%s'%(Att_Frequency_Y),            
                     'version': '%s'%(Att_Version),
                     'title': '%s'%(Att_title),
                     'reference': '%s'%(Att_reference),
                     'creation_date': '%s'%(Att_creation_date),
                     'convention': '%s'%(Att_convention),
                     'tracking_id': '%s'%(Att_tracking_id_Y),
                     'history': '%s'%(Att_history),
                     'disclaimer': '%s'%(Att_Disclaimer),
                     'cmip6_source_id': '%s'%(Att_cmip6_source_id),
                     'cmip6_institution_id': '%s'%(Att_cmip6_institution_id),
                     'cmip6_license': '%s'%(Att_cmip6_license),
                     }
        )

        print(SummaryXR)
        SummaryXR.to_netcdf(path=outfn_fwi_y, mode='w')


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
    # Running order: historical run for baseline => historical run => and other future scenarios
    # as calFWI needs to run historical run first, 
    # we need to first run historical only then run other scenarios parallelly. 
    

    inputid = int(sys.argv[1])
    
    histfirstid = False
    # True: inputid can range from 0 to 34 covering all models for historical run
    # False: inputid can range from 0 to 36*4-1=143  covering all models and two simulation scenarios
    # Historical: 65 years  = ~10 hours recommended for each
    # Scenarios: 86 years * 10min = ~14 hours recommended for each 

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

