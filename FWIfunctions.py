# This is a collection of functions including functions for the Canadian Fire Weather Index

# Taejin Park
# 2021-10-06
# taejin1392@gmail.com


"""
# The following functions are adopted from 
# https://github.com/buckinha/pyfwi.git
# Main modification is changing point based calculation to 2d array based calculation using numpy array

# FWI Functions:
# FFMC - takes temperature, relative humidity, wind speed, rain, and a previous FFMC value to produce the current FFMC value
# DMC - takes temperature, relative humidity, rainfall, previous DMC value, latitude, and current month to produce the current DMC value
# DC - takes temperature, rain, the previous DC value, latititude, and current month to produce the current DC value
# ISI - takes the wind speed and current FFMC value to produce the current ISI value
# BUI - takes the current DMC and DC values to produce the current BUI value
# FWI - takes the current ISI and BUI values to produce the current FWI value
"""


import numpy as np

def calFFMC(day_tas,day_rh,day_prec,day_wind,ffmc0): 
    """
    Calculates today's Fine Fuel Moisture Code

    PARAMETERS
    ----------
    TEMP (day_tas) is the 12:00 LST temperature in degrees celsius
    RH (day_rh) is the 12:00 LST relative humidity in %
    WIND (day_wind) is the 12:00 LST wind speed in kph
    RAIN (day_prec) is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
    FFMCPrev is the previous day's FFMC    
    """

    mo = (147.2*(101.0 - ffmc0))/(59.5 + ffmc0)   

    rf = day_prec-0.5
    mr = np.full(day_tas.shape,250)
    mr1 = mo + 42.5 * rf * np.exp(-100.0 / (251.0 - mo)) * (1.0-np.exp(-6.93 / rf))
    mr2 = mo + 42.5 * rf * np.exp(-100.0 / (251.0 - mo)) * (1.0-np.exp(-6.93 / rf)) + 0.0015 * pow(mo - 150.0, 2) * pow(rf, .5)
    mr = np.where(mo<=150, mr1, mr2)
    mr = np.where(mr>250, 250, mr)
    
    mo = np.where(day_prec>0.5, mr, mo)

    ed = 0.942 * pow(day_rh, 0.679) + 11.0 * np.exp((day_rh - 100.0) / 10.0) + 0.18 * (21.1 - day_tas) * (1.0 - np.exp(-0.115 * day_rh))
  
    ko = 0.424 * (1.0 - pow(day_rh / 100.0, 1.7)) + 0.0694 * pow(day_wind, 0.5) * (1.0 - pow(day_rh / 100.0, 8))
    kd = ko * 0.581 * np.exp(0.0365 * day_tas)
    m1 = ed + (mo - ed) * pow(10.0,-kd)

    ew = 0.618 * pow(day_rh,0.753) + 10.0 * np.exp((day_rh - 100.0) / 10.0) + 0.18 * (21.1 - day_tas) * (1.0 - np.exp(-0.115 * day_rh))
    k1 = 0.424 * (1.0 - pow((100.0 - day_rh) / 100.0, 1.7)) + 0.0694 * pow(day_wind, 0.5) * (1.0 - pow((100.0 - day_rh) / 100.0, 8))
    kw = k1 * 0.581 * np.exp(0.0365 * day_tas)
    m2 = ew - (ew - mo) * pow(10.0, -kw)
    
    m = mo
    m = np.where(mo>ed, m1, m)
    m = np.where(mo<ed, m2, m)

    FFMC = 59.5 * (250.0 - m) / (147.2 + m)
    return FFMC
# end: calFFMC



def calDMC(day_tas,day_rh,day_prec,lat_var,iday,dmc0):
    """
    Calculates today's Duff Moisture Code   

    PARAMETERS
    ----------
    TEMP is the 12:00 LST temperature in degrees celsius
    RH is the 12:00 LST relative humidity in %
    RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
    DMCPrev is the prevvious day's DMC
    Lat is the latitude in decimal degrees of the location for which calculations are being made
    Month is the month of Year (1..12) for the current day's calculations.
    """
    
    re = 0.92 * day_prec - 1.27
    mo = 20.0 + np.exp(5.6348 - dmc0 / 43.43)
    b1 = 100.0 / (0.5 + 0.3 * dmc0)
    b2 = 14.0 - 1.3 * np.log(dmc0)
    b3 = 6.2 * np.log(dmc0) - 17.2
    b = b3
    b = np.where(dmc0<=33.0, b1, b)
    b = np.where((dmc0>33.0) & (dmc0<=65.0), b2, b)
    mr = mo + 1000.0 * re / (48.77 + b * re)
    pr = 244.72 - 43.43 * np.log(mr - 20.0)
    dmc1 = np.where(pr>0.0, pr, 0)
    dmc0 = np.where(day_prec>1.5,dmc1,dmc0)

    d1 = daylength(iday,lat_var)
    k = 1.894 * (day_tas + 1.1) * (100.0 - day_rh) * d1 * 0.000001
    k = np.where(day_tas>-1.1,k,0.0)
    DMC = dmc0 + 100.0 * k
    return DMC
# end: calDMC



def calDC(day_tas,day_prec,lat_var,imonth,dc0):
    """
    Calculates today's Drought Code Parameters:
    TEMP is the 12:00 LST temperature in degrees celsius
    RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
    DCPrev is the previous day's DC
    LAT is the latitude in decimal degrees of the location for which calculations are being made
    MONTH is the month of Year (1..12) for the current day's calculations.
    """
    rd = 0.83 * day_prec - 1.27
    Qo = 800.0 * np.exp(-dc0 / 400.0)
    Qr = Qo + 3.937 * rd
    Dr = 400.0 * np.log(800.0 / Qr)
    dc1 = np.where(Dr>0.0,Dr,0.0)
    dc0 = np.where(day_prec>2.8,dc1,dc0)

    Lf = DryingFactor(lat_var, imonth-1)
    V = 0.36 * (day_tas+2.8) + Lf
    V = np.where(day_tas>-2.8,V,Lf)
    V = np.where(V<0.0,0.0, V)
    DC = dc0 + 0.5 * V
    return DC
# end: calDC



def calISI(day_wind,FFMC):
    """
    Calculates today's Initial Spread Index

    PARAMETERS
    ----------
    WIND is the 12:00 LST wind speed in kph
    FFMC is the current day's FFMC
    """
    fWIND = np.exp(0.05039 * day_wind)
    m = 147.2 * (101.0 - FFMC) / (59.5 + FFMC)
    fF = 91.9 * np.exp(-0.1386 * m) * (1.0 + pow(m, 5.31) / 49300000.0)
    ISI = 0.208 * fWIND * fF
    return ISI
# end: calISI



def calBUI(DMC,DC):
    """
    Calculates today's Buildup Index

    PARAMETERS
    ----------
    DMC is the current day's Duff Moisture Code
    DC is the current day's Drought Code
    """
    U1 = 0.8 * DMC * DC / (DMC + 0.4 * DC)
    U2 = DMC - (1.0 - 0.8 * DC / (DMC + 0.4 * DC)) * (0.92 + pow(0.0114 * DMC, 1.7))
    U = np.where(DMC <= 0.4 * DC, U1, U2)
    U = np.where(U  < 0.0, 0.0, U)
    BUI = U
    return BUI
# end: calBUI


def calFWI(ISI,BUI):
    """
    Calculates today's Fire Weather Index

    PARAMETERS
    ----------
    ISI is the current day's ISI
    BUI is the current day's BUI
    """

    fD1 = 0.626 * pow(BUI, 0.809) + 2.0
    fD2 = 1000.0 / (25.0 + 108.64 * np.exp(-0.023 * BUI))
    fD = np.where(BUI<=80.0,fD1,fD2)
    B = 0.1 * ISI * fD
    S1 = np.exp(2.72 * pow(0.434 * np.log(B), 0.647))
    S2 = B
    S = np.where(B>1.0,S1,S2)
    FWI = S
    return FWI
# end: calFWI



def DryingFactor(lat, Month):
    LfN = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6]
    LfS = [6.4, 5.0, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8]

    LfN_month = LfN[Month]
    LfS_month = LfS[Month]
    retVal = np.where(lat>0,LfN_month,LfS_month) # revised based on GFWED code
    retVal = np.where(np.logical_and(lat>-15,lat<=15),1.39,retVal)  # Equatorial numbers

    return retVal
# end: DryingFactor


def daylength(dayOfYear, lat):
    """
    Forsythe et al., "A model comparison for daylength as a
    function of latitude and day of year", Ecological Modelling,
    1995.
    Parameters
    ----------
    dayOfYear : int
        The day of the year. 1 corresponds to 1st of January
        and 365 to 31st December (on a non-leap year).
    lat : float
        Latitude of the location in degrees. Positive values
        for north and negative for south.
    Returns
    -------
    d : float
        Daylength in hours.
    """
    latInRad = np.deg2rad(lat)
    declinationOfEarth = 23.45*np.sin(np.deg2rad(360.0*(283.0+dayOfYear)/365.0))
    hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))))
    hourAngle = 2.0*hourAngle/15.0
    ckidx = -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))
    hourAngle = np.where(ckidx<=-1.0,24.0,hourAngle)
    hourAngle = np.where(ckidx>=1.0,0.0,hourAngle)  

    return hourAngle
# end: daylength