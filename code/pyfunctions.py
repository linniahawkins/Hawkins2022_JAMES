#!/usr/bin/env python2.7
 
import sys
import os, glob
import csv 
import numpy as np
import xarray as xr
import calendar
import matplotlib
import matplotlib.pyplot as plt 
from datetime import datetime, date, timedelta
import pandas as pd
import scipy
from scipy.optimize import curve_fit

#####################################
def round_time(dt=None, round_to=60*30):
   if dt == None: 
       dt = datetime.datetime.now()
   seconds = (dt - dt.min).seconds
   rounding = (seconds+round_to/2) // round_to * round_to
   return dt + timedelta(0,rounding-seconds,-dt.microsecond)

def SPA_readin_ecosystemfluxes(in_dir,start,end,freq):
	'''readin ecosystemfluxes data from SPA model simulations
	in_dir = directory of file full path to file
	start = start day; 
	end=end day
	freq=frequency of data: '30min','H'
	'''
	trng_hh = pd.date_range(start, end, freq=freq)
	filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
	df = pd.read_csv(filename,header=0)
	df.index = trng_hh
	df['mon'] = df.index.month
	return df

def find_quadrants(month):

	return month

def readin_leaflayers_Gs(in_dir,startdate,shour,ehour):
	filename = os.path.join(in_dir+'layer_sun_01.csv')
	df1 = pd.read_csv(filename,header=0)
	# get date time                                     
	dates = []
	for ii in range(len(df1)):   
		delta = timedelta(df1['Time (days)'][ii]) 
		dt = startdate + delta
		rdate = round_time(dt,round_to=60*30)     
		dates.append(rdate)
	df1.index = dates
	df1 = df1.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_sun_02.csv')
	df2 = pd.read_csv(filename,header=0)
	df2.index = dates
	df2 = df2.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_sun_03.csv')
	df3 = pd.read_csv(filename,header=0)
	df3.index = dates
	df3 = df3.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_sun_04.csv')
	df4 = pd.read_csv(filename,header=0)
	df4.index = dates
	df4 = df4.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_sun_05.csv')
	df5 = pd.read_csv(filename,header=0)
	df5.index = dates
	df5 = df5.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_sun_06.csv')
	df6 = pd.read_csv(filename,header=0)
	df6.index = dates
	df6 = df6.between_time(shour,ehour)

	# readin shade layers
	filename = os.path.join(in_dir+'layer_shade_01.csv')
	df7 = pd.read_csv(filename,header=0)
	# get date time                                     
	dates_shade = []
	for ii in range(len(df7)):   
		delta = timedelta(df7['Time (days)'][ii]) 
		dt = startdate + delta
		rdate = round_time(dt,round_to=60*30)     
		dates_shade.append(rdate)
	df7.index = dates_shade
	df7 = df7.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_shade_02.csv')
	df8 = pd.read_csv(filename,header=0)
	df8.index = dates_shade
	df8 = df8.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_shade_03.csv')
	df9 = pd.read_csv(filename,header=0)
	df9.index = dates_shade
	df9 = df9.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_shade_04.csv')
	df10 = pd.read_csv(filename,header=0)
	df10.index = dates_shade
	df10 = df10.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_shade_05.csv')
	df11 = pd.read_csv(filename,header=0)
	df11.index = dates_shade
	df11 = df11.between_time(shour,ehour)
	filename = os.path.join(in_dir+'layer_shade_06.csv')
	df12 = pd.read_csv(filename,header=0)
	df12.index = dates_shade
	df12 = df12.between_time(shour,ehour)

	gs_sun = np.empty([len(df1),6])
	gs_sun[:,0] = (df1['gs(mmol.m-2.s-1)'].values*df1['la(m2/m2)'].values)
	gs_sun[:,1] = (df2['gs(mmol.m-2.s-1)'].values*df2['la(m2/m2)'].values)
	gs_sun[:,2] = (df3['gs(mmol.m-2.s-1)'].values*df3['la(m2/m2)'].values)
	gs_sun[:,3] = (df4['gs(mmol.m-2.s-1)'].values*df4['la(m2/m2)'].values)
	gs_sun[:,4] = (df5['gs(mmol.m-2.s-1)'].values*df5['la(m2/m2)'].values)
	gs_sun[:,5] = (df6['gs(mmol.m-2.s-1)'].values*df6['la(m2/m2)'].values)

	Gs_sun = pd.DataFrame(gs_sun,index=df1.index,columns=['layer1','layer2','layer3','layer4','layer5','layer6'])

	gs_sha = np.empty([len(df7),6])
	gs_sha[:,0] =(df7['gs(mmol.m-2.s-1)'].values*df7['la(m2/m2)'].values)
	gs_sha[:,1] =(df8['gs(mmol.m-2.s-1)'].values*df8['la(m2/m2)'].values)
	gs_sha[:,2] =(df9['gs(mmol.m-2.s-1)'].values*df9['la(m2/m2)'].values)
	gs_sha[:,3] =(df10['gs(mmol.m-2.s-1)'].values*df10['la(m2/m2)'].values)
	gs_sha[:,4] =(df11['gs(mmol.m-2.s-1)'].values*df11['la(m2/m2)'].values)
	gs_sha[:,5] =(df12['gs(mmol.m-2.s-1)'].values*df12['la(m2/m2)'].values)

	Gs_shade = pd.DataFrame(gs_sha,index=df7.index,columns=['layer1','layer2','layer3','layer4','layer5','layer6'])


	return Gs_sun, Gs_shade


def readin_leaflayers(in_dir,variable,startdate,factor=1):
	filename = os.path.join(in_dir+'layer_sun_01.csv')
	df1 = pd.read_csv(filename,header=0,usecols=[variable,'Time (days)'])
	# get date time                                     
	dates = []
	for ii in range(len(df1)):   
		delta = timedelta(df1['Time (days)'][ii]) 
		dt = startdate + delta
		rdate = round_time(dt,round_to=60*30)     
		dates.append(rdate)
	df1.index = dates
	filename = os.path.join(in_dir+'layer_sun_02.csv')
	df2 = pd.read_csv(filename,header=0,usecols=[variable])
	df2.index = dates
	filename = os.path.join(in_dir+'layer_sun_03.csv')
	df3 = pd.read_csv(filename,header=0,usecols=[variable])
	df3.index = dates
	filename = os.path.join(in_dir+'layer_sun_04.csv')
	df4 = pd.read_csv(filename,header=0,usecols=[variable])
	df4.index = dates
	filename = os.path.join(in_dir+'layer_sun_05.csv')
	df5 = pd.read_csv(filename,header=0,usecols=[variable])
	df5.index = dates
	filename = os.path.join(in_dir+'layer_sun_06.csv')
	df6 = pd.read_csv(filename,header=0,usecols=[variable])
	df6.index = dates

	if (factor==1):
		total_sun = np.sum([df1[variable].values,df2[variable].values,df3[variable].values,df4[variable].values,df5[variable].values,df6[variable].values], axis=0)
	elif (factor==2):
		total_sun = np.mean([df1[variable].values,df2[variable].values,df3[variable].values,df4[variable].values,df5[variable].values,df6[variable].values], axis=0)
	
	canopy_sun = pd.DataFrame(total_sun,index=df1.index,columns=['total'])

		# readin shade layers
	filename = os.path.join(in_dir+'layer_shade_01.csv')
	df7 = pd.read_csv(filename,header=0,usecols=[variable,'Time (days)'])
	# get date time                                     
	dates_shade = []
	for ii in range(len(df7)):   
		delta = timedelta(df7['Time (days)'][ii]) 
		dt = startdate + delta
		rdate = round_time(dt,round_to=60*30)     
		dates_shade.append(rdate)
	df7.index = dates_shade
	filename = os.path.join(in_dir+'layer_shade_02.csv')
	df8 = pd.read_csv(filename,header=0,usecols=[variable])
	df8.index = dates_shade
	filename = os.path.join(in_dir+'layer_shade_03.csv')
	df9 = pd.read_csv(filename,header=0,usecols=[variable])
	df9.index = dates_shade
	filename = os.path.join(in_dir+'layer_shade_04.csv')
	df10 = pd.read_csv(filename,header=0,usecols=[variable])
	df10.index = dates_shade
	filename = os.path.join(in_dir+'layer_shade_05.csv')
	df11 = pd.read_csv(filename,header=0,usecols=[variable])
	df11.index = dates_shade
	filename = os.path.join(in_dir+'layer_shade_06.csv')
	df12 = pd.read_csv(filename,header=0,usecols=[variable])
	df12.index = dates_shade

	if (factor==1):
		total_sha = np.sum([df7[variable].values,df8[variable].values,df9[variable].values,df10[variable].values,df11[variable].values,df12[variable].values], axis=0)
	elif(factor==2):
		total_sha = np.mean([df7[variable].values,df8[variable].values,df9[variable].values,df10[variable].values,df11[variable].values,df12[variable].values], axis=0)
	
	canopy_sha = pd.DataFrame(total_sha,index=df7.index,columns=['total'])
	canopy_sun = pd.DataFrame(total_sun,index=df1.index,columns=['total'])

	tmp = canopy_sha.reindex_like(canopy_sun)
	canopy_total = pd.DataFrame(canopy_sun.total+tmp.total,index=df1.index,columns=['total'])

	return canopy_total


def weighted_leaflayers(in_dir,variable,weight,startdate,enddate):
	trng_hh = pd.date_range(startdate,enddate,freq='30min')
	filename = os.path.join(in_dir+'layer_sun_01.csv')
	df1 = pd.read_csv(filename,header=0,usecols=[variable, weight,'Time (days)'])
	canopy = np.empty([len(df1),12])
	# get date time                                     
	dates = []
	for ii in range(len(df1)):   
		delta = timedelta(df1['Time (days)'][ii]) 
		dt = startdate + delta
		rdate = round_time(dt,round_to=60*30)     
		dates.append(rdate)
	df1.index = dates
	filename = os.path.join(in_dir+'layer_sun_02.csv')
	df2 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df2.index = dates
	filename = os.path.join(in_dir+'layer_sun_03.csv')
	df3 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df3.index = dates
	filename = os.path.join(in_dir+'layer_sun_04.csv')
	df4 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df4.index = dates
	filename = os.path.join(in_dir+'layer_sun_05.csv')
	df5 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df5.index = dates
	filename = os.path.join(in_dir+'layer_sun_06.csv')
	df6 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df6.index = dates

	# readin shade layers
	filename = os.path.join(in_dir+'layer_shade_01.csv')
	df7 = pd.read_csv(filename,header=0,usecols=[variable, weight,'Time (days)'])
	# get date time                                     
	dates_shade = []
	for ii in range(len(df7)):   
		delta = timedelta(df7['Time (days)'][ii]) 
		dt = startdate + delta
		rdate = round_time(dt,round_to=60*30)     
		dates_shade.append(rdate)
	df7.index = dates_shade
	df7 = df7.reindex_like(df1)
	filename = os.path.join(in_dir+'layer_shade_02.csv')
	df8 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df8.index = dates_shade
	df8 = df8.reindex_like(df1)
	filename = os.path.join(in_dir+'layer_shade_03.csv')
	df9 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df9.index = dates_shade
	df9 = df9.reindex_like(df1)
	filename = os.path.join(in_dir+'layer_shade_04.csv')
	df10 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df10.index = dates_shade
	df10 = df10.reindex_like(df1)
	filename = os.path.join(in_dir+'layer_shade_05.csv')
	df11 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df11.index = dates_shade
	df11 = df11.reindex_like(df1)
	filename = os.path.join(in_dir+'layer_shade_06.csv')
	df12 = pd.read_csv(filename,header=0,usecols=[variable, weight])
	df12.index = dates_shade
	df12 = df12.reindex_like(df1)

	canopy = np.empty([len(df1),12])
	canopy[:,0] = (df1[variable].values*df1[weight].values)
	canopy[:,1] = (df2[variable].values*df2[weight].values)
	canopy[:,2] = (df3[variable].values*df3[weight].values)
	canopy[:,3] = (df4[variable].values*df4[weight].values)
	canopy[:,4] = (df5[variable].values*df5[weight].values)
	canopy[:,5] = (df6[variable].values*df6[weight].values)
	canopy[:,6] =(df7[variable].values*df7[weight].values)
	canopy[:,7] =(df8[variable].values*df8[weight].values)
	canopy[:,8] =(df9[variable].values*df9[weight].values)
	canopy[:,9] =(df10[variable].values*df10[weight].values)
	canopy[:,10] =(df11[variable].values*df11[weight].values)
	canopy[:,11] =(df12[variable].values*df12[weight].values)

	canopy_total = np.sum(canopy,axis=1)
	weight_total = df1[weight]+df2[weight]+df3[weight]+df4[weight]+df5[weight]+df6[weight]+df7[weight]+df8[weight]+df9[weight]+df10[weight]+df11[weight]+df12[weight]

	fluxweight_var = np.divide(canopy_total,weight_total)

	return pd.DataFrame({variable:fluxweight_var})


def vpd_from_rh(T,RH):
    """Calculate Vapor Pressure Deficit (kPa) from Temperature (C) and RH (%)
    Parameters : T   : Temperature (C)
                 RH  : Relative Humidity (%)
    Returns    : VPD : vapor pressure deficit (hPa)"""

    for x in range(0,len(RH)):
        if RH[x]<0:
            RH[x]=40

    # convert input variables
    TK = T+273.15 # convert temperature from celcius to Kelvin

    # Calculate VPD using Emanuel 1994 equations for saturated vapor pressure (es)
    es = np.exp(53.67957 - (6743.769/TK) - 4.8451*np.log(TK))
    ea = (RH/100)*es
    VPD = es-ea
    VPD = VPD/10 # convert hPa to kPa

    return VPD


def fit_expo_ab(xarray,yarray,upper):
    """Calculate exponential fit for dry-down events
    Parameters : xarray   : temperature (C)
                 yarray   : Net incoming solar radiation (W/m2)
    Returns    : PET : Penman montieth potential evapotranspiration (mm/day)"""
    xarray2 = np.arange(0,upper,0.05)

    def exponential(x, a, b):
        return a*np.exp(x*b)

    popt_exp, pcov_exp = scipy.optimize.curve_fit(exponential, xarray, yarray, p0=[1, -0.5], maxfev=5000)
    fit = (popt_exp[0]*np.exp(popt_exp[1]*xarray2))
    tau = (-1/popt_exp[1])

    return [fit, tau]
