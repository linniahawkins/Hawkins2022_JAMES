#!/usr/bin/env python2.7
 
import sys
import os, glob
import csv 
import numpy as np
import calendar
import matplotlib
import matplotlib.pyplot as plt 
from datetime import datetime, date, timedelta
import pandas as pd
import scipy
from scipy.optimize import curve_fit
from pyfunctions import *


def fit_model_tau(in_dir, obs_daily, startdate, enddate, soilm_ptile_low, soilm_ptile_hi,modid=1): 

	Gs_sun, Gs_shade = readin_leaflayers_Gs(in_dir,startdate,'8:00','16:00')
	Gs = weighted_leaflayers(in_dir,'gs(mmol.m-2.s-1)','la(m2/m2)',startdate,enddate)

	Gs_d = Gs.between_time('8:00','16:00').resample('D').mean()
	s = Gs_d.index[0]; e = Gs_d.index[-1]

	######### readin ecosystem fluxes ################ 
	trng_hh = pd.date_range(startdate, enddate, freq='30min')
	 
	filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
	data_eco = pd.read_csv(filename,header=0)

	data_eco.index = trng_hh
	data_eco = data_eco.between_time('8:00','16:00').resample('D').mean()

	########## build dataframe #############

	trng_d = pd.date_range(s,e,freq='D')

	data = pd.DataFrame({'Gs' : Gs_d['gs(mmol.m-2.s-1)'].values,  'vpd' : obs_daily['VPD'][s:e].values, 'swp' : obs_daily['SWP'][s:e].values, 'trans':data_eco['trans(W/m2)'][s:e].values/28.94,'gpp':data_eco['gpp(umol/m2/s)'][s:e].values},index = trng_d)
	data['mon'] = data.index.month
	data = data[data.mon > 5]
	data = data[data.mon < 9]
	if modid > 1:
		Gs_max = 2*np.nanmax(data.Gs)
	else:
		Gs_max = np.nanmax(data.Gs)

	Gs_scaled = data.Gs/Gs_max
	data['Gs_scaled'] = Gs_scaled
	data[data['vpd']<1] = np.NaN
	tmp_dat = data.dropna()

	ptile_low = np.nanpercentile(tmp_dat['swp'],soilm_ptile_low)
	ptile_hi = np.nanpercentile(tmp_dat['swp'],soilm_ptile_hi)

	tmp_dat = tmp_dat.loc[tmp_dat['swp']<ptile_hi]
	tmp_dat = tmp_dat.loc[tmp_dat['swp']>ptile_low]

	x_array = tmp_dat['vpd']
	y_array = tmp_dat['Gs_scaled']
	y_array[tmp_dat['vpd']>3] = 0

	model_fit, model_tau = fit_expo_ab(x_array,y_array,3.2)
	xarray = np.arange(0,3.2,0.05)

	out_data = pd.DataFrame({'fit':model_fit,'tau':model_tau,'xarray':xarray})
	out_data = out_data.sort_values('xarray')

	return out_data


def fit_sperry_tau(in_file,obs_daily,startdate,enddate, soilm_ptile_low,soilm_ptile_hi):

	index_h = pd.date_range(startdate, enddate, freq='H') 
	index_d = pd.date_range(startdate, enddate, freq='D') 
	dfc2 = pd.read_csv(in_file,header=0,index_col=None,skiprows=None)
	dfc2.index = index_h
	dfc2[dfc2['Gw, mmol m-2s-1']>1000] = np.nan

	dfc2 = dfc2.between_time('8:00','16:00').resample('D').mean()
	data_sperry = pd.DataFrame({'Gs' : dfc2['Gw, mmol m-2s-1'].values,  'vpd' : obs_daily['VPD'].values,'swp' : obs_daily['SWP'].values},index = index_d)
	data_sperry[data_sperry['vpd']<.1]=np.nan
	data_sperry['mon'] = data_sperry.index.month

	data_sperry = data_sperry[data_sperry.mon > 5]
	data_sperry = data_sperry[data_sperry.mon < 9]

	data_sperry['Gs_scaled'] = data_sperry['Gs']/(np.nanmax(data_sperry['Gs']))
	data_sperry['Gs_scaled'] = data_sperry['Gs'] #

	tmp_dat = data_sperry.dropna()

	ptile_low = np.nanpercentile(tmp_dat['swp'],soilm_ptile_low)
	ptile_hi = np.nanpercentile(tmp_dat['swp'],soilm_ptile_hi)

	tmp_dat = tmp_dat.loc[tmp_dat['swp']<ptile_hi]
	tmp_dat = tmp_dat.loc[tmp_dat['swp']>ptile_low]

	x_array = tmp_dat['vpd']
	y_array = tmp_dat['Gs_scaled']

	out = fit_expo_ab(x_array,y_array,3.2)
	xarray = np.arange(0,3.2,0.05)

	out_data = pd.DataFrame({'fit':out[0],'tau':out[1],'xarray':xarray})
	out_data = out_data.sort_values('xarray')

	return out_data

def fit_obs_tau(obs,startdate,enddate,soilm_ptile_low, soilm_ptile_hi,error='mean'):

	vpd = obs['VPD'].values
	Ta = obs['TA'].values
	T = obs['Trans(mm/day)'].values/3.2 # per unit LAI 
	T = T/86400*1000 # convert mm/day to g/m2/s
	Kg = 115.8 + 0.4236*Ta # (KPa m3 / Kg)

	Gs = Kg*(T/vpd)*40.1 # convert mm/s to mmol/m2/s
	Gs_lo = Kg*(T*0.6/vpd)*40.1
	Gs_hi = Kg*(T*1.4/vpd)*40.1

	index_hh = pd.date_range(startdate, enddate, freq='30min') 

	df_Gs = pd.DataFrame({'Gs': Gs,'Gs_lo':Gs_lo,'Gs_hi':Gs_hi,'vpd':obs['VPD'].values,'swp':obs['SWP'].values},index = index_hh)
	df_Gs = df_Gs.between_time('8:00','16:00').resample('D').mean()
	df_Gs['mon'] = df_Gs.index.month
	df_Gs = df_Gs[df_Gs.mon > 5]
	df_Gs = df_Gs[df_Gs.mon < 9]
	Gs_max = np.nanmax(df_Gs.Gs)
	df_Gs['Gs_scaled'] = df_Gs.Gs/Gs_max
	df_Gs['Gs_lo_scaled'] = df_Gs.Gs_lo/Gs_max
	df_Gs['Gs_hi_scaled'] = df_Gs.Gs_hi/Gs_max

	df_Gs[df_Gs['Gs_scaled']<.01] = np.NaN
	tmp_dat = df_Gs.dropna()

	ptile_low = np.nanpercentile(tmp_dat['swp'],soilm_ptile_low)
	ptile_hi = np.nanpercentile(tmp_dat['swp'],soilm_ptile_hi)

	tmp_dat = tmp_dat.loc[tmp_dat['swp']<ptile_hi]
	tmp_dat = tmp_dat.loc[tmp_dat['swp']>ptile_low]

	x_array = tmp_dat['vpd']

	if (error == 'mean'):
		y_array = tmp_dat['Gs_scaled']
	elif (error == 'upper'):
		y_array = tmp_dat['Gs_hi_scaled']
	elif (error == 'lower'):
		y_array = tmp_dat['Gs_lo_scaled']
	else:
		print('error')

	obs_fit, obs_tau = fit_expo_ab(x_array,y_array,3.2)
	xarray = np.arange(0,3.2,0.05)

	out_data = pd.DataFrame({'fit':obs_fit,'tau':obs_tau,'xarray':xarray})
	out_data = out_data.sort_values('xarray')

	return out_data


if __name__ == "__main__":

	st = datetime(2006,1,1,0,0,0) 
	en = datetime(2018,12,31,23,59,59)

	##############################################
	########## readin observed data ##############
	stdate=datetime(2006,1,1,0,0,0); endate = datetime(2018,12,31,23,59,59)
	trng_hh = pd.date_range(stdate, endate, freq='30min') 
	in_file = '../obs_data/US-Me2_AmeriFlux_2006-01-01_2018-12-31_30min.csv'
	obs = pd.read_csv(in_file,index_col=0,header=0, parse_dates=True, squeeze=True)
	obs_daily = obs.between_time('8:00','16:00').resample('D').mean()

	soilm_ptile_low = 75
	soilm_ptile_hi = 100
	in_dir='../model_data/SPA-WUEi/'  
	df_wuei = fit_model_tau(in_dir,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi,modid=2)

	in_dir='../model_data/SPA-WUE/'  
	df_wue = fit_model_tau(in_dir,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi)

	in_dir='../model_data/SPA-BB-H/'  
	df_bb = fit_model_tau(in_dir,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi)

	in_dir='../model_data/SPA-MED-H/'  
	df_med = fit_model_tau(in_dir,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi)

	in_file = '../model_data/gain-risk/standard_OUTPUT_timesteps.csv'
	df_sperry = fit_sperry_tau(in_file,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi)

	df_obs = fit_obs_tau(obs,stdate,endate,soilm_ptile_low,soilm_ptile_hi,error='mean')
	df_obs_lower = fit_obs_tau(obs,stdate,endate,soilm_ptile_low,soilm_ptile_hi,error='lower')
	df_obs_upper = fit_obs_tau(obs,stdate,endate,soilm_ptile_low,soilm_ptile_hi,error='upper')

	####### plot ################
	plt.figure(num=None, figsize=(13, 5), dpi=100, facecolor='w', edgecolor='k')
	plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9, wspace=0.3, hspace=0.3)
	plt.rcParams.update({'font.size': 20})

	ax1 = plt.subplot(1,2,1)
	ax1.fill_between(df_obs_lower['xarray'], df_obs_lower['fit'].values, df_obs_upper['fit'].values, facecolor='grey', alpha=0.5)
	ax1.plot(df_obs['xarray'],df_obs['fit'],linewidth=4,color='k',label='Observed',linestyle='dotted')
	ax1.plot(df_wuei['xarray'],df_wuei['fit'],color='sienna',linewidth=4,label='WUEi')
	ax1.plot(df_wue['xarray'],df_wue['fit'],color='darkgoldenrod',linewidth=4,label='WUE',linestyle='--')
	ax1.plot(df_bb['xarray'],df_bb['fit'],color='olivedrab',linewidth=4,label='BB-H')
	ax1.plot(df_med['xarray'],df_med['fit'],color='seagreen',linewidth=4,label='MED-H',linestyle='--')
	ax1.plot(df_sperry['xarray'], df_sperry['fit'],color='darkcyan',linewidth=4,label='Gain-Risk')
	ax1.set_xlabel('VPD (kPa)')
	ax1.set_ylabel(r'$G_{c}/G_{c,max}$')
	ax1.set_ylim([0,0.7])
	ax1.text(0.45, 0.95, '(a) SWP>75th', transform=ax1.transAxes,fontsize=18, va='top', ha='right')

	soilm_ptile_low = 0
	soilm_ptile_hi = 25
	in_dir='../model_data/SPA-WUEi/'  
	df_wuei = fit_model_tau(in_dir,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi,modid=2)

	in_dir='../model_data/SPA-WUE/'  
	df_wue = fit_model_tau(in_dir,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi)

	in_dir='../model_data/SPA-BB-H/'  
	df_bb = fit_model_tau(in_dir,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi)

	in_dir='../model_data/SPA-MED-H/'  
	df_med = fit_model_tau(in_dir,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi)

	in_file = '../model_data/gain-risk/standard_OUTPUT_timesteps.csv'
	df_sperry = fit_sperry_tau(in_file,obs_daily,stdate,endate,soilm_ptile_low,soilm_ptile_hi)

	df_obs = fit_obs_tau(obs,stdate,endate,soilm_ptile_low,soilm_ptile_hi,error='mean')
	df_obs_lower = fit_obs_tau(obs,stdate,endate,soilm_ptile_low,soilm_ptile_hi,error='lower')
	df_obs_upper = fit_obs_tau(obs,stdate,endate,soilm_ptile_low,soilm_ptile_hi,error='upper')

	ax2 = plt.subplot(1,2,2)
	ax2.fill_between(df_obs_lower['xarray'], df_obs_lower['fit'].values, df_obs_upper['fit'].values, facecolor='grey', alpha=0.5)
	ax2.plot(df_obs['xarray'],df_obs['fit'],linewidth=4,color='k',label='Observed',linestyle='dotted')
	ax2.plot(df_wuei['xarray'],df_wuei['fit'],color='sienna',linewidth=4,label='WUEi')
	ax2.plot(df_wue['xarray'],df_wue['fit'],color='darkgoldenrod',linewidth=4,label='WUE',linestyle='--')
	ax2.plot(df_bb['xarray'],df_bb['fit'],color='olivedrab',linewidth=4,label='BB-H')
	ax2.plot(df_med['xarray'],df_med['fit'],color='seagreen',linewidth=4,label='MED-H',linestyle='--')
	ax2.plot(df_sperry['xarray'], df_sperry['fit'],color='darkcyan',linewidth=4,label='Gain-Risk')
	ax2.set_xlabel('VPD (kPa)')
	ax2.set_ylabel(r'$G_{c}/G_{c,max}$')
	ax2.set_ylim([0,0.7])
	plt.legend(fontsize=14) 
	ax2.text(0.45, 0.95, '(b) SWP<25th', transform=ax2.transAxes,fontsize=18, va='top', ha='right')

	plt.show()
