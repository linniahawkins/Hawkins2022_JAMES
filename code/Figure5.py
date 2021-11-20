#!/usr/bin/env python2.7
 
import sys
import os, glob
import csv 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
from datetime import datetime, date, timedelta
import pandas as pd
from pyfunctions import *

mon = 8

st=datetime(2006,1,1,0,0,0)
en=datetime(2018,12,31,23,59,59)
index_hh = pd.date_range(st, en, freq='30min') #'30min','H'

stdate = datetime(2015,8,1,0,0,0)
endate = datetime(2015,8,30,0,0,0)
m = '5:00'
a = '17:00'

######################################
######### Readin AmeriFlux ###########

in_dir='../obs_data/' 
filename = os.path.join(in_dir+'US-Me2_AmeriFlux_2006-01-01_2018-12-31_30min.csv')
df_obs = pd.read_csv(filename,header=0,index_col=0, parse_dates=True, squeeze=True)
df_obs[df_obs<-9000]=np.NaN

######################################
############ READIN TCAN #############

in_dir='../obs_data/'  
filename = os.path.join(in_dir+'US-Me2_Tcan_2014-01-01_2015-12-31_30min.csv')
df = pd.read_csv(filename,header=0,parse_dates=False)
index_h = pd.date_range(datetime(2015,1,1,0,0,0),datetime(2016,1,1,0,0,0), freq='30min')
df.index = index_h
varlist = ['Tair','Tcan_Avg_corr','RH']
obs = df[varlist][stdate:endate]
obs.index = obs.index - pd.Timedelta(minutes=30)
obs['Tdiff'] = obs['Tcan_Avg_corr'] - obs['Tair']
obs['VPDleaf'] = vpd_from_rh(obs['Tcan_Avg_corr'].values,obs['RH'].values)
obs_h = obs.resample('H').mean()
obs_h['hour'] = obs_h.index.hour

ds_obs = xr.Dataset(obs_h)
dat = ds_obs.groupby("dim_0.hour").mean()
obs_dd = dat.to_dataframe()

##################################
############ readin Sperry #######

in_file = '../model_data/gain-risk/standard_OUTPUT_timesteps.csv'
index_h = pd.date_range(datetime(2006,1,1,0,0,0), datetime(2018,12,31,23,59,59), freq='H') 
data_sperry = pd.read_csv(in_file,header=0,index_col=None,skiprows=None)
data_sperry.index = index_h
data_sperry['Tdiff'] = data_sperry['leaftempt'] - data_sperry['T air, C']
data_sperry = data_sperry.between_time(m,a).reindex(index_h, fill_value=np.nan)

data_sperry_h = data_sperry[stdate:endate].between_time(m,a)

ds_data = xr.Dataset(data_sperry_h)
dat = ds_data.groupby("dim_0.hour").mean()
data_sperry_dd = dat.to_dataframe()

#####################################
############# readin WUE ############
in_dir='../model_data/SPA-WUE/'  
wue_canopy = weighted_leaflayers(in_dir,'tempdf(oC)','agr(umol.m-2.s-1)',st,en)
wue_canopyT = wue_canopy.reindex(index_hh, fill_value=np.nan)
wue_canopyT['Tcan'] = wue_canopyT['tempdf(oC)']+df_obs['TA']
                                
data_wue = wue_canopyT[stdate:endate].between_time(m,a)
ds_data = xr.Dataset(data_wue)
dat = ds_data.groupby("dim_0.hour").mean()
data_wue_dd = dat.to_dataframe()

#####################################
############## readin WUEi ##########
in_dir='../model_data/SPA-WUEi/'   
wuei_canopy = weighted_leaflayers(in_dir,'tempdf(oC)','agr(umol.m-2.s-1)',st,en)
wuei_canopyT = wuei_canopy.reindex(index_hh, fill_value=np.nan)
wuei_canopyT['Tcan'] = wuei_canopyT['tempdf(oC)']+df_obs['TA']
                               
data_wuei = wuei_canopyT[stdate:endate].between_time(m,a)
ds_data = xr.Dataset(data_wuei)
dat = ds_data.groupby("dim_0.hour").mean()
data_wuei_dd = dat.to_dataframe()

#####################################
######## readin Ballberry ###########
in_dir='../model_data/SPA-BB-H/' 
bb_canopy = weighted_leaflayers(in_dir,'tempdf(oC)','agr(umol.m-2.s-1)',st,en)
bb_canopyT = bb_canopy.reindex(index_hh, fill_value=np.nan)
bb_canopyT['Tcan'] = bb_canopyT['tempdf(oC)']+df_obs['TA']

data_bb = bb_canopyT[stdate:endate].between_time(m,a)
ds_data = xr.Dataset(data_bb)
dat = ds_data.groupby("dim_0.hour").mean()
data_bb_dd = dat.to_dataframe()

#####################################
########## readin Medlyn ############
in_dir='../model_data/SPA-MED-H/'  
med_canopy = weighted_leaflayers(in_dir,'tempdf(oC)','agr(umol.m-2.s-1)',st,en)
med_canopyT = med_canopy.reindex(index_hh, fill_value=np.nan)
med_canopyT['Tcan'] = med_canopyT['tempdf(oC)']+df_obs['TA']
                            
data_med = med_canopyT[stdate:endate].between_time(m,a)
ds_data = xr.Dataset(data_med)
dat = ds_data.groupby("dim_0.hour").mean()
data_med_dd = dat.to_dataframe()

######## plot ##############

plt.figure(num=None, figsize=(16, 6), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.15, right=.9, top=.95, wspace=0.3, hspace=0.2)
plt.rcParams.update({'font.size': 18})
ax1 = plt.subplot(1,2,1)
ax1.plot(obs_dd['Tdiff'],color='k',linestyle='dotted',linewidth=3,label='Measured')
ax1.plot(data_wuei_dd['tempdf(oC)'],color='sienna',linewidth=3,label='WUEi')
ax1.plot(data_wue_dd['tempdf(oC)'],color='darkgoldenrod',linewidth=3,label='WUE',linestyle='--')
ax1.plot(data_bb_dd['tempdf(oC)'],color='olivedrab',linewidth=3,label='BB-H')
ax1.plot(data_med_dd['tempdf(oC)'],color='seagreen',linewidth=3,label='MED-H',linestyle='--')
ax1.plot(data_sperry_dd['Tdiff'],color='darkcyan',linewidth=3,label='Gain-Risk')
ax1.set_xlim([2,22])
ax1.set_ylabel('Tcan-Tair (C)')
ax1.legend(loc='lower right',fontsize=14)
ax1.set_xlabel('hour of day')
ax1.text(0.1, 0.95, '(a)', transform=ax1.transAxes,fontsize=20, va='top', ha='right')

ax2 = plt.subplot(1,2,2)
ax2.plot([-6,6],[-6,6],linestyle='-', color='k',alpha=0.7,label='1:1 line')
ax2.scatter(obs['Tdiff'].between_time(m,a),data_wuei['tempdf(oC)'],s=40,color='sienna',label='WUEi',alpha=0.5)
ax2.scatter(obs['Tdiff'].between_time(m,a),data_wue['tempdf(oC)'],s=40,color='darkgoldenrod',label='WUE',alpha=0.5)
ax2.scatter(obs['Tdiff'].between_time(m,a),data_bb['tempdf(oC)'],s=40,color='olivedrab',label='BB-H',alpha=0.5)
ax2.scatter(obs['Tdiff'].between_time(m,a),data_med['tempdf(oC)'],s=40,color='seagreen',label='MED-H',alpha=0.5)
ax2.scatter(obs_h['Tdiff'].between_time(m,a).dropna(),data_sperry_h['Tdiff'],s=40,color='darkcyan',label='Gain-Risk',alpha=0.5)
ax2.set_ylim([-3,3])
ax2.set_ylabel('modeled Tcan-Tair (C)')
ax2.legend(loc='lower right',fontsize=14)
ax2.set_xlabel('measured Tcan-Tair (C)')
ax2.text(0.1, 0.95, '(b)', transform=ax2.transAxes,fontsize=20, va='top', ha='right')

plt.show()
