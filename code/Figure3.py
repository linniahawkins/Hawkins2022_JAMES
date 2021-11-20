#!/usr/bin/env python2.7

# Readin daily output from SPA
 
import sys
import os, glob
import csv 
import numpy as np
import calendar
import matplotlib
import matplotlib.pyplot as plt 
from datetime import datetime, timedelta
import pandas as pd
import seaborn as sns
import xarray as xr
from pyfunctions import *

st = datetime(2006,1,1,0,0,0)
en = datetime(2018,12,31,23,59,59)

##############################################
########## readin observed data ##############
stdate=datetime(2006,1,1,0,0,0); endate = datetime(2018,12,31,23,59,59)
trng_hh = pd.date_range(stdate, endate, freq='30min') 
in_file = '../obs_data/US-Me2_AmeriFlux_2006-01-01_2018-12-31_30min.csv'
df_obs = pd.read_csv(in_file,index_col=0,header=0, parse_dates=True, squeeze=True)
df_obs['mon'] = df_obs.index.month
df_obs = df_obs[df_obs['mon']==7]
##########################################
# find days where daily maximum VPD is >90th percentile

VPD = pd.DataFrame(df_obs['VPD'].resample('D').max())
vpd_90 = np.nanquantile(VPD['VPD'],.9)
vpd_75 = np.nanquantile(VPD['VPD'],.75)
vpd_50 = np.nanquantile(VPD['VPD'],.50)

VPD_high = pd.DataFrame(VPD[VPD['VPD']>vpd_75])
VPD_low = pd.DataFrame(VPD[VPD['VPD']<vpd_50])

soilm_75 = np.nanquantile(df_obs['SWP'],.75)
soilm_50 = np.nanquantile(df_obs['SWP'],.5)
soilm_25 = np.nanquantile(df_obs['SWP'],.25)

swp_obs = df_obs['SWP']
soilm_low = pd.DataFrame(swp_obs[swp_obs<soilm_25])
soilm_high = pd.DataFrame(swp_obs[swp_obs>soilm_50])

# vpd soilm
dry_dry = soilm_low.reindex(VPD_high.index).dropna()
print(np.shape(dry_dry))
dry_wet = soilm_high.reindex(VPD_high.index).dropna()
print(np.shape(dry_wet))
wet_wet = soilm_high.reindex(VPD_low.index).dropna()
print(np.shape(wet_wet))
wet_dry = soilm_low.reindex(VPD_low.index).dropna()
print(np.shape(wet_dry))


sel_dates = dry_dry.index
sv_ids = pd.DataFrame(columns = ['dates'])
for ij in range(len(sel_dates)):
    day = sel_dates[ij]
    day_range = pd.date_range(datetime(day.year,day.month,day.day,0,0,0),datetime(day.year,day.month,day.day,23,59,59),freq='30min')
    sv_ids = pd.concat([sv_ids,pd.DataFrame(day_range)])
sv_ids.columns=['dates','value']
sv_dd_ids = sv_ids.set_index('dates')

sel_dates = dry_wet.index
sv_ids = pd.DataFrame(columns = ['dates'])
for ij in range(len(sel_dates)):
    day = sel_dates[ij]
    day_range = pd.date_range(datetime(day.year,day.month,day.day,0,0,0),datetime(day.year,day.month,day.day,23,59,59),freq='30min')
    sv_ids = pd.concat([sv_ids,pd.DataFrame(day_range)])
sv_ids.columns=['dates','value']
sv_dw_ids = sv_ids.set_index('dates')

sel_dates = wet_wet.index
sv_ids = pd.DataFrame(columns = ['dates'])
for ij in range(len(sel_dates)):
    day = sel_dates[ij]
    day_range = pd.date_range(datetime(day.year,day.month,day.day,0,0,0),datetime(day.year,day.month,day.day,23,59,59),freq='30min')
    sv_ids = pd.concat([sv_ids,pd.DataFrame(day_range)])
sv_ids.columns=['dates','value']
sv_ww_ids = sv_ids.set_index('dates')

sel_dates = wet_dry.index
sv_ids = pd.DataFrame(columns = ['dates'])
for ij in range(len(sel_dates)):
    day = sel_dates[ij]
    day_range = pd.date_range(datetime(day.year,day.month,day.day,0,0,0),datetime(day.year,day.month,day.day,23,59,59),freq='30min')
    sv_ids = pd.concat([sv_ids,pd.DataFrame(day_range)])
sv_ids.columns=['dates','value']
sv_wd_ids = sv_ids.set_index('dates')

df_obs = df_obs.dropna(subset=['Trans(mm/day)'])
df_obs_dd = df_obs.reindex(sv_dd_ids.index)
df_obs_dw = df_obs.reindex(sv_dw_ids.index)
df_obs_ww = df_obs.reindex(sv_ww_ids.index)
df_obs_wd = df_obs.reindex(sv_wd_ids.index)

######### create data sets ############
# high VPD low soilm
data1 = df_obs_dd
ds_data1 = xr.Dataset(data1)
dat = ds_data1.groupby("dates.hour").mean()
data_dd = dat.to_dataframe()

# high VPD high soilm
data1 = df_obs_dw
ds_data1 = xr.Dataset(data1)
dat = ds_data1.groupby("dates.hour").mean()
data_dw = dat.to_dataframe()

# low VPD high soilm
data1 = df_obs_ww
ds_data1 = xr.Dataset(data1)
dat = ds_data1.groupby("dates.hour").mean()
data_ww = dat.to_dataframe()

# low VPD low soilm
data1 = df_obs_wd
ds_data1 = xr.Dataset(data1)
dat = ds_data1.groupby("dates.hour").mean()
data_wd = dat.to_dataframe()

##################################
# Readin Models
##################################

min_syr = 2006; max_eyr = 2018
start = datetime(min_syr,1,1,0,0,0); end=datetime(max_eyr,12,31,23,59,59); freq='30min'

#readin WUEi
in_dir='../model_data/SPA-WUEi/'  
data_wuei = SPA_readin_ecosystemfluxes(in_dir,start,end,freq)
ds_data1 = xr.Dataset(data_wuei.loc[sv_dd_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_wuei_dd = dat.to_dataframe()


ds_data1 = xr.Dataset(data_wuei.loc[sv_dw_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_wuei_dw = dat.to_dataframe()


ds_data1 = xr.Dataset(data_wuei.loc[sv_ww_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_wuei_ww = dat.to_dataframe()


ds_data1 = xr.Dataset(data_wuei.loc[sv_wd_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_wuei_wd = dat.to_dataframe()


#readin WUE
in_dir='../model_data/SPA-WUE/'
data_wue = SPA_readin_ecosystemfluxes(in_dir,start,end,freq)
ds_data1 = xr.Dataset(data_wue.loc[sv_dd_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_wue_dd = dat.to_dataframe()

ds_data1 = xr.Dataset(data_wue.loc[sv_dw_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_wue_dw = dat.to_dataframe()

ds_data1 = xr.Dataset(data_wue.loc[sv_ww_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_wue_ww = dat.to_dataframe()

ds_data1 = xr.Dataset(data_wue.loc[sv_wd_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_wue_wd = dat.to_dataframe()



# readin BB
in_dir='../model_data/SPA-BB-H/'
data_bb = SPA_readin_ecosystemfluxes(in_dir,start,end,freq)
ds_data1 = xr.Dataset(data_bb.loc[sv_dd_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_bb_dd = dat.to_dataframe()

ds_data1 = xr.Dataset(data_bb.loc[sv_dw_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_bb_dw = dat.to_dataframe()

ds_data1 = xr.Dataset(data_bb.loc[sv_ww_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_bb_ww = dat.to_dataframe()

ds_data1 = xr.Dataset(data_bb.loc[sv_wd_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_bb_wd = dat.to_dataframe()


#readin med
in_dir='../model_data/SPA-MED-H/'  
data_med = SPA_readin_ecosystemfluxes(in_dir,start,end,freq)
ds_data1 = xr.Dataset(data_med.loc[sv_dd_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_med_dd = dat.to_dataframe()

ds_data1 = xr.Dataset(data_med.loc[sv_dw_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_med_dw = dat.to_dataframe()

ds_data1 = xr.Dataset(data_med.loc[sv_ww_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_med_ww = dat.to_dataframe()

ds_data1 = xr.Dataset(data_med.loc[sv_wd_ids.index.values])
dat = ds_data1.groupby("dim_0.hour").mean()
data_med_wd = dat.to_dataframe()



# readin gain-risk
in_dir='../model_data/gain-risk/'   
start = datetime(2006,1,1,0,0,0); end=datetime(2018,12,31,23,59,59)
trng_h = pd.date_range(start, end, freq='H') #'30min','H'
filename = os.path.join(in_dir+'standard_OUTPUT_timesteps.csv')
dfc2 = pd.read_csv(filename,header=0,index_col=None,skiprows=None)
dfc2.index = trng_h

ds_data1 = xr.Dataset(dfc2.reindex(sv_dd_ids.index))
dat = ds_data1.groupby("dates.hour").mean()
data_sperry_dd = dat.to_dataframe()

ds_data1 = xr.Dataset(dfc2.reindex(sv_dw_ids.index))
dat = ds_data1.groupby("dates.hour").mean()
data_sperry_dw = dat.to_dataframe()

ds_data1 = xr.Dataset(dfc2.reindex(sv_ww_ids.index))
dat = ds_data1.groupby("dates.hour").mean()
data_sperry_ww = dat.to_dataframe()

ds_data1 = xr.Dataset(dfc2.reindex(sv_wd_ids.index))
dat = ds_data1.groupby("dates.hour").mean()
data_sperry_wd = dat.to_dataframe()

################# plot #######################
plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.1, right=.9, top=.95, wspace=0.2, hspace=0.2)
plt.rcParams.update({'font.size': 20})

plt.subplot(2,2,1)
plt.fill_between(data_dd['Trans(mm/day)'].index, data_dd['Trans(mm/day)']*0.6, data_dd['Trans(mm/day)']*1.4, facecolor='k', alpha=0.3)
plt.plot(data_dd['Trans(mm/day)'],color='k',label='Sapflow',linewidth=4,linestyle='dotted')
plt.plot(data_wuei_dd['trans(W/m2)']/28.94,color='sienna',linewidth=4,label='WUEi')
plt.plot(data_wue_dd['trans(W/m2)']/28.94,color='darkgoldenrod',linewidth=4,label='WUE',linestyle='--')
plt.plot(data_bb_dd['trans(W/m2)']/28.94,color='olivedrab',linewidth=4,label='BB-H',linestyle='-')
plt.plot(data_med_dd['trans(W/m2)']/28.94,color='seagreen',linewidth=4,label='MED-H',linestyle='--')
plt.plot(data_sperry_dd['E (mm timestep-1)']*24,color='darkcyan',alpha=0.7,linewidth=4,label='Gain-Risk')
plt.ylim([0,6.5])
plt.xlim([4.5,20])
plt.ylabel('Transpiration\n(mm/day)')
plt.legend(loc='upper right',fontsize=14)

plt.subplot(2,2,2)
plt.fill_between(data_dw['Trans(mm/day)'].index, data_dw['Trans(mm/day)']*0.6, data_dw['Trans(mm/day)']*1.4, facecolor='k', alpha=0.3)
plt.plot(data_dw['Trans(mm/day)'],color='k',linewidth=4,linestyle='dotted')
plt.plot(data_wuei_dw['trans(W/m2)']/28.94,linewidth=4,color='sienna')
plt.plot(data_wue_dw['trans(W/m2)']/28.94,linewidth=4,color='darkgoldenrod',linestyle='--')
plt.plot(data_bb_dw['trans(W/m2)']/28.94,linewidth=4,color='olivedrab',linestyle='-')
plt.plot(data_med_dw['trans(W/m2)']/28.94,linewidth=4,color='seagreen',linestyle='--')
plt.plot(data_sperry_dw['E (mm timestep-1)']*24,linewidth=4,alpha=0.7,color='darkcyan')
plt.ylim([0,6.5])
plt.xlim([4.5,20])

plt.subplot(2,2,3)
plt.fill_between(data_wd['Trans(mm/day)'].index, data_wd['Trans(mm/day)']*0.6, data_wd['Trans(mm/day)']*1.4, facecolor='k', alpha=0.3)
plt.plot(data_wd['Trans(mm/day)'],color='k',label='Sapflow',linewidth=4,linestyle='dotted')
plt.plot(data_wuei_wd['trans(W/m2)']/28.94,linewidth=4,color='sienna',label='WUEi')
plt.plot(data_wue_wd['trans(W/m2)']/28.94,linewidth=4,color='darkgoldenrod',label='WUE',linestyle='--')
plt.plot(data_bb_wd['trans(W/m2)']/28.94,linewidth=4,color='olivedrab',label='Ball-Berry',linestyle='-')
plt.plot(data_med_wd['trans(W/m2)']/28.94,linewidth=4,color='seagreen',label='Medlyn',linestyle='--')
plt.plot(data_sperry_wd['E (mm timestep-1)']*24,linewidth=4,color='darkcyan',label='Gain-Risk',alpha=0.7)
plt.ylim([0,6.5])
plt.xlim([4.5,20])
plt.ylabel('Transpiration\n(mm/day)')
plt.xlabel('hour of day')

plt.subplot(2,2,4)
plt.fill_between(data_ww['Trans(mm/day)'].index, data_ww['Trans(mm/day)']*0.6, data_ww['Trans(mm/day)']*1.4, facecolor='k', alpha=0.3)
plt.plot(data_ww['Trans(mm/day)'],color='k',linestyle='dotted',linewidth=4)
plt.plot(data_wuei_ww['trans(W/m2)']/28.94,linewidth=4,color='sienna')
plt.plot(data_wue_ww['trans(W/m2)']/28.94,linewidth=4,color='darkgoldenrod',linestyle='--')
plt.plot(data_bb_ww['trans(W/m2)']/28.94,linewidth=4,color='olivedrab',linestyle='-')
plt.plot(data_med_ww['trans(W/m2)']/28.94,linewidth=4,color='seagreen',linestyle='--')
plt.plot(data_sperry_ww['E (mm timestep-1)']*24,linewidth=4,color='darkcyan',alpha=0.7)
plt.ylim([0,6.5])
plt.xlim([4.5,20])
plt.xlabel('hour of day')


plotname = 'Figure3'
plt.savefig(plotname, dpi=300)

plt.show()

