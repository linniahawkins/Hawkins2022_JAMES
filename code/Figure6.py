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

st=datetime(2015,1,1,0,30,0)
en=datetime(2015,12,31,23,59,59)

stdate = datetime(2015,8,1,0,0,0)
endate = datetime(2015,8,30,0,0,0)
m = '8:00'
a = '16:00'

######################################
########## Readin AmeriFlux ##########

in_dir='../obs_data/' 
filename = os.path.join(in_dir+'US-Me2_AmeriFlux_2006-01-01_2018-12-31_30min.csv')
df_obs = pd.read_csv(filename,header=0,index_col=0, parse_dates=True, squeeze=True)
df_obs[df_obs<-9000]=np.NaN
df_obs = df_obs[stdate:endate]

######################################
############ READIN TCAN #############
in_dir='../obs_data/'  
filename = os.path.join(in_dir+'US-Me2_Tcan_2014-01-01_2015-12-31_30min.csv')
df = pd.read_csv(filename,header=0,parse_dates=False)
index_h = pd.date_range(datetime(2015,1,1,0,0,0),datetime(2016,1,1,0,0,0), freq='30min')
df.index = index_h
obs = df[stdate:endate]
obs.index = obs.index - pd.Timedelta(minutes=30)

############################
#readin WUE
############################
in_dir='../model_data/test_tleaf_med0/'  
filename = os.path.join(in_dir+'layer_sun_01.csv')
data_tcan = pd.read_csv(filename,header=0) 
                                 
startdate = datetime(2015,1,1,0,0,0) 
dates = []
for ii in range(len(data_tcan)):   
	delta = timedelta(data_tcan['Time (days)'][ii]) 
	dt = startdate + delta
	rdate = round_time(dt,round_to=60*30)     
	dates.append(rdate)    
data_tcan.index = dates
data_tcan = data_tcan[stdate:endate]
data_tcan['Tcan'] = data_tcan ['tempdf(oC)']+df_obs['TA'][data_tcan.index]
data_tcan = data_tcan[stdate:endate].between_time(m,a)

T1m = np.cumsum(data_tcan['etr(W.m-2)'].between_time('8:00','12:00').values/28.94) 
T1a = np.cumsum(data_tcan['etr(W.m-2)'].between_time('13:00','16:00').values/28.94) 

T1 = np.cumsum(data_tcan['etr(W.m-2)'].values/28.94) 
print(T1)

in_dir='../model_data/test_tleaf_default_med/'  
filename = os.path.join(in_dir+'layer_sun_01.csv')
data_wue = pd.read_csv(filename,header=0) 

# get date time                                     
startdate = datetime(2015,1,1,0,0,0) 
dates = []
for ii in range(len(data_wue)):   
	delta = timedelta(data_wue['Time (days)'][ii]) 
	dt = startdate + delta
	rdate = round_time(dt,round_to=60*30)     
	dates.append(rdate)    
data_wue.index = dates
data_wue = data_wue[stdate:endate]
data_wue['Tcan'] = data_wue['tempdf(oC)']+df_obs['TA'][data_wue.index]
data_wue2 = data_wue[stdate:endate].between_time(m,a)

T2m = np.cumsum(data_wue2['etr(W.m-2)'].between_time('8:00','12:00').values/28.94)
T2a = np.cumsum(data_wue2['etr(W.m-2)'].between_time('13:00','16:00').values/28.94) 

T2 = np.cumsum(data_wue2['etr(W.m-2)'].values/28.94) 

##################################

plt.figure(num=None, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.17, right=.9, top=.95, wspace=0.3, hspace=0.2)
plt.rcParams.update({'font.size': 16})
ax2 = plt.subplot(1,2,1)
ax2.plot([12,37],[12,37],color='k',label='1:1 line')
plt.scatter(data_tcan['Tcan'],data_wue2['Tcan'],c=data_wue2.index.hour,cmap='RdYlBu_r',alpha=1,label='MED-H')
ax2.text(0.12, 0.95, '(a)', transform=ax2.transAxes,fontsize=18, va='top', ha='right')
cb = plt.colorbar()
cb.ax.set_ylabel('hour of day', fontsize=16)
plt.xlabel(r'$T_{can}$ prescribed')
plt.ylabel(r'$T_{can}$ modeled')
plt.legend(loc='lower right',fontsize=16)

ax2 = plt.subplot(1,2,2)
ax2.plot([0,0.8],[0,0.8],color='k',label='1:1 line')
ax2.text(0.12, 0.95, '(b)', transform=ax2.transAxes,fontsize=18, va='top', ha='right')
plt.scatter(data_tcan['etr(W.m-2)']/28.94,data_wue2['etr(W.m-2)']/28.94,c=data_wue2.index.hour,cmap='RdYlBu_r',alpha=1,label='MED-H')
cb = plt.colorbar()
cb.ax.set_ylabel('hour of day', fontsize=16)
plt.xlabel('Transpiration (mm/day) \n($T_{can}$ prescribed)')
plt.ylabel('Transpiration (mm/day) \n($T_{can}$ modeled)')
plt.legend(loc='lower right',fontsize=16)

plt.figure(num=None, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.1, right=.9, top=.95, wspace=0.2, hspace=0.2)
plt.rcParams.update({'font.size': 16})
plt.subplot(1,2,1)
plt.plot(T1m,color='darkkhaki',alpha=1,label='MED-H-Tcan')
plt.plot(T2m,color='k',alpha=1,label='MED-H')
plt.ylabel('Transpiration')
plt.legend(loc='lower right',fontsize=16)

plt.subplot(1,2,2)
plt.plot(T1a,color='darkkhaki',alpha=1,label='MED-H-Tcan')
plt.plot(T2a,color='k',alpha=1,label='MED-H')
plt.ylabel('Transpiration')
plt.legend(loc='lower right',fontsize=16)

plt.show()
