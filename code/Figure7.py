#!/usr/bin/env python2.7

# Readin daily output from SPA
 
import sys
import os, glob
import csv 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
from datetime import datetime, timedelta
import pandas as pd
import xarray as xr
from pyfunctions import *
from it_metrics import *
import squarify 
import math

st = datetime(2006,1,1,0,0,0)
en = datetime(2018,12,31,23,59,59)

m = '8:00'
a = '16:00'

##############################################
########## readin observed data ##############

startdate=datetime(2006,1,1,0,0,0); enddate = datetime(2018,12,31,23,59,59)
trng_hh = pd.date_range(startdate, enddate, freq='30min') 
in_file = '../obs_data/US-Me2_AmeriFlux_2006-01-01_2018-12-31_30min.csv'
obs_data = pd.read_csv(in_file,index_col=0,header=0, parse_dates=True, squeeze=True)
obs_data = obs_data.between_time(m,a).resample('D').mean()

##################################
# Readin Models
##################################
trng_hh = pd.date_range(datetime(2006,1,1,0,0,0), end=datetime(2018,12,31,23,59,59), freq='30min') #'30min','H'

#####################################
#readin wuei
in_dir='../model_data/SPA-WUEi/'  
filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
data_wuei = pd.read_csv(filename,header=0)
data_wuei.index = trng_hh
data_wuei = data_wuei.between_time(m,a).resample('D').mean()
data_wuei[data_wuei['gpp(umol/m2/s)']<1] = np.nan

###################################
# calculate predictive accuracy
###################################

data = pd.DataFrame({'trans_wuei':data_wuei['trans(W/m2)']/28.94,'trans':obs_data['Trans(mm/day)'].values,'vpd':obs_data['VPD'],'soilm':obs_data['SWP']},index=data_wuei.index)
data_mj = data[(data.index.month>4) & (data.index.month<9)]
ptile = np.nanpercentile(data_mj['soilm'],75)
data_mj = data[data['soilm']>ptile]
mod_var = 'trans_wuei'
out_mj = cal_it_performance(data_mj,mod_var,'trans','soilm','vpd',15,1,0,1)

data_ja = data[(data.index.month>4) & (data.index.month<9)]
ptile = np.nanpercentile(data_ja['soilm'],25)
data_ja = data[data['soilm']<ptile]
mod_var = 'trans_wuei'
out_ja = cal_it_performance(data_ja,mod_var,'trans','soilm','vpd',15,1,0,1)

############### plot ##############

colors = ['grey','darkslategrey','cadetblue','lightblue','steelblue']

plt.figure(figsize=(12,8), dpi= 80)
plt.rcParams.update({'font.size': 20})

plt.subplot(1,2,1)
sizes = [float(out_mj[0] - out_mj[1]),float(out_mj[2]),float(out_mj[3]),float(out_mj[4]),float(out_mj[5])]
percentages = [float((out_mj[0]-out_mj[1])/out_mj[0]),out_mj[2]/out_mj[0],out_mj[3]/out_mj[0],out_mj[4]/out_mj[0],out_mj[5]/out_mj[0]]
percentages = [int(100*round(percentages[n],2)) for n in range(len(percentages))]
labels = ['Missing Information'+"\n"+str(percentages[0])+"%",'Unique SWP'+"\n"+str(percentages[1])+"%",'Unique VPD'+"\n"+str(percentages[2])+"%",'R'"\n" + "2%",'Synergistic'"\n" + str(percentages[4]) + "%"]
squarify.plot(sizes=sizes, label=labels, color=colors, alpha=.8)
plt.title('SWP>75th percentile')
plt.axis('off')

plt.subplot(1,2,2)
sizes = [float(out_ja[0] - out_ja[1]),float(out_ja[2]),float(out_ja[3]),float(out_ja[4]),float(out_mj[5])]
percentages = [float((out_ja[0]-out_ja[1])/out_ja[0]),out_ja[2]/out_ja[0],out_ja[3]/out_ja[0],out_ja[4]/out_ja[0],out_ja[5]/out_ja[0]]
percentages = [int(100*round(percentages[n],2)) for n in range(len(percentages))]
labels = ['Missing Information'+"\n"+str(percentages[0])+"%",'Unique SWP'+"\n"+str(percentages[1])+"%",'Unique VPD'+"\n"+str(percentages[2])+"%",'R'"\n" + "2%",'Synergistic'"\n" + str(percentages[4]) + "%"]
squarify.plot(sizes=sizes, label=labels, color=colors, alpha=.8)
plt.title('SWP<25th percentile')
plt.axis('off')

plt.show()
