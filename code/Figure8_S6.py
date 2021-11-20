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
from it_metrics import *


st = datetime(2006,1,1,0,0,0)
en = datetime(2018,12,31,23,59,59)

m = '8:00'
a = '16:00'

##################################
# Readin FLUXnet data
##################################

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

#readin wue
in_dir='../model_data/SPA-WUE/'  
filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
data_wue = pd.read_csv(filename,header=0)
data_wue.index = trng_hh
data_wue = data_wue.between_time(m,a).resample('D').mean()
data_wue[data_wue['gpp(umol/m2/s)']<1] = np.nan

#readin bb
in_dir='../model_data/SPA-BB-H/'  
filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
data_bb = pd.read_csv(filename,header=0)
data_bb.index = trng_hh
data_bb = data_bb.between_time(m,a).resample('D').mean()
data_bb[data_bb['gpp(umol/m2/s)']<1] = np.nan

#readin medlyn
in_dir='../model_data/SPA-MED-H/'  
filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
data_med = pd.read_csv(filename,header=0)
data_med.index = trng_hh
data_med = data_med.between_time(m,a).resample('D').mean()
data_med[data_med['gpp(umol/m2/s)']<1] = np.nan

# readin sperry
in_dir='../model_data/gain-risk/'   
trng_h = pd.date_range(datetime(2006,1,1,0,0,0), datetime(2018,12,31,23,59,59), freq='H') 
filename = os.path.join(in_dir+'standard_OUTPUT_timesteps.csv')
dfc2 = pd.read_csv(filename,header=0,index_col=None,skiprows=None)
dfc2.index = trng_h

data_sperry = dfc2.between_time(m,a).resample('D').mean()

###################################
# calculate predictive accuracy
###################################

out_data = np.empty([5,7])

# build dataframe 
data = pd.DataFrame({'trans_wue':data_wue['trans(W/m2)']/28.94,'trans_wuei':data_wuei['trans(W/m2)']/28.94,'trans_bb':data_bb['trans(W/m2)']/28.94,'trans_med':data_med['trans(W/m2)']/28.94,'trans_sperry':data_sperry['E (mm timestep-1)']*24,'trans':obs_data['Trans(mm/day)'].values,'vpd':obs_data['VPD'],'soilm':obs_data['SWP']},index=data_wue.index)
data_mjja = data[(data.index.month>4) & (data.index.month<9)]
ptile25 = np.percentile(data_mjja['soilm'],25)
ptile75 = np.percentile(data_mjja['soilm'],75)
#data = data_mjja[data_mjja['soilm']<ptile25]
data = data_mjja[data_mjja['soilm']>ptile75]

mod_var = ['trans_wuei','trans_wue','trans_bb','trans_med','trans_sperry']

for mod in range(len(mod_var)):
	out = cal_it_performance(data,mod_var[mod],'trans','soilm','vpd',15,1,0)
	out_data[mod,:] = out

#np.savetxt('out_data_swp_JJ.csv',out_data,fmt='%4.4f',delimiter=',')

#################################
# bootstrap
#################################

wuei_it = np.empty([10000,7])
ss = int(len(data)*0.8)
for ij in range(10000):
	subset = random.sample(range(len(data)),ss)
	data_subset = data.iloc[subset]
	wuei_it[ij,:] = cal_it_performance(data_subset,'trans_wuei','trans','soilm','vpd',15,1,0)


wue_it = np.empty([10000,7])
for ij in range(10000):
	subset = random.sample(range(len(data)),ss)
	data_subset = data.iloc[subset]
	wue_it[ij,:] = cal_it_performance(data_subset,'trans_wue','trans','soilm','vpd',15,1,0)


bb_it = np.empty([10000,7])
for ij in range(10000):
	subset = random.sample(range(len(data)),ss)
	data_subset = data.iloc[subset]
	bb_it[ij,:] = cal_it_performance(data_subset,'trans_bb','trans','soilm','vpd',15,1,0)


med_it = np.empty([10000,7])
for ij in range(10000):
	subset = random.sample(range(len(data)),ss)
	data_subset = data.iloc[subset]
	med_it[ij,:] = cal_it_performance(data_subset,'trans_med','trans','soilm','vpd',15,1,0)


sperry_it = np.empty([10000,7])
for ij in range(10000):
	subset = random.sample(range(len(data)),ss)
	data_subset = data.iloc[subset]
	sperry_it[ij,:] = cal_it_performance(data_subset,'trans_sperry','trans','soilm','vpd',15,1,0)

######### Figures ##################

colors=['sienna','goldenrod','olivedrab','seagreen','darkcyan']
model_ap=np.array([wuei_it[:,0],wue_it[:,0],bb_it[:,0],med_it[:,0],sperry_it[:,0]])
model_aft=np.array([wuei_it[:,5],wue_it[:,5],bb_it[:,5],med_it[:,5],sperry_it[:,5]])
model_afp=np.array([wuei_it[:,6],wue_it[:,6],bb_it[:,6],med_it[:,6],sperry_it[:,6]])
positions = np.array([1,2,3,4,5])
mod_labels=['WUEi','WUE','BB-H','MED-H','Gain-Risk']

plt.figure(num=None, figsize=(13, 11), dpi=100, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9, wspace=0.23, hspace=0.45)
plt.rcParams.update({'font.size': 14})
plt.suptitle("SWP > 75th percentile", fontsize=16)
ax = plt.subplot(4,2,1)
for i,c in enumerate(colors):
	ax.boxplot([model_ap[i,:], ],positions = [positions[i], ],
	showfliers=False,widths=0.5,patch_artist=True,
	boxprops=dict(facecolor=c, color=c),
    medianprops=dict(color='w', lw=0.5),
    capprops=dict(color=c, lw=2), whiskerprops=dict(color=c, lw=2))
ax.text(0.08, 0.95, '(a)', transform=ax.transAxes,fontsize=14, va='top', ha='right')
ax.set_ylabel(r'$A_{P}$', fontsize=14)
ax.set_xlim([0.5, len(colors) + 0.5])
ax.set_ylim([0.6, 1])
ax.set_xticks(positions)
ax.set_xticklabels(mod_labels, fontsize=14)
ax.set_title('predictive performance', fontsize=14)

ax = plt.subplot(4,2,3)
for i,c in enumerate(colors):
	ax.boxplot([model_aft[i,:], ],positions = [positions[i], ],
	showfliers=False,widths=0.5,patch_artist=True,
    boxprops=dict(facecolor=c, color=c),
    medianprops=dict(color='w', lw=0.5),
    capprops=dict(color=c,lw=2), whiskerprops=dict(color=c, lw=2))

ax.plot([0,len(colors) + 0.5],[0,0],linestyle='--',color='darkgrey')
ax.text(0.08, 0.95, '(b)', transform=ax.transAxes,fontsize=14, va='top', ha='right')
ax.set_ylabel(r'$A_{f,t}$', fontsize=14)
ax.set_xlim([0.5, len(colors) + 0.5])
ax.set_ylim([-0.6, 0.6])
ax.set_xticks(positions)
ax.set_xticklabels(mod_labels, fontsize=14)
ax.set_title('functional performance', fontsize=14)

ax = plt.subplot(4,2,4)
for i,c in enumerate(colors):
	ax.boxplot([model_afp[i,:], ],positions = [positions[i], ],
	showfliers=False,widths=0.5,patch_artist=True,
    boxprops=dict(facecolor=c, color=c),
    medianprops=dict(color='w', lw=0.5),
    capprops=dict(color=c,lw=2), whiskerprops=dict(color=c, lw=2))

ax.plot([0,len(colors) + 0.5],[0,0],linestyle='--',color='darkgrey')
ax.text(0.08, 0.95, '(c)', transform=ax.transAxes,fontsize=14, va='top', ha='right')
ax.set_ylabel(r'$A_{f,p}$', fontsize=14)
ax.set_xlim([0.5, len(colors) + 0.5])
ax.set_ylim([-0.05, 1.05])
ax.set_xticks(positions)
ax.set_xticklabels(mod_labels, fontsize=14)
ax.set_title('functional accuracy', fontsize=14)

####### unique information ##################

model_afu1=np.array([wuei_it[:,1],wue_it[:,1],bb_it[:,1],med_it[:,1],sperry_it[:,1]])
model_afu2=np.array([wuei_it[:,2],wue_it[:,2],bb_it[:,2],med_it[:,2],sperry_it[:,2]])

ax = plt.subplot(4,2,5)
for i,c in enumerate(colors):
	ax.boxplot([model_afu1[i,:], ],positions = [positions[i], ],
	showfliers=False,widths=0.5,patch_artist=True,
	boxprops=dict(facecolor=c, color=c),
    medianprops=dict(color='w', lw=0.5),
    capprops=dict(color=c, lw=2), whiskerprops=dict(color=c, lw=2))

ax.plot([0,len(colors) + 0.5],[0,0],linestyle='--',color='darkgrey')
ax.text(0.08, 0.95, '(d)', transform=ax.transAxes,fontsize=14, va='top', ha='right')
ax.set_ylabel(r'$A_{f,u,swp}$', fontsize=14)
ax.set_xlim([0.5, len(colors) + 0.5])
ax.set_ylim([-0.35, 0.35])
ax.set_xticks(positions)
ax.set_xticklabels(mod_labels, fontsize=14)
ax.set_title('unique-SWP',fontsize='14')

ax = plt.subplot(4,2,6)
for i,c in enumerate(colors):
	ax.boxplot([model_afu2[i,:], ],positions = [positions[i], ],
	showfliers=False,widths=0.5,patch_artist=True,
    boxprops=dict(facecolor=c, color=c),
    medianprops=dict(color='w', lw=0.5),
    capprops=dict(color=c,lw=2), whiskerprops=dict(color=c, lw=2))

ax.plot([0,len(colors) + 0.5],[0,0],linestyle='--',color='darkgrey')
ax.text(0.08, 0.95, '(e)', transform=ax.transAxes,fontsize=14, va='top', ha='right')
ax.set_ylabel(r'$A_{f,u,vpd}$', fontsize=14)
ax.set_xlim([0.5, len(colors) + 0.5])
ax.set_ylim([-0.35, 0.35])
ax.set_xticks(positions)
ax.set_xticklabels(mod_labels, fontsize=14)
ax.set_title('unique-VPD',fontsize='14')


####### synergystic & redundant information ##################

model_afs=np.array([wuei_it[:,3],wue_it[:,3],bb_it[:,3],med_it[:,3],sperry_it[:,3]])
model_afr=np.array([wuei_it[:,4],wue_it[:,4],bb_it[:,4],med_it[:,4],sperry_it[:,4]])

ax = plt.subplot(4,2,7)
for i,c in enumerate(colors):
	ax.boxplot([model_afs[i,:], ],positions = [positions[i], ],
	showfliers=False,widths=0.5,patch_artist=True,
	boxprops=dict(facecolor=c, color=c),
    medianprops=dict(color='w', lw=0.5),
    capprops=dict(color=c, lw=2), whiskerprops=dict(color=c, lw=2))
ax.plot([0,len(colors) + 0.5],[0,0],linestyle='--',color='darkgrey')
ax.text(0.08, 0.95, '(f)', transform=ax.transAxes,fontsize=14, va='top', ha='right')
ax.set_ylabel(r'$A_{f,s}$', fontsize=14)
ax.set_xlim([0.5, len(colors) + 0.5])
ax.set_ylim([-0.35, 0.35])
ax.set_xticks(positions)
ax.set_xticklabels(mod_labels, fontsize=14)
ax.set_title('synergistic',fontsize=14)

ax = plt.subplot(4,2,8)
for i,c in enumerate(colors):
	ax.boxplot([model_afr[i,:], ],positions = [positions[i], ],
	showfliers=False,widths=0.5,patch_artist=True,
    boxprops=dict(facecolor=c, color=c),
    medianprops=dict(color='w', lw=0.5),
    capprops=dict(color=c,lw=2), whiskerprops=dict(color=c, lw=2))
ax.plot([0,len(colors) + 0.5],[0,0],linestyle='--',color='darkgrey')
ax.text(0.08, 0.95, '(g)', transform=ax.transAxes,fontsize=14, va='top', ha='right')
ax.set_ylabel(r'$A_{f,r}$', fontsize=14)
ax.set_xlim([0.5, len(colors) + 0.5])
ax.set_ylim([-0.35, 0.35])
ax.set_xticks(positions)
ax.set_xticklabels(mod_labels, fontsize=14)
ax.set_title('redundant',fontsize=14)

plt.savefig("Figure_S6.png", dpi=300)

plt.show()
