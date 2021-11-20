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

##############################################
########## readin observed data ##############

startdate=datetime(2006,1,1,0,0,0); enddate = datetime(2018,12,31,23,59,59)
trng_hh = pd.date_range(startdate, enddate, freq='30min') 
in_file = '../obs_data/US-Me2_AmeriFlux_2006-01-01_2018-12-31_30min.csv'
obs = pd.read_csv(in_file,index_col=0,header=0, parse_dates=True, squeeze=True)
obs_mon = obs.between_time('9:00','12:00').resample('M').mean()
obs_mon['mon'] = obs_mon.index.month

############################
#readin WUEi
############################
in_dir='../model_data/SPA-WUEi/'   
filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
data_wuei = pd.read_csv(filename,header=0)
data_wuei.index = trng_hh
data_wuei = data_wuei.between_time('9:00','12:00').resample('M').mean()
data_wuei['mon'] = data_wuei.index.month

d13c_wuei = weighted_leaflayers(in_dir,'D13C','agr(umol.m-2.s-1)',startdate,enddate)
d13c_wuei = d13c_wuei.between_time('9:00','12:00').resample('M').mean()
d13c_wuei['mon'] = d13c_wuei.index.month
lwp_wuei = weighted_leaflayers(in_dir,'psil(MPa)','agr(umol.m-2.s-1)',startdate,enddate)
lwp_wuei = lwp_wuei.between_time('9:00','12:00').resample('M').mean()
lwp_wuei['mon'] = lwp_wuei.index.month

############################
#readin WUE
############################
in_dir='../model_data/SPA-WUE/'  
filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
data_wue = pd.read_csv(filename,header=0)
data_wue.index = trng_hh
data_wue = data_wue.between_time('9:00','12:00').resample('M').mean()
data_wue['mon'] = data_wue.index.month

d13c_wue = weighted_leaflayers(in_dir,'D13C','agr(umol.m-2.s-1)',startdate,enddate)
d13c_wue = d13c_wue.between_time('9:00','12:00').resample('M').mean()
d13c_wue['mon'] = d13c_wue.index.month
lwp_wue = weighted_leaflayers(in_dir,'psil(MPa)','agr(umol.m-2.s-1)',startdate,enddate)
lwp_wue = lwp_wue.between_time('9:00','12:00').resample('M').mean()
lwp_wue['mon'] = lwp_wue.index.month

############################
#readin BB
############################
in_dir='../model_data/SPA-BB-H/'
filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
data = pd.read_csv(filename,header=0)
data.index = trng_hh
data_bb = data.between_time('9:00','12:00').resample('M').mean()
data_bb['mon'] = data_bb.index.month

d13c_bb = weighted_leaflayers(in_dir,'D13C','agr(umol.m-2.s-1)',startdate,enddate)
d13c_bb = d13c_bb.reindex(data.index)
d13c_bb[data['gpp(umol/m2/s)']<2] = np.nan
d13c_bb = d13c_bb.between_time('9:00','12:00').resample('M').mean()
d13c_bb['mon'] = d13c_bb.index.month
lwp_bb = weighted_leaflayers(in_dir,'psil(MPa)','agr(umol.m-2.s-1)',startdate,enddate)
lwp_bb = lwp_bb.between_time('9:00','12:00').resample('M').mean()
lwp_bb['mon'] = lwp_bb.index.month

############################
#readin Medlyn
############################
in_dir='../model_data/SPA-MED-H/'  
filename = os.path.join(in_dir+'ecosystem_fluxes.csv')
data_med = pd.read_csv(filename,header=0)
data_med.index = trng_hh
data_med = data_med.between_time('9:00','12:00').resample('M').mean()
data_med['mon'] = data_med.index.month

d13c_med = weighted_leaflayers(in_dir,'D13C','agr(umol.m-2.s-1)',startdate,enddate)
d13c_med = d13c_med.between_time('9:00','12:00').resample('M').mean()
d13c_med['mon'] = d13c_med.index.month
lwp_med = weighted_leaflayers(in_dir,'psil(MPa)','agr(umol.m-2.s-1)',startdate,enddate)
lwp_med = lwp_med.between_time('9:00','12:00').resample('M').mean()
lwp_med['mon'] = lwp_med.index.month

############################
# readin Sperry
############################
in_dir='../model_data/gain-risk/'   
trng_h = pd.date_range(datetime(2006,1,1,0,0,0), datetime(2018,12,31,23,59,59), freq='H') 
filename = os.path.join(in_dir+'standard_OUTPUT_timesteps.csv')
dfc2 = pd.read_csv(filename,header=0,index_col=None,skiprows=None)
dfc2.index = trng_h
dfc2['D13C'] = 4.4 + 22.6*(dfc2['ci, Pa']/28.17)
dfc2['mon'] = dfc2.index.month

data_sperry = dfc2.between_time('9:00','12:00').resample('M').mean()

##################################################
############ plot LWP vs D13C ####################
##################################################

plt.figure(num=None, figsize=(14, 9), dpi=100, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9, wspace=0.2, hspace=0.4)
plt.rcParams.update({'font.size': 16})
ax1 = plt.subplot(1,2,1)

for ij in range(5,9):
	yvar = d13c_wuei['D13C'].loc[d13c_wuei['mon']==ij]
	xvar = lwp_wuei['psil(MPa)'].loc[lwp_wuei['mon']==ij]
	ax1.scatter(xvar,yvar,40,color='sienna',alpha = 1,label='WUEi')

for ij in range(5,9):
	yvar = d13c_wue['D13C'].loc[d13c_wue['mon']==ij]
	xvar = lwp_wue['psil(MPa)'].loc[lwp_wue['mon']==ij]
	ax1.scatter(xvar,yvar,40,color='darkgoldenrod',alpha = 1,label='WUE')

for ij in range(5,9):
	yvar = d13c_bb['D13C'].loc[d13c_bb['mon']==ij]
	xvar = lwp_bb['psil(MPa)'].loc[lwp_bb['mon']==ij]
	ax1.scatter(xvar,yvar,40,color='olivedrab',alpha = 1,label='BB-H')

for ij in range(5,9):
	yvar = d13c_med['D13C'].loc[d13c_med['mon']==ij]
	xvar = lwp_med['psil(MPa)'].loc[lwp_med['mon']==ij]
	ax1.scatter(xvar,yvar,40,color='seagreen',alpha = 1,label='MED-H')

for ij in range(5,9):
	yvar = data_sperry['D13C'].loc[data_sperry['mon']==ij]
	xvar = -data_sperry['P, Mpa'].loc[data_sperry['mon']==ij]
	ax1.scatter(xvar,yvar,40,color='darkcyan',alpha=1,label='gain-risk')

ax1.text(0.09, 0.95, '(a)', transform=ax1.transAxes,fontsize=18, va='top', ha='right')

plt.xlabel(r'$P (MPa)$',fontsize=22)
plt.ylabel(r'$D13C$',fontsize=22)
plt.ylim([10,28])
plt.xlim([-2.8,-.25])

ax2 = plt.subplot(1,2,2)

for ij in range(5,9):
	yvar = data_wuei['gpp(umol/m2/s)']/(data_wuei['trans(W/m2)']/28.94)
	xvar = lwp_wuei['psil(MPa)'].loc[lwp_wuei['mon']==ij]
	ax2.scatter(xvar,yvar.loc[data_wuei['mon']==ij],40,color='sienna',alpha = 1,label='WUEi')

for ij in range(5,9):
	yvar = data_wue['gpp(umol/m2/s)']/(data_wue['trans(W/m2)']/28.94)
	xvar = lwp_wue['psil(MPa)'].loc[lwp_wue['mon']==ij]
	ax2.scatter(xvar,yvar.loc[data_wue['mon']==ij],40,color='darkgoldenrod',alpha = 1,label='WUE')

for ij in range(5,9):
	yvar = data_bb['gpp(umol/m2/s)']/(data_bb['trans(W/m2)']/28.94)
	xvar = lwp_bb['psil(MPa)'].loc[lwp_bb['mon']==ij]
	ax2.scatter(xvar,yvar.loc[data_bb['mon']==ij],40,color='olivedrab',alpha = 1,label='SPA-BB-H')

for ij in range(5,9):
	yvar = data_med['gpp(umol/m2/s)']/(data_med['trans(W/m2)']/28.94)
	xvar = lwp_med['psil(MPa)'].loc[lwp_med['mon']==ij]
	ax2.scatter(xvar,yvar.loc[data_med['mon']==ij],40,color='seagreen',alpha = 1,label='SPA-MED-H')

for ij in range(5,9):
	yvar = (data_sperry['Anet, umol s-1m-2']*3.4)/(data_sperry['E (mm timestep-1)']*24)
	xvar = -data_sperry['P, Mpa'].loc[data_sperry['mon']==ij]
	ax2.scatter(xvar,yvar.loc[data_sperry['mon']==ij],40,color='darkcyan',alpha = 1,label='gain-risk')

#for ij in range(5,9):
#	yvar = (obs_mon['GPP_PI'])/(obs_mon['Trans(mm/day)'])
#	xvar = -data_sperry['P, Mpa'].loc[data_sperry['mon']==ij]
#	plt.scatter(xvar,yvar.loc[data_sperry['mon']==ij],80,color='k',alpha = 0.7,label='Sperry')

ax2.text(0.09, 0.95, '(b)', transform=ax2.transAxes,fontsize=18, va='top', ha='right')

plt.xlabel(r'$P (MPa)$',fontsize=22)
plt.ylabel(r'$A/T$', fontsize=22)
plt.xlim([-2.8,-.25])

#plt.legend(loc='upper left')

plt.show()
