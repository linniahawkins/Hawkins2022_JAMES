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

m = '4:00'; a = '19:00'

##############################################
########## readin observed data ##############
stdate=datetime(2006,1,1,0,0,0); endate = datetime(2018,12,31,23,59,59)
trng_hh = pd.date_range(stdate, endate, freq='30min')
in_file = '../obs_data/US-Me2_AmeriFlux_2006-01-01_2018-12-31_30min.csv'
obs_data = pd.read_csv(in_file,index_col=0,header=0, parse_dates=True, squeeze=True)
                            
#############################################
###### calculate canopy conductance Gc ######
vpd = obs_data['VPD'].values
Ta = obs_data['TA'].values
T = obs_data['Trans(mm/day)'].values/3.4 # divide by LAI (m2 leaf area)
T = T/86400*1000 # convert mm/day to g/m2/s
Kg = 115.8 + 0.4236*Ta # (KPa m3 / Kg)

Gs = Kg*(T/vpd)*40.1 # convert mm/s to mmol/m2/s
Gs_lo = Kg*(T*0.6/vpd)*40.1
Gs_hi = Kg*(T*1.4/vpd)*40.1

df_Gs = pd.DataFrame({'Gs': Gs,'Gs_lo':Gs_lo,'Gs_hi':Gs_hi},index=obs_data.index)
df_Gs = df_Gs.between_time(m,a)

##############################################
################# readin WUEi ################
in_dir='../model_data/SPA-WUEi/' 
startdate = datetime(2006,1,1,0,0,0) 
enddate = datetime(2018,12,31,23,59,59)

tr_canopy = readin_leaflayers(in_dir,'etr(W.m-2)',startdate,1)
wuei_tr = tr_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

gs_canopy = weighted_leaflayers(in_dir,'gs(mmol.m-2.s-1)','agr(umol.m-2.s-1)',startdate,enddate)
wuei_gs = gs_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

ci_canopy = weighted_leaflayers(in_dir,'ci(ppm)','agr(umol.m-2.s-1)',startdate,enddate)
wuei_ci = ci_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

ca_canopy = weighted_leaflayers(in_dir,'coa(ppm)','agr(umol.m-2.s-1)',startdate,enddate)
wuei_ca = ca_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

lwp_canopy = weighted_leaflayers(in_dir,'psil(MPa)','agr(umol.m-2.s-1)',startdate,enddate)
wuei_lwp = lwp_canopy.between_time('5:30','19:30').reindex(trng_hh, fill_value=np.nan)

agr_canopy = readin_leaflayers(in_dir,'agr(umol.m-2.s-1)',startdate,1)
wuei_agr = agr_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)


##############################################
################# readin WUE #################
in_dir='../model_data/SPA-WUE/' 
filename = os.path.join(in_dir+'layer_sun_02.csv')

tr_canopy = readin_leaflayers(in_dir,'etr(W.m-2)',startdate,1)
wue_tr = tr_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

gs_canopy = weighted_leaflayers(in_dir,'gs(mmol.m-2.s-1)','agr(umol.m-2.s-1)',startdate,enddate)
wue_gs = gs_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

ci_canopy = weighted_leaflayers(in_dir,'ci(ppm)','agr(umol.m-2.s-1)',startdate,enddate)
wue_ci = ci_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

ca_canopy = weighted_leaflayers(in_dir,'coa(ppm)','agr(umol.m-2.s-1)',startdate,enddate)
wue_ca = ca_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

lwp_canopy = weighted_leaflayers(in_dir,'psil(MPa)','agr(umol.m-2.s-1)',startdate,enddate)
wue_lwp = lwp_canopy.between_time('5:30','19:30').reindex(trng_hh, fill_value=np.nan)

agr_canopy = readin_leaflayers(in_dir,'agr(umol.m-2.s-1)',startdate,1)
wue_agr = agr_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

###############################################
############### readin BB-H ###################
in_dir='../model_data/SPA-BB-H/'

tr_canopy = readin_leaflayers(in_dir,'etr(W.m-2)',startdate,1)
bb_tr = tr_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

gs_canopy = weighted_leaflayers(in_dir,'gs(mmol.m-2.s-1)','agr(umol.m-2.s-1)',startdate,enddate)
bb_gs = gs_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

ci_canopy = weighted_leaflayers(in_dir,'ci(ppm)','agr(umol.m-2.s-1)',startdate,enddate)
bb_ci = ci_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

ca_canopy = weighted_leaflayers(in_dir,'coa(ppm)','agr(umol.m-2.s-1)',startdate,enddate)
bb_ca = ca_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

lwp_canopy = weighted_leaflayers(in_dir,'psil(MPa)','agr(umol.m-2.s-1)',startdate,enddate)
bb_lwp = lwp_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

agr_canopy = readin_leaflayers(in_dir,'agr(umol.m-2.s-1)',startdate,1)
bb_agr = agr_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

###############################################
################ readin Medlyn ################
in_dir='../model_data/SPA-MED-H/'  

tr_canopy = readin_leaflayers(in_dir,'etr(W.m-2)',startdate,1)
med_tr = tr_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

gs_canopy = weighted_leaflayers(in_dir,'gs(mmol.m-2.s-1)','agr(umol.m-2.s-1)',startdate,enddate)
med_gs = gs_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

ci_canopy = weighted_leaflayers(in_dir,'ci(ppm)','agr(umol.m-2.s-1)',startdate,enddate)
med_ci = ci_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

ca_canopy = weighted_leaflayers(in_dir,'coa(ppm)','agr(umol.m-2.s-1)',startdate,enddate)
med_ca = ca_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

lwp_canopy = weighted_leaflayers(in_dir,'psil(MPa)','agr(umol.m-2.s-1)',startdate,enddate)
med_lwp = lwp_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

agr_canopy = readin_leaflayers(in_dir,'agr(umol.m-2.s-1)',startdate,1)
med_agr = agr_canopy.between_time(m,a).reindex(trng_hh, fill_value=np.nan)

###############################################
################# readin sperry ###############
in_file = '../model_data/gain-risk/standard_OUTPUT_timesteps.csv'
index_h = pd.date_range(datetime(2006,1,1,0,0,0), datetime(2018,12,31,23,59,59), freq='H') 
data_sperry = pd.read_csv(in_file,header=0,index_col=None,skiprows=None)
data_sperry.index = index_h
data_sperry = data_sperry.between_time(m,a).reindex(index_h, fill_value=np.nan)

###############################################
##################### plot ####################
plt.style.use('seaborn-pastel')

st = datetime(2010,8,10,0,0,0)
end = datetime(2010,8,15,0,0,0)

plt.figure(num=None, figsize=(11, 14), dpi=80, facecolor='w', edgecolor='k')
ax1 =plt.subplot(6,1,1)
ax1.plot(obs_data['SWP'][st:end],color='navy',linestyle='dotted',label='SWP')
ax1.set_ylabel(r'$SWP (MPa)$',color='navy',fontsize=14)
ax1.tick_params(axis='y', colors='navy',labelsize=14)
ax1.tick_params(axis='x',labelsize=14)
ax2 = ax1.twinx()
ax2.plot(obs_data['VPD'][st:end],color='sienna',linestyle='dotted',label='VPD (kPa)')
ax2.set_ylabel(r'$VPD (kPa)$',color='sienna',fontsize=14)
ax2.tick_params(axis='y', colors='sienna',labelsize=14)
ax2.text(0.03, 0.95, '(a)', transform=ax1.transAxes,fontsize=12, va='top', ha='right')
ax2.get_xaxis().set_ticks([])
plt.legend(loc='upper right')
plt.tick_params(labelsize=14)

ax1=plt.subplot(6,1,2)
ax1.fill_between(obs_data[st:end].index,obs_data['Trans(mm/day)'][st:end].values*0.6, obs_data['Trans(mm/day)'][st:end].values*1.4, facecolor='grey', alpha=0.3)
ax1.plot(obs_data['Trans(mm/day)'][st:end],color='k',linestyle='dotted',label='Observed')
ax1.plot(wuei_tr[st:end]/28.94,color='sienna',label='WUEi') # convert W/m2 to mm/day 
ax1.plot(wue_tr[st:end]/28.94,color='darkgoldenrod',label='WUE',linestyle='--')
ax1.plot(bb_tr[st:end]/28.94,color='olivedrab',label='BB-H')
ax1.plot(med_tr[st:end]/28.94,color='seagreen',label='MED-H',linestyle='--')
ax1.plot(data_sperry['E (mm timestep-1)'][st:end]*24,color='darkcyan',label='Gain-Risk')
ax1.text(0.03, 0.95, '(b)', transform=ax1.transAxes,fontsize=12, va='top', ha='right')
ax1.set_ylabel(r'$T (mm/day)$',fontsize=14)
ax1.tick_params(labelsize=14)
ax1.set_xticklabels([])
plt.legend(loc='upper right',framealpha=1)

ax1 = plt.subplot(6,1,3)
ax1.fill_between(df_Gs[st:end].index,df_Gs['Gs_lo'][st:end].values, df_Gs['Gs_hi'][st:end].values, facecolor='grey', alpha=0.3)
ax1.plot(df_Gs['Gs'][st:end],color='k',linestyle='dotted',label='Observed')
ax1.plot(wuei_gs[st:end],color='sienna',label='WUEi')
ax1.plot(wue_gs[st:end],color='darkgoldenrod',label='WUE',linestyle='--')
ax1.plot(bb_gs[st:end],color='olivedrab',label='Ball-Berry')
ax1.plot(med_gs[st:end],color='seagreen',label='Medlyn',linestyle='--')
ax1.plot(data_sperry['Gw, mmol m-2s-1'][st:end],color='darkcyan',label='Sperry')
ax1.text(0.03, 0.95, '(c)', transform=ax1.transAxes,fontsize=12, va='top', ha='right')
ax1.set_ylabel(r'$G_{c} (mmol/{m^2}/s)$',fontsize=14)
#ax1.set_ylim([0,65])
ax1.tick_params(labelsize=14)
ax1.get_xaxis().set_ticks([])

ax1 = plt.subplot(6,1,4)
ax1.plot(wuei_lwp[st:end],color='sienna',label='WUEi')
ax1.plot(wue_lwp[st:end],color='darkgoldenrod',label='WUE',linestyle='--')
ax1.plot(bb_lwp[st:end],color='olivedrab',label='Ball-Berry')
ax1.plot(med_lwp[st:end],color='seagreen',label='Medlyn',linestyle='--')
plt.plot(-data_sperry['P, Mpa'][st:end],color='darkcyan',label='Sperry')
ax1.text(0.03, 0.95, '(d)', transform=ax1.transAxes,fontsize=12, va='top', ha='right')
plt.ylabel(r'$P (MPa)$',fontsize=14)
ax1.tick_params(labelsize=14)
ax1.get_xaxis().set_ticks([])

ax1 = plt.subplot(6,1,5)
ax1.fill_between(obs_data[st:end].index,obs_data['GPP_PI'][st:end].values*0.8, obs_data['GPP_PI'][st:end].values*1.2, facecolor='grey', alpha=0.3)
ax1.plot(obs_data['GPP_PI'][st:end],color='k',linestyle='dotted',label='Observed')
ax1.plot(wuei_agr[st:end],color='sienna',label='WUEi')
ax1.plot(wue_agr[st:end],color='darkgoldenrod',label='WUE',linestyle='--')
ax1.plot(bb_agr[st:end],color='olivedrab',label='Ball-Berry')
ax1.plot(med_agr[st:end],color='seagreen',label='Medlyn',linestyle='--')
plt.plot(data_sperry['Anet, umol s-1m-2'][st:end]*3.2,color='darkcyan',label='Sperry')
ax1.text(0.03, 0.95, '(e)', transform=ax1.transAxes,fontsize=12, va='top', ha='right')
plt.ylabel(r'$GPP (umol/{m^2}/s)$',fontsize=14)
ax1.tick_params(labelsize=14)
ax1.get_xaxis().set_ticks([])

ax1 = plt.subplot(6,1,6)
ax1.plot(wuei_ci['ci(ppm)'][st:end]/wuei_ca['coa(ppm)'][st:end],color='sienna',label='WUEi')
ax1.plot(wue_ci['ci(ppm)'][st:end]/wue_ca['coa(ppm)'][st:end],color='darkgoldenrod',label='WUE',linestyle='--')
ax1.plot(bb_ci['ci(ppm)'][st:end]/bb_ca['coa(ppm)'][st:end],color='olivedrab',label='Ball-Berry')
ax1.plot(med_ci['ci(ppm)'][st:end]/med_ca['coa(ppm)'][st:end],color='seagreen',label='Medlyn',linestyle='--')
ax1.plot(data_sperry['ci, Pa'][st:end]/28.17,color='darkcyan',label='Sperry')
ax1.text(0.03, 0.95, '(f)', transform=ax1.transAxes,fontsize=12, va='top', ha='right')
ax1.set_ylabel(r'$C_{i}/C_{a}$',fontsize=14)
ax1.tick_params(labelsize=14)


plt.savefig('0.diurnal_leaflayers.png', dpi=300)
plt.show()
