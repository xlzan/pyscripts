#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from qharv.reel import mole
from qharv.sieve import mean_df
from qharv.plantation import kyrt
from common_func import drummond_dv

def write_sc(fout, fregex='dat.json', cache_dir='cache'):
  flist = mole.findall(fregex, cache_dir)
  dfl = []
  for floc in flist:
    df1 = pd.read_json(floc)
    # label
    df1['temp'] = df1['spath'].apply(lambda x: int(x.replace('t_', '')))
    path = os.path.dirname(floc)
    meta = meta_from_path(path)
    for key, val in meta.items():
      df1[key] = val
    dfl.append(df1)
  df = pd.concat(dfl, sort=False).reset_index(drop=True)
  df.to_json(fout)

def meta_from_path(path):
  folder = os.path.basename(path)
  # e.g. rs25-pol_cryst-56-WCGS-nodal
  rsp, cns = folder.split('_')
  rst, pol = rsp.split('-')
  sl, nt, method, typ = cns.split('-')
  if method == 'FS':  # why FS? sounds like Finite-Size
    method = 'SF'  # SignFul (SF)
  rs = int(rst.replace('rs', ''))
  npart = int(nt)
  meta = {'rs': rs, 'pol': pol, 'sl': sl, 'npart': npart, 'method': method, 'typ': typ}
  return meta

def get_DMC():
  # get DMC data from Drummond's PRL (2009).
  DMC_drummond = pd.DataFrame({
      'rs' : [25, 30, 35, 40, 45, 50, 75, 100],
      'temp' : [0, 0, 0, 0, 0, 0, 0, 0],
      'pol_cryst' : [-3.7731, -3.1917, -2.7669, -2.4432, -2.1881, -1.9817, -1.34905, -1.0243],
      'pol_cryst_error' : [0.0002, 0.0002, 0.0001, 0.0001, 0.0001, 0.0002, 0.0001, 0.0001],
      'unp_cryst' : [-3.7751, -3.1922, -2.7672, -2.4431, -2.1875, -1.9814, -1.34905, -1.0243],
      'unp_cryst_error' : [0.0003, 0.0002, 0.0001, 0.0002, 0.0002, 0.0002, 0.0001, 0.0001],
      'unp_liq' : [-3.7774, -3.1926, -2.7665, -2.4416, np.nan,np.nan, -1.3466, np.nan],
      'unp_liq_error' : [0.0002, 0.0001, 0.0001, 0.0001, np.nan,np.nan, 0.00015, np.nan],
      'pol_liq' : [-3.7740, -3.1913, -2.7657, -2.4416, np.nan,np.nan,np.nan,np.nan],
      'pol_liq_error' : [0.0002, 0.0001, 0.0001, 0.0001, np.nan,np.nan,np.nan,np.nan],
      })
  pol_all = ['pol', 'unp']
  # Drummond used WCGS nodes to calculate crystal DMC energies and GS nodes for liquid.
  sl_all = ['cryst', 'liq']
  method_all = ['WCGS', 'GS']
  typ = 'nodal'
  flag = 0
  for pol in pol_all:
      for sl, method in zip(sl_all,method_all):
          df1 = DMC_drummond.copy()[['rs','temp','%s_%s'%(pol,sl),'%s_%s_error'%(pol,sl)]]
          df1.loc[:,'pol'] = pol
          df1.loc[:,'sl'] = sl
          # Drummond's results had already been fsc-ed, large number of particle is used here
          df1.loc[:,'npart'] = 1e8
          df1.loc[:,'method'] = method
          df1.loc[:,'typ'] = typ
          df1 = df1.rename(columns={'%s_%s'%(pol,sl):'E_tot_mean','%s_%s_error'%(pol,sl):'E_tot_error'})
          #transfer to Ry unit.
          df1['E_tot_mean'] = df1['E_tot_mean']*2*1e-2
          df1['E_tot_error'] = df1['E_tot_error']*2*1e-2
          if (flag==0):
              dfl = df1.copy()
              flag = 1
          else:
              dfl = pd.concat([dfl,df1])
  return dfl 


def paul_plot(rs, pol, highest_temp):
  from qharv.plantation import sugar
  # step 1: gather data
  fjson = 'cache/all-eos.json'
  sugar.cache(write_sc)(fjson)

  # step 2: clean data
  df = pd.read_json(fjson)
  # step 2.1: add DMC data
  dfl = get_DMC()
  dfl.reset_index(drop=True, inplace=True) 

  df = pd.concat([df,dfl],sort=True)
  df.reset_index(drop=True, inplace=True) 
  # remove entries with Ground State node
  #df = df.query('rs==%d and pol==%s and temp<%d'%(rs,pol,highest_temp))
  sel = (df.rs == rs) & (df.pol == pol) & (df.temp < highest_temp)
  df = df[sel].reset_index(drop=True)
  df = df.sort_values(['typ','temp'],ascending=False).reset_index(drop=True)
  
  # step 3: plot
  fig_loc = 'rs%d-%s-e-temp.png'%(rs,pol)
  kyrt.set_style()
  fig, ax = plt.subplots(1, 1)
  ymult = 1e6  # Ry -> uRy

  # assign one color to each type (nodal or NOLEAK)
  typs = df.typ.unique()
  colors = kyrt.dark8
  # assign one marker to each method and solidity
  methods = df.method.unique()
  marker_map = {'FP': 'x', 'WCGS': 'o', 'GS': 's'}
  sls = df.sl.unique()
  ms_map = {'cryst': 7, 'liq': 10}
  mew = 1
  ls = ''
  style0 = {'markeredgecolor': 'k', 'c': 'none', 'marker': 's', 'ms': 10, 'mew': mew, 'ls': ls}
  high_temp_all = []
  for typ, myc in zip(typs, colors):
    for method in methods:
      for sl in sls:

        # all GS results are in liquid phase.
        #if (method=='GS'):
        #    sl = 'liq'
        # all WCGS results are in solid phase.
        #if ((method=='WCGS') or (method=='NOLEAK')):
        #    sl = 'cryst'
        sel = (df.typ == typ) & (df.method == method) & (df.sl == sl)
        dft = df[sel].reset_index(drop=True)
        if (dft.empty):
            continue
        # use drummond's fsc formula (2008)
        x, np = mean_df.xyye(df, 'temp', 'npart', sel=sel, yerr=False, sort=True)
        dy = drummond_dv(np,rs)
        x, ym, ye = mean_df.xyye(df, 'temp', 'E_tot', sel=sel, sort=True)
        # get right xlims
        if (method=='FP'):
            high_temp = max(x,default=0)
            high_temp_all.append(high_temp)



        marker = marker_map[method]
        ms = ms_map[sl]
        style = style0.copy()
        style.update({'markeredgecolor': myc, 'marker': marker, 'ms': ms})

        ax.errorbar(x, (ym+2*dy)*ymult, ye*ymult, **style)

  # add legends
  styles = []
  for typ, myc in zip(typs, colors):
    style = style0.copy()
    style.update({'markeredgecolor': myc})
    styles.append(style)
  typ_leg = kyrt.create_legend(ax, styles, typs, loc='upper left', bbox_to_anchor=(1.04,1))
  ax.add_artist(typ_leg)
  styles = []
  labels = []
  typ='nodal'
  for method in methods:
    for sl in sls:
      # skip following cases cause they don't exist
      #if ((method=='GS' and sl=='cryst') or (method=='WCGS' and sl=='liq') or (method=='NOLEAK' and sl=='liq')):
      #  continue
      sel = (df.typ == typ) & (df.method == method) & (df.sl == sl)
      dft = df[sel].reset_index(drop=True)
      if (dft.empty):
          continue
      style = style0.copy()
      label = method+' '+sl
      labels.append(label)
      marker = marker_map[method]
      ms = ms_map[sl]
      style.update({'marker': marker, 'c': 'none', 'ms': ms})
      styles.append(style)
  #method_leg = kyrt.create_legend(ax, styles, labels, loc='lower right')
  method_leg = kyrt.create_legend(ax, styles, labels, loc='lower left', bbox_to_anchor=(1.04,0))
  ax.add_artist(method_leg)

  highest_temp = min(highest_temp,max(high_temp_all))
  ax.set_xlim(-10, highest_temp)
  ax.set_xlabel('Temperature (uRy)')
  ax.set_ylabel('Total Energy (uRy/e)')
  # touch up
  fig.tight_layout()
  fig.savefig(fig_loc, dpi=320)
  #plt.show()

def paul_excel():
  import xlwt
  from tempfile import TemporaryFile
  from qharv.plantation import sugar
  # step 1: gather data
  fjson = 'data/all-eos.json'
  sugar.cache(write_sc)(fjson)

  # step 2: clean data
  df = pd.read_json(fjson)
  # step 2.1: add DMC data
  dfl = get_DMC()
  df = pd.concat([df,dfl],sort=True)
  # remove entries with too large errorbar
  bsel = df.E_tot_error > 100*1e-6
  df.drop(df.loc[bsel].index, inplace=True)
  df = df.sort_values(['rs','temp'],ascending=True).reset_index(drop=True)
  
  # step 3: create excel
  file_loc = 'eos_table.xlsx'
  writer = pd.ExcelWriter(file_loc)

  pols = df.pol.unique()
  methods = df.method.unique()
  sls = df.sl.unique()
  for pol in pols:
    for method in methods:
      for sl in sls:
        sel = (df.pol == pol) & (df.method == method) & (df.sl == sl)
        dft = df[sel].reset_index(drop=True)
        # use drummond's fsc formula (2008)
        dft['fsc'] = 2*drummond_dv(dft.npart,dft.rs)
        dft['E_tot_mean'] = dft['E_tot_mean']+dft['fsc']
        data = dft[['rs','temp','E_tot_mean','E_tot_error']]
        data.to_excel(writer,'%s_%s_%s'%(pol,sl,method))
  writer.save()


if __name__ == '__main__':
  #paul_excel()
  rs_all = [25, 30, 40, 50]
  pol_all = ['pol', 'unp']
  #rs_all = [50]
  for rs in rs_all:
    for pol in pol_all:
      print(rs,pol)
      paul_plot(rs,pol,1000)
# end __main__
