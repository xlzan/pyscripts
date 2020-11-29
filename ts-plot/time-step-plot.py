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
    df1['tau'] = df1['spath'].apply(lambda x: float(x.split('/')[0].replace('tau_', '')))
    path = os.path.dirname(floc)
    meta = meta_from_path(path)
    for key, val in meta.items():
      df1[key] = val
    dfl.append(df1)
  df = pd.concat(dfl, sort=False).reset_index(drop=True)
  df.to_json(fout)

def meta_from_path(path):
  folder = os.path.basename(path)
  # e.g. rs20-t250-pol-N56-WCGS-nodal
  rst, tt, pol, nt, method, typ = folder.split('-')
  rs = int(rst.replace('rs', ''))
  t = int(tt.replace('t', ''))
  npart = int(nt.replace('N', ''))
  meta = {'rs': rs, 'temp': t, 'pol': pol, 'npart': npart, 'method': method, 'typ': typ}
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


def paul_plot(rs,etype,temp):
  from qharv.plantation import sugar
  # step 1: gather data
  fjson = 'cache/all-ts.json'
  sugar.cache(write_sc)(fjson)

  # step 2: clean data
  df = pd.read_json(fjson)
  # step 2.1: add DMC data
  #dfl = get_DMC()
  #dfl.reset_index(drop=True, inplace=True) 
  #df = pd.concat([df,dfl],sort=True)
  #df.reset_index(drop=True, inplace=True) 
  # remove entries with rs
  #df = df.query('rs==%d and pol==%s and temp<%d'%(rs,pol,highest_temp))
  sel = (df.rs == rs) & (df.temp == temp)
  df = df[sel].reset_index(drop=True)
  if (df.empty):
      return
  if (etype == 'E_tot' or etype == 'E_virial'):
      prefactor = 2.0
  else:
      prefactor = 1.0
  df = df.sort_values(['npart','method'],ascending=True).reset_index(drop=True)
  
  # step 3: plot
  fig_loc = 'rs%d-t%d-%s.png'%(rs,temp,etype)
  kyrt.set_style()
  fig, ax = plt.subplots(1, 1)
  ymult = 1e6  # Ry -> uRy

  # assign one color to each N (56, 80, 120 and 168)
  nparts = df.npart.unique()
  colors = kyrt.dark8
  # assign one marker to each method
  methods = df.method.unique()
  marker_map = {'WCGS': 's', 'GS': 'o', 'FP': 'x'}
  # assign marker sizes to different type
  typs = df.typ.unique()
  ms_map = {'nodal': 7, 'NOLEAK': 10}
  mew = 1
  style0 = {'markeredgecolor': 'k', 'c': 'none', 'marker': 's', 'ms': 10, 'mew': mew, 'linestyle': ''}
  for npart, myc in zip(nparts, colors):
    for method in methods:
       for typ in typs:

        sel = (df.typ == typ) & (df.method == method) & (df.npart == npart)
        dft = df[sel].reset_index(drop=True)
        if (dft.empty):
            continue
        # use drummond's fsc formula (2008)
        x, np = mean_df.xyye(df, 'tau', 'npart', sel=sel, yerr=False, sort=True)
        dy = drummond_dv(np,rs)
        x, ym, ye = mean_df.xyye(df, 'tau', etype, sel=sel, sort=True)

        marker = marker_map[method]
        ms = ms_map[typ]
        style = style0.copy()
        style.update({'markeredgecolor': myc, 'marker': marker, 'ms': ms})

        ax.errorbar(x, (ym+prefactor*dy)*ymult, ye*ymult, **style)

  # add legends
  styles = []
  for npart, myc in zip(nparts, colors):
    style = style0.copy()
    style.update({'markeredgecolor': myc})
    styles.append(style)
  npart_leg = kyrt.create_legend(ax, styles, nparts, loc='upper left', bbox_to_anchor=(1.04,1))
  ax.add_artist(npart_leg)
  styles = []
  labels = []
  for method in methods:
    for typ in typs:
      sel = (df.method == method) & (df.typ == typ)
      dft = df[sel].reset_index(drop=True)
      if (dft.empty):
          continue
      style = style0.copy()
      label = method+' '+typ
      labels.append(label)
      marker = marker_map[method]
      ms = ms_map[typ]
      style.update({'marker': marker, 'c': 'none', 'ms': ms})
      styles.append(style)
  method_leg = kyrt.create_legend(ax, styles, labels, loc='lower left', bbox_to_anchor=(1.04,0))
  ax.add_artist(method_leg)

  ax.set_xlabel('Tau (Ry-1)')
  ax.set_ylabel('Energy (uRy/e)')
  # touch up
  fig.tight_layout()
  fig.savefig(fig_loc, dpi=320)
  plt.cla()
  plt.close("all")
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
  #rs_all = [25, 30, 35, 40, 45, 50]
  rs_all = [45]
  #etype_all = ['E_tot', 'PE_tot', 'KE_tot']
  etype_all = ['E_virial']
  #temp_all = [125, 250]
  temp_all = [125]
  for rs in rs_all:
   for etype in etype_all:
    for temp in temp_all:
      print(rs,etype,temp)
      paul_plot(rs,etype,temp)
# end __main__
