#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from qharv.reel import mole
from qharv.sieve import mean_df
from qharv.plantation import kyrt
from common_func import drummond_dv

def write_skk(fout, rs, fregex='skk.json', cache_dir='cache'):
  flist = mole.findall('rs%d*/%s'%(rs,fregex), cache_dir)
  dfl = []
  for floc in flist:
    df1 = pd.read_json(floc)
    # label
    df1['temp'] = df1['spath'].apply(lambda x: int(x.replace('t_', '')))
    df1_groupby = df1[['temp','sk0_mean_mean']].groupby(by='temp',as_index=False).max()
    df2 = pd.merge(df1_groupby,df1,on=['temp','sk0_mean_mean'],how='left')
    path = os.path.dirname(floc)
    meta = meta_from_path(path)
    for key, val in meta.items():
      df2[key] = val
    dfl.append(df2)
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

def paul_plot(rs, pol, highest_temp):
  from qharv.plantation import sugar
  # step 1: gather data
  fjson = 'cache/all-skk-rs%s.json'%(rs)
  sugar.cache(write_skk)(fjson,rs)

  # step 2: clean data
  df = pd.read_json(fjson)
  sel = (df.pol == pol) & (df.temp < highest_temp)
  df = df[sel].reset_index(drop=True)
  df = df.sort_values(['typ','temp'],ascending=False).reset_index(drop=True)
  
  # step 3: plot
  fig_loc = 'rs%d-%s-e-skk.png'%(rs,pol)
  kyrt.set_style()
  fig, ax = plt.subplots(1, 1)

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
        sel = (df.typ == typ) & (df.method == method) & (df.sl == sl)
        df1 = df[sel].reset_index(drop=True)
        if (df1.empty):
            continue

        # get right xlims
        if (method=='FP'):
            high_temp = int(max(df1.temp,default=0))
            high_temp_all.append(high_temp)

        marker = marker_map[method]
        ms = ms_map[sl]
        style = style0.copy()
        style.update({'markeredgecolor': myc, 'marker': marker, 'ms': ms})

        ax.plot(df1.temp,df1.sk0_mean_mean, **style)
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

  method_leg = kyrt.create_legend(ax, styles, labels, loc='lower left', bbox_to_anchor=(1.04,0))
  ax.add_artist(method_leg)
  highest_temp = min(highest_temp,max(high_temp_all))
  ax.set_xlim(100, highest_temp)
  ax.set_xlabel('Temperature (uRy)')
  ax.set_ylabel('S(G)')
  # touch up
  fig.tight_layout()
  fig.savefig(fig_loc, dpi=320)
  #plt.show()


if __name__ == '__main__':
  rs_all = [25, 30, 40, 50]
  pol_all = ['pol', 'unp']
  for rs in rs_all:
    for pol in pol_all:
      print(rs,pol)
      paul_plot(rs,pol,700)
# end __main__
