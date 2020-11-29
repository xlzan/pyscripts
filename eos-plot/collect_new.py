#!/usr/bin/env python
import os 
import numpy as np
import pandas as pd

from qharv.reel import mole
from qharv.plantation import sugar

def read_folder(folder, obs, nskip=0):
  fregex = '*.s*.%s.json' % obs
  flist = mole.findall(fregex, folder)
  stl = [fname.split('.')[-3] for fname in flist]
  seriesl = [int(st.replace('s', '')) for st in stl]
  idx = np.argsort(seriesl)
  dfl = []
  for i in idx[nskip:]:
      df1 = pd.read_json(flist[i])
      dfl.append(df1)
  df = pd.concat(dfl, sort=False).reset_index(drop=True)
  return df

@sugar.cache
def collect_ld(fout, folders, suf, ld, **kwargs):
  mylabels = ld['labels']
  mydrops = ld['drops']
  dfl = []
  for folder in folders:
    df1 = read_folder(folder, suf, **kwargs)
    dfl.append(df1)
  df = pd.concat(dfl, sort=False)
  # averge over node, core
  df['spath'] = df.path.apply(lambda x:'/'.join(x.split('/')[:-2]))
  # divide columns into labels, data, and drops
  labels0 = ['spath', 'series']
  drops0 = ['node', 'core']
  labels = labels0+mylabels
  drops = drops0+mydrops
  # groupby labels, average data, drop drops
  mdf = df.groupby(labels).mean().drop(columns=drops)
  edf = df.groupby(labels).std(ddof=1).apply(lambda x: x/len(x)**0.5)
  # combine mean and error
  def sufcols(df, suf):
      ymap = {}
      for col in df.columns:
          ymap[col] = col+suf
      return ymap
  mdf.rename(columns=sufcols(mdf, '_mean'), inplace=True)
  edf.rename(columns=sufcols(edf, '_error'), inplace=True)
  medf = pd.concat([mdf, edf], axis=1, sort=False)
  medf.reset_index().to_json(fout)

def edf_func_err(df, labels, ecols):
  edf = df.groupby(labels)[ecols].apply(
          lambda x: sum(x**2)**0.5/len(x))
  return edf

@sugar.cache
def collect_dat(fout, folders, suf, **kwargs):
  fproc = kwargs.pop('fproc', None)
  if fproc is None:  # no intermediate processing
    fproc = lambda x: x
  dfl = []
  for folder in folders:
    df1 = read_folder(folder, suf, **kwargs)
    dfl.append(df1)
  df = pd.concat(dfl, sort=False)
  fproc(df)
  # averge over node, core
  df['spath'] = df.path.apply(lambda x:'/'.join(x.split('/')[:-2]))
  # divide columns into labels, data, and drops
  labels = ['spath', 'series']
  drops = ['node', 'core']
  cols = np.setdiff1d(df.columns, labels+drops+['path'])
  mcols = [col for col in cols if col.endswith('_mean')]
  ecols = [col for col in cols if col.endswith('_error')]
  assert len(mcols) == len(ecols)
  # groupby labels, average data, drop drops
  mdf = df.groupby(labels)[mcols].mean()
  edf = df.groupby(labels)[ecols].apply(lambda x: (x**2).sum()**0.5/len(x))
  ## why not use standard error?
  #sdf = df.groupby(labels)[mcols].apply(lambda x: np.std(x, ddof=1)/len(x)**0.5)
  #col_map = {}
  #for mc, ec in zip(mcols, ecols):
  #  col_map[mc] = ec
  #edf = sdf.rename(columns=col_map)
  medf = pd.concat([mdf, edf], axis=1, sort=False)
  medf.reset_index().to_json(fout)

def ratio_sc_df(df):
  ncols = [col for col in df.columns if col.startswith('n_')]
  for nc in ncols:
    if nc.endswith('error'): continue
    name = nc[2:].replace('mean', '').strip('_')
    ymn = '%s_mean' % name
    yen = '%s_error' % name
    norms = df[nc].values
    df[ymn] /= norms
    df[yen] /= norms
  df.drop(columns=ncols, inplace=True)

def collect_grsk(fout, folders, suf):
  from qharv.sieve import mean_df
  drops = ['node', 'core']
  dfl = []
  for folder in folders:
    df = read_folder(folder, suf)
    # drop 'node', 'core' from path
    df['spath'] = df.path.apply(lambda x:'/'.join(x.split('/')[:-2]))
    df.drop(columns='path', inplace=True)
    # average
    labels, mcols, ecols = mean_df.categorize_columns(df)
    for drop in drops:
      labels.remove(drop)
    cols = [mc.replace('_mean', '') for mc in mcols]
    mdf = df.groupby(labels).apply(mean_df.dfme, cols).reset_index()
    dfl.append(mdf)
  mydf = pd.concat(dfl).reset_index(drop=True)
  mydf.to_json(fout)

def collect_grsk_new(fout, folders, suf):
  from qharv.sieve import mean_df
  drops = ['node', 'core']
  dfl = []
  for folder in folders:
    df = read_folder(folder, suf)
    # drop 'node', 'core' from path
    df['spath'] = df.path.apply(lambda x:'/'.join(x.split('/')[:-2]))
    df.drop(columns='path', inplace=True)
    # average
    labels, mcols, ecols = mean_df.categorize_columns(df)
    for drop in drops:
      labels.remove(drop)
    mdfm = df.groupby(labels).mean().add_suffix('_mean').reset_index(drop=False)
    mdfe = df.groupby(labels).std().add_suffix('_error').reset_index(drop=False)
    mdf = pd.merge(mdfm,mdfe,how='right',on=labels)
    mdf.drop(["core_mean","core_error","node_mean","node_error"],axis=1,inplace=True)
    dfl.append(mdf)
  mydf = pd.concat(dfl).reset_index(drop=True)
  mydf.to_json(fout)

if __name__ == '__main__':
  folder = './data/'
  rs_all = os.listdir(folder)
  typ = 'NOLEAK'
  for rs in rs_all:
      print(rs)
      #rs25-pol_cryst-WCGS
      rsp,suffix,node=rs.split('-')
      if (suffix=='pol_liq'):
          N=57
      elif (suffix=='unp_liq'):
          N=58
      else:
          N=56
      flist = mole.findall('t_*', folder+rs, ftype='d')
      #rs25-pol_cryst-56-FP-nodal
      fout = 'cache/%s-%s-%d-%s-%s/'%(rsp,suffix,N,node,typ)
      for suf, fproc in zip(['dat', 'sc'], [None, None]):
        jout = fout+'%s.json' % suf
        collect_dat(jout, flist, suf, fproc=fproc)
      for suf in ['grr', 'skk']:
        jout = fout+'%s.json' % suf
        collect_grsk_new(jout, flist, suf)
# end __main__
