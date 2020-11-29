#!/usr/bin/env python3
import os
import h5py
import subprocess as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from qharv.reel import scalar_dat


if __name__ == '__main__':
    
    rs_all = [25, 30]
    # time step for different rs
    tau_inv_all_all = [[0.4, 0.32, 0.256, 0.2, 0.16, 0.128, 0.1, 0.08, 0.064, 0.05, 0.032],[0.2, 0.16, 0.128, 0.1, 0.08, 0.064, 0.05, 0.032]]
    # change N and corresponding clong, cwide. a is the length of cell.
    N_all = [57, 81, 121]
    cutk_ratio_all = [10.0, 10.0, 14.0]
    for i in range(0,len(rs_all)):
      rs = rs_all[i]
      print('rs=',rs)
      tau_inv_all = tau_inv_all_all[i]
      lam = 1.0/rs/rs
      eps = 2.0/rs
      cmd = 'cd ..; mkdir rs%d' % (rs)
      sp.check_call(cmd,shell=True)
      for j in range(0,len(N_all)):
        N = N_all[j]
        print('N=',N)
        cutk_ratio = cutk_ratio_all[j]
        L = (N*np.pi)**0.5
        diagonal = 2**0.5*L
        cutk = cutk_ratio*np.pi/L
        cmd = 'cd ../rs%d; mkdir N_%d' % (rs,N)
        sp.check_call(cmd,shell=True)
        for k in range(0,len(tau_inv_all)):
          tau_inv = tau_inv_all[k]
          tau = 1.0/tau_inv
          print('tau=',tau)
          cutr = L/2.0-(lam*tau)**0.5
          cmd = 'cd ../rs%d/N_%d; mkdir tau_%.3f' % (rs, N, tau)
          sp.check_call(cmd,shell=True)
          cmd = 'cd ../rs%d/N_%d/tau_%.3f; mkdir gendm' % (rs, N, tau)
          sp.check_call(cmd,shell=True)
          path = '../rs%d/N_%d/tau_%.3f/gendm/rs%d.%d.in'
          pid = 'rs%d.%d' % (rs,N)
          f = open(path%(rs,N,tau,rs,N),'w')
          f.write('UNITS R a0\n')
          f.write('TYPE e %.8e\n' % (lam))
          f.write('TYPE e %.8e\n' % (lam))
          f.write('GRID 450 LOG .0001 %.8f\n' % (L/2.0))
          f.write('SQUARER %f 3 2 3 40 10\n' % (tau_inv))
          f.write('POT COULOPT %.8f %.6f %.6f 2 0.D0 1.0 %.7f %.7f' % (cutr, cutk, eps, L, L))
          f.close()
          cmd = 'cd ../rs%d/N_%d/tau_%.3f/gendm; potgen_lr %s > %s.plr.out' % (rs,N,tau,pid,pid)
          sp.check_call(cmd,shell=True)
          cmd = 'cd ../rs%d/N_%d/tau_%.3f/gendm; squarer %s > %s.sq.out' % (rs,N,tau,pid,pid)
          sp.check_call(cmd,shell=True)


# end __main__
