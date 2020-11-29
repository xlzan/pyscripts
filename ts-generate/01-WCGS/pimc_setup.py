#!/usr/bin/env python
import os
import numpy as np
import subprocess as sp


if __name__ == '__main__':
   
    #pol is polarization of the system: pol --> polarized case; unp --> unpolarized case
    pol = 'pol'
    #initial configuration: s --> solid; l --> liquid
    initial = 's'
    # node name: w --> WCGS node; f --> Free Particle node; g --> Ground State node
    node_nm = 'w'
    # wcgs ratio from Drummond's PRL2009, 0.071 --> pol; 0.082 --> unp 
    wcgs_ratio = 0.071
    # NREF value: w --> 2; f --> 2; g --> 0
    in_nref = 2
    # leakage: 1 --> nodal action; 0 --> NOLEAK
    leakage = 1
    # specify rs
    rs_all = [25, 30]
    tau_inv_all_all = [[0.4, 0.32, 0.256, 0.2, 0.16, 0.128, 0.1, 0.08, 0.064, 0.05, 0.032],[0.2, 0.16, 0.128, 0.1, 0.08, 0.064, 0.05, 0.032]]
    # gamma is the 3-cycle permutation prefactor
    gamma = 2.0
    # specify temperature
    temp_all_all = [[200], [200]]
    # level_shift to calculate level for very (rs,temp,tau) point
    # level = int(max(np.floor(np.log2(Nslices)-level_shift),1))
    level_shift_all = [5, 5]
    # nstep per block = nsteps_all/N*2**(4-level)
    nsteps_all = [420000, 420000]
    ## npp = npp_ori/(2**(max(np.ceil(np.log2(Nslices)),9)-9))
    #npp_ori = 16
    # npp_all_all for [56, 80 ,120] based on memory/upi, [12,24]: nslices>512 --> npp=12; nslices<=512 --> npp=24
    npp_all_all = [[12,24], [8,18], [10,10]]
    # change N and corresponding clong, cwide.
    N_all = [56, 80, 120]
    # ncell_type is the cell type: 1 --> simple cube(square); 3 --> triangular lattice
    ncell_type = 3
    # ncell_displace is the displacement percent in crys_coord: 0 --> perfect lattice; 1 --> liquid
    ncell_displace = 0.3
    ncell_long = [4, 5, 6]
    ncell_wide = [7, 8, 10]
    a = (2.0*np.pi/3**0.5)**0.5
    cutk_ratio_all = [10.0, 10.0, 14.0]
    # path % (rs,N,tau,temp,initial,node_nm)
    path = '../rs%d/N_%d/tau_%.3f/t_%d/%s_%s.sy'

    for i in range(0,len(rs_all)):
     rs = rs_all[i]
     print('rs=',rs)
     tau_inv_all = tau_inv_all_all[i]
     temp_all = temp_all_all[i]
     level_shift = level_shift_all[i]
     nsteps = nsteps_all[i]
     lamb_e = 1.0/rs**2
     for j in range(0,len(N_all)):
       N = N_all[j]
       print('N=',N)
       clong = ncell_long[j]
       cwide = ncell_wide[j]
       cutk_ratio = cutk_ratio_all[j]
       npp_all = npp_all_all[j]
       L1 = clong*3**0.5*a
       L2 = cwide*a
       L = min(L1,L2)
       cutk = cutk_ratio*np.pi/L
       for k in range(0,len(tau_inv_all)):
         tau_inv = tau_inv_all[k]
         tau = 1.0/tau_inv
         print('tau=',tau)
         cmd = 'cd ../rs%d/N_%d/tau_%.3f; mkdir conf' % (rs,N,tau)
         sp.check_call(cmd,shell=True)
         cmd = '''cd ../rs%d/N_%d/tau_%.3f/conf; 
                  echo -e %s%d'\n' %d'\n' 1.0'\n' 2'\n' %d %d %d'\n' %.1f'\n' | crys_coord''' % (rs,N,tau,initial,N,N,ncell_type,clong,cwide,ncell_displace)
         sp.check_call(cmd,shell=True)
         for m in range(0,len(temp_all)):
           temp = temp_all[m]
           print('temp=',temp)
           cmd = 'cd ../rs%d/N_%d/tau_%.3f; mkdir t_%d' % (rs,N,tau,temp)
           sp.check_call(cmd,shell=True)
           cmd = 'cd ../rs%d/N_%d/tau_%.3f/t_%d; ln -s ../gendm/*.dm ../gendm/*.yk ./' % (rs,N,tau,temp)
           sp.check_call(cmd,shell=True)
           cmd = 'cd ../rs%d/N_%d/tau_%.3f/t_%d; ln -s ../conf/*.ic ./' % (rs,N,tau,temp)
           sp.check_call(cmd,shell=True)
           f = open(path%(rs,N,tau,temp,initial,node_nm),'w')
           f.write('UNITS R a0\n')
           f.write('ENORM %.4e\n' % (N))
           f.write('BOXSIZE %.7f %.7f\n' % (L1, L2))
           beta = 1e6/temp
           nslices_ori = beta/tau
           nslices = round(beta/tau)
           if np.mod(nslices,2)!=0:
              flag = nslices_ori-nslices
              if flag > 0:
                  nslices = nslices+1
              else:
                  nslices = nslices-1
           beta_now = tau*nslices
           f.write('BETA %.7e\n' % (beta_now))
           f.write('NSLICES %d\n' % (nslices))
           if (pol == 'pol'):
               f.write('TYPE e %.6e 1 %d %d %s%d.e.ic\n' % (lamb_e, N, N, initial, N))
           else:
               f.write('TYPE e %.6e 2 %d %d %d %s%d.e.ic\n' % (lamb_e, N, N/2.0, N/2.0, initial, N))
           f.write('CUTK %.8e 1\n' % (cutk))
           f.write('POT PAIR e e rs%s.%d.dm 10 10\n' % (rs,N))
           f.write('NREF %d\n' %(in_nref))
           if (not leakage):
               f.write('NOLEAK\n')
           if (node_nm == 'w'):
               f.write('BESTGAUSSIAN %.3f %d 0.5\n' %(wcgs_ratio,rs))
           f.write('FREESAMP e\n')
           f.write('VWINDOW 8 16\n')
           if cutk <6.0 :
              f.write('ADDKRING %.6f 6.0\n' % (cutk*1.01))
           f.write('GAMMA 1 e 1.0\n')
           f.write('GAMMA 2 e 0.0\n')
           f.write('GAMMA 3 e %.1f\n' % gamma)
           f.write('TIMER 800 432000\n')
           f.write('RESTART \n')
           level = int(max(np.floor(np.log2(nslices)-level_shift),1))
           nstep = nsteps/N*2**(4-level)
           f.write('OMOVE   140   %d  %d    \n' % (nstep,level))
           f.write('SELECT  200   %d  %d  %d\n' % (int(nstep/2.0),level,2*N))
           f.write('SELECT  200   %d  %d  %d\n' % (int(nstep/2.0),level,2*N))
           f.close()
           #npp = npp_ori/(2**(max(np.ceil(np.log2(nslices)),9)-9))
           if (nslices > 512):
               npp = npp_all[0]
           else:
               npp = npp_all[1]
           cmd = 'cd ../rs%d/N_%d/tau_%.3f; paraset -n %d t_%d -e' % (rs, N, tau, npp, temp)
           sp.check_call(cmd,shell=True)


# end __main__
