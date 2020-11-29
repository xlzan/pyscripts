#!/usr/bin/env python
import os
import numpy as np
import subprocess as sp


if __name__ == '__main__':
   
    #initial configuration: s --> solid; l --> liquid
    initial = 'l'
    # node name: w --> WCGS node; f --> Free Particle node; g --> Ground State node
    node_nm = 'g'
    # specify rs
    rs_all = [25, 30]
    tau_inv_all_all = [[0.4, 0.32, 0.256, 0.2, 0.16, 0.128, 0.1, 0.08, 0.064, 0.05, 0.032],[0.2, 0.16, 0.128, 0.1, 0.08, 0.064, 0.05, 0.032]]
    # specify temperature
    temp_all_all = [[200,500], [200]]
    # npp_all_all for [56, 80 ,120] based on memory/upi, [12,24]: nslices>512 --> npp=12; nslices<=512 --> npp=24
    npp_all_all = [[12,24], [8,18], [10,10]]
    # change N and corresponding clong, cwide.
    N_all = [57, 81, 121]
    NN_all = [58, 81, 122]
    # path % (rs,N,tau,temp,initial,node_nm)
    path = '../rs%d/N_%d/submit_shell/tau_%.3f-t_%d.sh'
    BIN_path = '/vol8/home/nudtyuanjianmin/zxl/software/upi13-2DEG/upi13-N%d-%d/upi2d'

    for i in range(0,len(rs_all)):
     rs = rs_all[i]
     print('rs=',rs)
     tau_inv_all = tau_inv_all_all[i]
     temp_all = temp_all_all[i]
     for j in range(0,len(N_all)):
       N = N_all[j]
       NN = NN_all[j]
       print('N=',N)
       npp_all = npp_all_all[j]
       cmd = 'cd ../rs%d/N_%d/; mkdir submit_shell' % (rs,N)
       sp.check_call(cmd,shell=True)
       for k in range(0,len(tau_inv_all)):
         tau_inv = tau_inv_all[k]
         tau = 1.0/tau_inv
         print('tau=',tau)
         for m in range(0,len(temp_all)):
           temp = temp_all[m]
           print('temp=',temp)
           beta = 1e6/temp
           nslices_ori = beta/tau
           nslices = round(beta/tau)
           if np.mod(nslices,2)!=0:
              flag = nslices_ori-nslices
              if flag > 0:
                  nslices = nslices+1
              else:
                  nslices = nslices-1
           if (nslices > 512):
               npp = npp_all[0]
           else:
               npp = npp_all[1]
           BIN = BIN_path%(NN,max(nslices,512))
           f = open(path%(rs,N,tau,temp),'w')
           f.write('#!/bin/bash\n')
           f.write('cwd=`pwd`\n')
           f.write('echo "number of processors: %d"\n' % npp)
           f.write('BIN=%s\n' % BIN)
           f.write('prefix=%s_%s.s000  # start this run\n' % (initial,node_nm))
           f.write('date\n')
           f.write('NODE_DIR=../tau_%.3f/t_%d/node0 \n' % (tau,temp))
           f.write('for (( icopy=0; icopy<%d; icopy++   )); do\n' % npp)
           f.write('    cdir=$NODE_DIR/core$icopy\n')
           f.write('    cd $cdir\n')
           f.write('    echo $prefix | $BIN >& $prefix.err &\n')
           f.write('    sleep 1\n')
           f.write('    cd $cwd\n')
           f.write('done\n')
           f.write('wait\n')
           f.write('date\n')
           f.close()
           cmd = 'cd ../rs%d/N_%d/submit_shell; chmod +x tau_%.3f-t_%d.sh' % (rs, N, tau, temp)
           sp.check_call(cmd,shell=True)





# end __main__
