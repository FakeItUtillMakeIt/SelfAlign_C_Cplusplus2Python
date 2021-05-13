import numpy as np 
import os 
import json
from collections import defaultdict
# from Bio import SeqIO
# import pandas as pd 
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import shape 
#from ont_fast5_api.fast5_interface import get_fast5_file
from scipy.io import loadmat,savemat
#from statsmodels import robust
import h5py 

from align import double_array,double_p
import align 
from ctypes import *
import subprocess
import time



def smatrixACBF(x1,K1,x2,K2,pfbad = -1000000,selfalignment = True,lookback = 10):
    # pfbad = -1000000
    # selfalignment = True
    # lookback = 10
    D1 = np.array([np.diag(x) for x in K1]).T
    D2 = np.array([np.diag(x) for x in K2]).T
    s=-1*np.array([np.inf]*x1.shape[1]*(lookback+1)).reshape(x1.shape[1],lookback+1)
    
    def eig_(x):
        if x.shape[0] == 1:
            return [x]
        else:
            return np.linalg.eig(x)

    def inv_(x):
        if x.shape[0] == 1:
            return 1/x 
        else:
            return np.linalg.inv(x)

    if selfalignment:
        for c1 in range(x1.shape[1]):
            for c2 in range(max(0,c1-lookback),c1+1):
                dx_sq = np.power((x1[:,c1]-x2[:,c2]),2)
                V = 1/D1[:,c1] + 1/D2[:,c2]
                good_components = pfbad < -0.5*( np.log(2*np.pi*V) + dx_sq/V )
                this_K1 = K1[c1][good_components,good_components]
                this_K2 = K2[c2][good_components,good_components]
                if sum(good_components) == 2:
                    this_K1=np.diag(this_K1)
                    this_K2=np.diag(this_K2)
                this_x1 = x1[good_components, c1]
                this_x2 = x2[good_components, c2]
                num_good = sum(good_components)
                Kx = np.dot(this_K1,this_x1) + np.dot(this_K2,this_x2)
                s[c1,c2-c1+lookback] = (x1.shape[0]-num_good)*pfbad + 0.5*(-num_good*np.log(2*np.pi) + sum(np.log(eig_(this_K1)[0])) + sum(np.log(eig_(this_K2)[0])) - sum(np.log(eig_(this_K1+this_K2)[0])) - np.dot(np.dot(this_x1.T,this_K1),this_x1) - np.dot(np.dot(this_x2.T,this_K2),this_x2) + np.dot(np.dot(Kx.T,inv_((this_K1 + this_K2))),Kx))
    return s
#remove toggle
read = h5py.File('/mnt/raid/lla/hectobio/basecalling/read_1_section_1.mat','r')
x_f = read['read/x_f'][()].T
stiffs = read['read/k_f'][()].T[0]
k_f = [read[stiffs[i]][()] for i in range(len(stiffs))]
npts_f = read['read/npts_f'][()].T[0]
error_floor = [0,100]
psa = -4
lookback = 10
#step_probabilities = [step | skip | skip+ | back | back+ | hold | bad | bad+ | optional p_sa/p_slip]
step_probabilities =  [0.6, 1e-2, 1e-2, .27, .2, .13, 0, 0,np.nan]
cov_floor = np.diag(np.power(error_floor,2))
filteredFeatures = x_f
filteredStiffs = k_f
filteredNumPts = npts_f
total_self_alignment = np.array(range(filteredFeatures.shape[1]))#1:size(filteredFeatures,2)
selfalignment = np.ones(filteredFeatures.shape[1])
stiffs_array = []
for cL in range(len(filteredStiffs)):
    covmat = np.linalg.inv(filteredStiffs[cL])
    covmat = covmat + cov_floor
    stiffs_array.append(np.linalg.inv(covmat))
x_meas = filteredFeatures
K_meas = stiffs_array
x_ref=[]
K_ref=[]
step_probabilities =  [0.6, 1e-2, 1e-2, .27, .2, .13, 0, 0, np.nan]
psa = 4
lookback = 10
step_penalties = [-0.2528,  
                  -2.0928,
                  -2.0928,
                  -3.8128,
                  -5.9128,
                  -2.9128,
                  -5.1328,
                  -1.8496,
                  -8.5100
                  ]
              
pfbad = -1000000*np.ones(x_meas.shape[0])
use_periodic_boundaries = False
modes = 1
slip_locations = np.zeros(x_meas.shape[1])
sa_lookback = min(lookback, x_meas.shape[1]-1)
alignment_type = 'S'
x_ref = x_meas
K_ref = K_meas
lookback = sa_lookback
step_penalties[-1] = psa
reorder = False
isdiag = False
isdc = False
isprincomp = False
score_matrix = smatrixACBF(x_meas,K_meas,x_meas,K_meas)

step_penalties = np.log(step_probabilities)
step_penalties[7] = step_penalties[6]
if alignment_type == 'S':
    step_penalties[8] = psa
step_penalties=step_penalties.reshape(1,-1)
if step_penalties.shape[1] == 9:
    if step_penalties.shape[0] == modes:
        step_penalties = np.tile(step_penalties, (1, score_matrix.shape[0]))
        step_penalties = step_penalties[0]

print(score_matrix.shape)


##调用C
##行小列大没问题，行大列小出现内存错误

def alignmentC(score_matrix,step_penalties,modes,lookback,alignment_type,slip_locations,use_periodic_boundaries):
    num_row = c_int(score_matrix.shape[0])
    num_col = c_int(score_matrix.shape[1])
    score_matrix2=double_array(num_row.value*num_col.value)
    for i in range(len(score_matrix.flatten('F'))):
        #print(score_matrix.flatten('F')[i])
        score_matrix2[i]=score_matrix.flatten('F')[i]



    modes = double_array(1)
    modes[0] = 1.0
    num_modes  = c_int(1)
    lookback=c_int(10)

    step_penalties2=double_array(step_penalties.shape[0])
    for i in range(len(step_penalties)):
        step_penalties2[i] = step_penalties[i]

    slip_locations2=double_array(slip_locations.shape[0])
    for i in range(len(slip_locations)):
        slip_locations2[i] = slip_locations[i]
    use_periodic_boundaries = c_bool()
    use_periodic_boundaries = False
    alignment_type = c_char()
    alignment_type = 'S'

    alignment_trace=double_array(num_row.value)
    bad_levels=double_array(num_row.value)
    cumulate_score=double_array(num_row.value)

    alignment_matrix=double_array(num_row.value*num_col.value)

    start_time=time.time()
    align.mexFunction(7,alignment_trace,bad_levels,cumulate_score,alignment_matrix,score_matrix2,num_row.value,num_col.value,step_penalties2,modes,num_modes.value,lookback.value,alignment_type,slip_locations2,use_periodic_boundaries)

    list_align_trace=[]
    list_bad_levels=[]
    list_cumulate_score=[]
    list_align_mat=[]

    for i in range(num_row.value):
        list_align_trace.append(alignment_trace[i])

    for i in range(num_row.value):
        list_bad_levels.append(bad_levels[i])
    
    for i in range(num_row.value):
        list_cumulate_score.append(cumulate_score[i])

    for i in range(num_row.value*num_col.value):
        list_align_mat.append(alignment_matrix[i])

    array_align_trace=np.array(list_align_trace)
    array_bad_levels=np.array(list_bad_levels)
    array_cumulate_score=np.array(list_cumulate_score)
    array_align_mat=np.array(list_align_mat)
    # print(array_align_mat.reshape(num_col.value,num_row.value).T)
    # print(time.time()-start_time)


    return array_align_trace,array_bad_levels,array_cumulate_score,array_align_mat.reshape(num_col.value,num_row.value).T

trace,levels,cumulate,align_mat=alignmentC(score_matrix,step_penalties,modes,lookback,alignment_type,slip_locations,use_periodic_boundaries)

print(trace,levels,cumulate)
