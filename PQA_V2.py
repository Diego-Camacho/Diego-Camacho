###########################################
#####     PQC score  computation      #####
#####    Versión 2.0, 20/01/2020      #####
##### Diego Camacho - Victor Nieto    #####
###########################################



import numpy as np
import matplotlib.pyplot as plt


def pear_cor(signal_der,signal_izq ):
    #Generation of correlated matrix : corr_mat
    corr_mat = np.corrcoef(signal_der,signal_izq)
    # Retorno del coeficiente de pearson [0,1]
    return corr_mat[0, 1]

def mean_auto(perf):
    # Arry that will hold 1000 autocorrelations
    auto_corr= np.empty(1000)
    for i in range(1000):
      shuffled= np.random.choice(perf,size=len(perf))
      auto_corr[i]= pear_cor(shuffled[1:],shuffled[:-1])
      return np.mean(auto_corr)

def PQC_score(VP):
    # Generation of the set vector, in perfect order
    l = [[x] *VP.count(x) for x in set(VP)]
    perf_VP = np.array([item for sublist in l for item in sublist])
    #  PQC score
    VP = np.array(VP)
    PQC= ( pear_cor(VP[1:],VP[:-1]) - mean_auto(perf_VP) )/ pear_cor(perf_VP[1:],perf_VP[:-1])
    return PQC

def Z_score(VP):
    
    PQCx=PQC_score(VP)
    rand_PQC= np.empty(1000)
    
    for i in range(1000):
        shuffled= np.random.choice(VP,size=len(VP), replace=False)
        shuffled= shuffled.tolist()
        rand_PQC[i]=PQC_score(shuffled)
        
    Zx= (PQCx - np.mean(rand_PQC))/np.std(rand_PQC)    
    return Zx

###############################################
##Intrinsec noise due to clasification number##
###############################################
## Noise proportions

import math

def Noise_proportions(VP):

    l = [[x] *VP.count(x) for x in set(VP)]
    sort_V = np.array([item for sublist in l for item in sublist])
    
    noise_p=[]
    #Randomized by shuffling different percentages of items
    v = [str(i) for i in sort_V]
    v = np.array(v)

    for x in range(0, 110, 10):
        n_items = (x/100.0)*len(v)
        index = np.random.choice(range(0,len(v)), size=int(n_items), replace=False)
        swapping = np.random.choice(index, size=len(index), replace=False)
    
        v_swap=[x for x in v]
        for y in range(0,int(n_items)):
            v_swap[index[y]], = v[swapping[y]]

        noise_p.append(Z_score([int(i) for i in v_swap]))

    return noise_p

def perf_clasnumb(classes=10,length=100):
  elements=int(round(float(length)/classes))
  if float(length)/classes < elements:
    rest= length/classes - 1
  elif  float(length)/classes == elements:
    rest=0
  else:
    rest= length/classes + 1
  print(rest)
  l=[[x]*elements for x in range(classes-1)]
  if rest:
    l.append(list([classes-1]*rest))
  perf=[item for sublist in l for item in sublist]
  print(perf)
  Zs=Z_score(perf)
  return(Zs)

