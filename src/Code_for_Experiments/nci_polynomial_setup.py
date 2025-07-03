from dataclasses import replace
import numpy as np
import random
import matplotlib.pyplot as plt
import networkx as nx
from math import log, ceil
import pandas as pd
import seaborn as sns
import nci_linear_setup as ncls
from scipy import interpolate, special
from itertools import combinations
import scipy.sparse

# Scale down the effects of higher order terms
a1 = 0.8      # for linear effects
a2 = 0.2    # for quadratic effects
a3 = 1   # for cubic effects
a4 = 1   # for quartic effects

# Define f(z)
f_linear = lambda alpha, z, gz: alpha + a1*z
f_quadratic = lambda alpha, z, gz: alpha + a1*z + a2*np.multiply(gz,gz)
f_cubic = lambda alpha, z, gz: alpha + a1*z + a2*np.multiply(gz,gz) + a3*np.power(gz,3)
f_quartic = lambda alpha, z, gz: alpha + a1*z + a2*np.multiply(gz,gz) + a3*np.power(gz,3) + a4*np.power(gz,4)

# Define f(z)_cov_adj
f_linear_cov_adj = lambda alpha, z, gz, X, b: alpha + a1*z + X.dot(b)
#f_quadratic_cov_adj = lambda alpha, z, gz, X, b: alpha + a1*z + a2*np.multiply(gz,gz) + X.dot(b)
f_quadratic_cov_adj = lambda alpha, z, gz, X, b, C_linear, C_quad: alpha + a1*C_linear.dot(z) + a2*np.multiply(gz,gz) + X.dot(b) - a2*(C_quad**2).dot(z) / (np.array(np.sum(C_quad,1)).flatten())**2
f_cubic_cov_adj = lambda alpha, z, gz, X, b: alpha + a1*z + a2*np.multiply(gz,gz) + a3*np.power(gz,3) + X.dot(b)
f_quartic_cov_adj = lambda alpha, z, gz, X, b: alpha + a1*z + a2*np.multiply(gz,gz) + a3*np.power(gz,3) + a4*np.power(gz,4) + X.dot(b)

def ppom(beta, C, alpha):
  '''
  Returns k-degree polynomial potential outcomes function fy
  
  f (function): must be of the form f(z) = alpha + z + a2*z^2 + a3*z^3 + ... + ak*z^k
  C (np.array): weighted adjacency matrix
  alpha (np.array): vector of null effects
  '''
  # n = C.shape[0]
  # assert np.all(f(alpha, np.zeros(n), np.zeros(n)) == alpha), 'f(0) should equal alpha'
  #assert np.all(np.around(f(alpha, np.ones(n)) - alpha - np.ones(n), 10) >= 0), 'f must include linear component'

  if beta == 0:
      return lambda z: alpha + a1*z
  elif beta == 1:
      f = f_linear
      return lambda z: alpha + a1*C.dot(z)
  else:
      g = lambda z : C.dot(z) / np.array(np.sum(C,1)).flatten()
      if beta == 2:
          f = f_quadratic
      elif beta == 3:
          f = f_cubic
      elif beta == 4:
          f = f_quartic
      else:
          print("ERROR: invalid degree")
      return lambda z: f(alpha, C.dot(z), g(z)) 

def ppom_cov_adj(beta, C_quad, alpha, X, b, C_linear):
  if beta == 0:
      return lambda z: alpha + a1*z
  elif beta == 1:
      f = f_linear_cov_adj
      return lambda z: alpha + a1*C_linear.dot(z)
  else:
      g = lambda z : C_quad.dot(z) / np.array(np.sum(C_quad,1)).flatten()
      if beta == 2:
          f = f_quadratic_cov_adj
      elif beta == 3:
          f = f_cubic_cov_adj
      elif beta == 4:
          f = f_quartic_cov_adj
      else:
          print("ERROR: invalid degree")
      return lambda z: f(alpha, z, g(z), X, b, C_linear, C_quad)



SNIPE_beta = lambda n,y,w : np.sum(y*w)/n 

def SNIPE_weights(n, p, A, z, beta):
  '''
  Compute the weights w_i(z) for each population unit

  n (int): population size
  p (float): treatment probability
  A (scipy csr array): adjacency matrix in scipy csr format
  z (numpy array): treatment vector
  beta (int): degree of the potential outcomes model

  Returns a numpy array W where the i-th entry is the weight w_i(z) associated to unit i
  '''
  treated_neighb = A.dot(z)
  control_neighb = A.dot(1-z)
  W = np.zeros(n)
  for i in range(n):
    w = 0
    a_lim = min(beta,int(treated_neighb[i]))
    for a in range(a_lim+1):
      b_lim = min(beta - a,int(control_neighb[i]))
      for b in range(b_lim+1):
        w = w + ((1-p)**(a+b) - (-p)**(a+b)) * p**(-a) * (p-1)**(-b) * special.binom(treated_neighb[i],a)  * special.binom(control_neighb[i],b)
    W[i] = w

  return W



def poly_interp_splines(n, P, sums, spltyp = 'quadratic'):
  '''
  Returns estimate of TTE using spline polynomial interpolation 
  via scipy.interpolate.interp1d

  n (int): popluation size
  P (numpy array): sequence of probabilities p_t
  sums (numpy array): sums of outcomes at each time step
  spltyp (str): type of spline, can be 'quadratic, or 'cubic'
  '''
  assert spltyp in ['quadratic', 'cubic'], "spltyp must be 'quadratic', or 'cubic'"
  f_spl = interpolate.interp1d(P, sums, kind=spltyp, fill_value='extrapolate')
  TTE_hat = (1/n)*(f_spl(1) - f_spl(0))
  return TTE_hat

def poly_interp_linear(n, P, sums):
  '''
  Returns two estimates of TTE using linear polynomial interpolation 
  via scipy.interpolate.interp1d
  - the first is with kind = 'linear' (as in... ?)
  - the second is with kind = 'slinear' (as in linear spline)

  n (int): popluation size
  P (numpy array): sequence of probabilities p_t
  sums (numpy array): sums of outcomes at each time step
  '''

  #f_lin = interpolate.interp1d(P, sums, fill_value='extrapolate')
  f_spl = interpolate.interp1d(P, sums, kind='slinear', fill_value='extrapolate')
  #TTE_hat1 = (1/n)*(f_lin(1) - f_lin(0))
  TTE_hat2 = (1/n)*(f_spl(1) - f_spl(0))
  #return TTE_hat1, TTE_hat2
  return TTE_hat2


def poly_regression_prop(beta, y, A, z):
  '''
  Returns an estimate of the TTE using polynomial regression using
  numpy.linalg.lstsq

  beta (int): degree of polynomial
  y (numpy array): observed outcomes
  A (square numpy array): network adjacency matrix
  z (numpy array): treatment vector
  '''
  n = A.shape[0]

  if beta == 0:
      X = np.ones((n,2))
      X[:,1] = z
  else:
      X = np.ones((n,2*beta+1))
      count = 1
      treated_neighb = (A.dot(z)-z)/(np.array(A.sum(axis=1)).flatten()-1+1e-10)
      for i in range(beta):
          X[:,count] = np.multiply(z,np.power(treated_neighb,i))
          X[:,count+1] = np.power(treated_neighb,i+1)
          count += 2

  v = np.linalg.lstsq(X,y,rcond=None)[0]
  return np.sum(v)-v[0]

def poly_regression_prop_cy(beta, y, A, z):
  n = A.shape[0]
  X = np.ones((n,2*beta+2))
  z = z.reshape((n,1))
  treated_neighb = (A.dot(z)-z)/(np.array(A.sum(axis=1)).flatten()-1+1e-10)
  # temp = 1
  # for i in range(beta+1):
  #     X[:,i] = np.multiply(z,temp)
  #     X[:,beta+1+i] = np.multiply(1-z,temp)
  #     temp = temp * treated_neighb
  treated_neighb = np.power(treated_neighb.reshape((n,1)), np.arange(beta+1).reshape((1,beta+1)))
  X[:,:beta+1] = z.dot(treated_neighb)
  X[:,beta+1:] = (1-z).dot(treated_neighb)

  v = np.linalg.lstsq(X,y,rcond=None)[0]
  return np.sum(v[:beta+1])-v[beta+1]

def poly_regression_num(beta, y, A, z):
  '''
  Returns an estimate of the TTE using polynomial regression using
  numpy.linalg.lstsq

  beta (int): degree of polynomial
  y (numpy array): observed outcomes
  A (square numpy array): network adjacency matrix
  z (numpy array): treatment vector
  '''
  n = A.shape[0]

  if beta == 0:
      X = np.ones((n,2))
      X[:,1] = z
  else:
      X = np.ones((n,2*beta+1))
      count = 1
      treated_neighb = (A.dot(z)-z)
      for i in range(beta):
          X[:,count] = np.multiply(z,np.power(treated_neighb,i))
          X[:,count+1] = np.power(treated_neighb,i+1)
          count += 2

  # least squares regression
  v = np.linalg.lstsq(X,y,rcond=None)[0]

  # Estimate TTE
  count = 1
  treated_neighb = np.array(A.sum(axis=1)).flatten()-1
  for i in range(beta):
      X[:,count] = np.power(treated_neighb,i)
      X[:,count+1] = np.power(treated_neighb,i+1)
      count += 2
  TTE_hat = np.sum((X @ v) - v[0])/n
  return TTE_hat

def poly_regression_num_cy(beta, y, A, z):
  n = A.shape[0]

  X = np.ones((n,2*beta+2))
  z = z.reshape((n,1))
  treated_neighb = (A.dot(z)-z)
  # temp = 1
  # for i in range(beta+1):
  #     X[:,i] = np.multiply(z,temp)
  #     X[:,beta+1+i] = np.multiply(1-z,temp)
  #     temp = temp * treated_neighb
  treated_neighb = np.power(treated_neighb.reshape((n,1)), np.arange(beta+1).reshape((1,beta+1)))
  X[:,:beta+1] = z.dot(treated_neighb)
  X[:,beta+1:] = (1-z).dot(treated_neighb)

  # least squares regression
  v = np.linalg.lstsq(X,y,rcond=None)[0]

  # Estimate TTE
  X = np.zeros((n,2*beta+2))
  deg = np.array(A.sum(axis=1)).flatten()-1
  # temp = 1
  # for i in range(beta+1):
  #     X[:,i] = temp
  #     temp = temp * deg
  X[:,:beta+1] = np.power(deg.reshape((n,1)), np.arange(beta+1).reshape((1,beta+1)))
  TTE_hat = np.sum((X @ v) - v[beta+1])/n

  
  return TTE_hat


def Reg_beta(n, y, X, w):
  Xest = X*w.reshape(-1,1)
  Yest = y*w
  thetahat = np.linalg.inv(Xest.T @ Xest) @ Xest.T @ Yest
  est = np.sum(Yest) -  np.sum(Xest @ thetahat)
  return est/n


def get_c_est_deg2(A, z, y, treatment_vec, p, i):
    treatment_vec_cp = np.copy(treatment_vec)
    neighbor_idx = np.nonzero(A[[i],:].toarray()[0])[0].tolist()
    total_num = min(2, len(neighbor_idx))
    choice_list = []
    for choose_num in range(total_num + 1):
        choice_sublist = list(combinations(neighbor_idx, choose_num))
        choice_sublist = [list(tup) + [-1] * (total_num - choose_num) for tup in choice_sublist]
        choice_list.extend(choice_sublist)
    choice_list = np.array(choice_list, dtype=int)
    tilde_zi = np.array([np.product(z[
                    choice_list[x][choice_list[x] != -1]
                ]) for x in range(len(choice_list))])
    E_inv = np.zeros((len(choice_list), len(choice_list)))
    for s_id, s_vec in enumerate(choice_list):
        s = s_vec[s_vec != -1]
        for t_id, t_vec in enumerate(choice_list):
            t = t_vec[t_vec != -1]
            if len(s) == 0 and len(t) == 0:
                val = 1 + p / (1 - p) * len(neighbor_idx) + (len(neighbor_idx) - 1) * len(neighbor_idx) / 2 * (p / (1 - p)) ** 2
            elif len(s) == 1 and len(t) == 0:
                val = -1 / p * (p / (1 - p) + (p / (1 - p)) ** 2 * (len(neighbor_idx) - 1))
            elif len(s) == 0 and len(t) == 1:
                val = -1 / p * (p / (1 - p) + (p / (1 - p)) ** 2 * (len(neighbor_idx) - 1))
            elif len(s) == 2 and len(t) == 0:
                val = 1 / (1 - p) ** 2
            elif len(s) == 0 and len(t) == 2:
                val = 1 / (1 - p) ** 2
            elif len(s) == 1 and len(t) == 1:
                if s[0] == t[0]:
                    val = 1 / p ** 2 * (p / (1 - p) + (p / (1 - p)) ** 2 * (len(neighbor_idx) - 1))
                else:
                    val = 1 / (1 - p) ** 2
            elif len(s) == 1 and len(t) == 2:
                if s[0] in t:
                    val = -1 / (p * (1 - p) ** 2)
                else:
                    val = 0
            elif len(s) == 2 and len(t) == 1:
                if t[0] in s:
                    val = -1 / (p * (1 - p) ** 2)
                else:
                    val = 0
            elif len(s) == 2 and len(t) == 2:
                if (s == t).all():
                    val = 1 / (p ** 2 * (1 - p) ** 2)
                else:
                    val = 0

            E_inv[s_id, t_id] = val
          
          
    return (E_inv @ tilde_zi[:,None] * y[i]).reshape(-1)


def compute_component_B_deg2(X, A, p):
    n = len(X)
    ret = np.zeros((X.shape[1], X.shape[1]))
    for i in range(n):
        for ip in np.unique(A[:,A[[i],:].indices].nonzero()[0]):
            temp_len = len(set(A[[i],:].indices) & set(A[[ip],:].indices))
            ret += X[i][:, None] @ X[ip][None, :] * temp_len / (p * (1 - p))
            ret += X[i][:, None] @ X[ip][None, :] * temp_len * (temp_len - 1) / 2 * ((1 - 2*p) / (p * (1 - p)))**2
    return ret / n ** 2

def compute_component_D_deg2(X, A, p, c_est_list):
    '''
    c_est_list: a list of length n, the i-th element of the list is a vector (\hat{c}_{i,S})_S
    '''
    n = len(X)
    ret = np.zeros((X.shape[1],))
    for i in range(n):
        l = len(A[[i],:].indices)
        sp_to_sp_id = {neighbor : index for index, neighbor in enumerate(A[[i],:].indices)}
        pair_to_pair_id = {neighbor : [] for neighbor in A[[i],:].indices}
        for comb_id, comb in enumerate(list(combinations(A[[i],:].indices,2))):
            pair_to_pair_id[comb[0]].append(comb_id + 1 + l)
            pair_to_pair_id[comb[1]].append(comb_id + 1 + l)
        
        beta1_id = len(set(A[[i],:].indices))
        c_est_null = c_est_list[i][0] + p * sum(c_est_list[i][1:(beta1_id+1)]) + sum(c_est_list[i][(beta1_id +1):]) * p ** 2

        for ip in np.unique(A[:,A[[i],:].indices].nonzero()[0]):
            commons = set(A[[i],:].indices) & set(A[[ip],:].indices)
            pair_to_pair_id_with_common_neighbors = {neighbor : [] for neighbor in A[[i],:].indices}
            pair_to_pair_id_no_common_neighbors = {neighbor : [] for neighbor in A[[i],:].indices}
            for comb_id, comb in enumerate(list(combinations(A[[i],:].indices,2))):
                if comb[0] in commons and comb[1] in commons:
                    pair_to_pair_id_with_common_neighbors[comb[0]].append(comb_id + 1 + l)
                    pair_to_pair_id_with_common_neighbors[comb[1]].append(comb_id + 1 + l)
                else:
                    pair_to_pair_id_no_common_neighbors[comb[0]].append(comb_id + 1 + l)
                    pair_to_pair_id_no_common_neighbors[comb[1]].append(comb_id + 1 + l)

            temp_len = len(commons)
            ret += X[ip] * c_est_null * temp_len / (p * (1 - p))
            
            ret += X[ip] * c_est_null * temp_len * (temp_len - 1) / 2 * ((1 - 2*p) / (p * (1 - p)))**2
            
            for sp in set(A[[i],:].indices):
                sp_id = sp_to_sp_id[sp]
                c_est_first_order = c_est_list[i][sp_id + 1] + p*np.sum(c_est_list[i][pair_to_pair_id[sp]])
                if sp in commons:
                    ret += X[ip] * c_est_first_order * (1-2*p) / (p * (1 - p))
                    ret += X[ip] * c_est_first_order * (1-2*p) / (p * (1 - p)) * 2 * (temp_len-1)
                    ret += X[ip] * c_est_first_order * (1-2*p)**3 / (p * (1 - p))**2 * (temp_len-1)
                    ret += X[ip] * (np.sum(c_est_list[i][pair_to_pair_id[sp]]))
                    ret += X[ip] * (np.sum(c_est_list[i][pair_to_pair_id_no_common_neighbors[sp]])) * (1-2*p)**2 / (p * (1 - p))
                    ret += X[ip] * (np.sum(c_est_list[i][pair_to_pair_id_with_common_neighbors[sp]])) * (1-2*p)**2 / (p * (1 - p)) * 2
                    ret += X[ip] * (np.sum(c_est_list[i][pair_to_pair_id_with_common_neighbors[sp]])) * (1-2*p)**4 / (p * (1 - p))**2 / 2
                    ret += X[ip] * (np.sum(c_est_list[i][pair_to_pair_id_no_common_neighbors[sp]])) * (1-2*p)**2 / (p * (1 - p)) * (temp_len-1)
                    ret += X[ip] * (np.sum(c_est_list[i][pair_to_pair_id_with_common_neighbors[sp]])) * (1-2*p)**2 / (p * (1 - p)) * (temp_len-2)
                else:
                    ret += X[ip] * c_est_first_order * (1-2*p) / (p * (1 - p)) * temp_len
                    
            
    return 2 * ret/ n ** 2

def compute_component_D_deg2_old(X, A, p, c_est_list):
    '''
    c_est_list: a list of length n, the i-th element of the list is a vector (\hat{c}_{i,S})_S
    '''
    n = len(X)
    ret = np.zeros((X.shape[1],))
    for i in range(n):
        l = len(A[[i],:].indices)
        sp_to_sp_id = {neighbor : index for index, neighbor in enumerate(A[[i],:].indices)}
        pair_to_pair_id = {neighbor : {neighbor_p : int((2*l-1-index) * index / 2 + index_p - index - 1) for index_p, neighbor_p in enumerate(A[[i],:].indices) if index_p > index} for index, neighbor in enumerate(A[[i],:].indices)}
        for ip in np.unique(A[:,A[[i],:].indices].nonzero()[0]):
            beta1_id = len(set(A[[i],:].indices))
            temp_len = len(set(A[[i],:].indices) & set(A[[ip],:].indices))
            c_est_null = c_est_list[i][0] + p * sum(c_est_list[i][1:(beta1_id+1)]) + sum(c_est_list[i][(beta1_id +1):]) * p ** 2
            ret += X[ip] * c_est_null * temp_len / (p * (1 - p))
            ret += X[ip] * c_est_null * temp_len * (temp_len - 1) / 2 * (1 / p ** 3 + 1 / (1 - p) ** 3) * (1 - 2*p) ** 2
            for sp in set(A[[i],:].indices) & set(A[[ip],:].indices):
                sp_id = sp_to_sp_id[sp]
                ret += X[ip] * (c_est_list[i][sp_id + 1] + p*(np.sum(c_est_list[i][list(pair_to_pair_id[sp].values())]))) * (1-2*p) / (p * (1 - p))
            ret += X[ip] * (np.sum(c_est_list[i][l + 1:])) * ((1-p)**2 / p ** 3 + p**2 / (1 - p) ** 3) * (1 - 2*p) ** 2
    #t2 = time.time()
    #print("component D", t2-t1)
    return 2 * ret/ n ** 2 


def simpleWeights_deg2(X, A, diag=5, offdiag=5, rand_diag=np.array([]), rand_offdiag=np.array([])):
    '''
    Returns weights generated from simpler model

    A (numpy array): adjacency matrix of the network
    diag (float): maximum norm of direct effects
    offidiag (float): maximum norm of the indirect effects
    '''
    n = A.shape[0]

    if rand_offdiag.size == 0:
        rand_offdiag = np.random.rand(n)
    C_offdiag = offdiag*rand_offdiag

    in_deg = scipy.sparse.diags(np.array(A.sum(axis=1)).flatten(),0)  # array of the in-degree of each node
    C = in_deg.dot(A - scipy.sparse.eye(n))
    col_sum = np.array(C.sum(axis=0)).flatten()
    col_sum[col_sum==0] = 1
    temp = scipy.sparse.diags(C_offdiag/col_sum)
    C = C.dot(temp)

    if rand_diag.size == 0:
        rand_diag = np.random.rand(n)
    C_diag = (np.sum(X, axis=1)+diag)*rand_diag
    C.setdiag(C_diag)

    return C


