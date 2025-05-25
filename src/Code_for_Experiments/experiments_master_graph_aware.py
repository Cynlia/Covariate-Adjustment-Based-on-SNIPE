'''
Script to run (three) experiments for different values of beta

Running this script will save the (new) results into outputFiles/new/
To view the original results from the paper, see outputFiles/graph_aware/
'''

# Setup
import numpy as np
import pandas as pd
import sys
import time
import scipy.sparse
import nci_linear_setup as ncls
import nci_polynomial_setup as ncps
import os
import pickle
from multiprocessing import Pool
from functools import partial
#num_cores = int(os.getenv('SLURM_CPUS_PER_TASK'))
import warnings
warnings.filterwarnings("ignore")



path_to_module = 'Code-for-Experiments/'
#save_path = 'outputFiles/new/'
save_path = 'outputFiles/graph_aware/'
save_path_graphs = 'graphs/'

# Read seed from command-line argument
# seed = int(sys.argv[1])

#############################################
# parameters for covariate adjustment
#############################################

covariate_adj = True
covariate_dim = 3


def main(argv):
    if len(argv) > 1:
        beta = int(argv[0])
    else:
        beta = 2

    G = 1         # number of graphs we want to average over (10)
    T = 1          # number of trials per graph (500)

    graphStr = "srgg"

    for beta in [2]:

        f = open(save_path+'experiments_output_deg'+str(beta)+'_SNIPE'+'.txt', 'w')
        startTime1 = time.time()

        ###########################################
        # Run Experiment 1: Varying Size of Network
        ###########################################
        if True:
            diag = 10       # controls magnitude of direct effects
            r = 2           # ratio between indirect and direct effects
            p = 0.2         # treatment probability
            pct = 1

            results = []
            sizes = np.array([5000, 6000, 7000, 8000, 9000, 10000]) #    [500, 600, 700, 800, 900]
            if beta == 1:
                sigmas = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02])
            elif beta == 2:
                sigmas = np.array([0.02, 0.018, 0.016, 0.016, 0.014, 0.014])

            for n, sigma in zip(sizes, sigmas):
                print("n = {}".format(n))
                startTime2 = time.time()

                results.extend(run_experiment(G,T,n,p,r,sigma,pct,graphStr,diag,beta))

                executionTime = (time.time() - startTime2)
                print('Runtime (in seconds) for n = {} step: {}'.format(n,executionTime),file=f)
                print('Runtime (in seconds) for n = {} step: {}'.format(n,executionTime))

            executionTime = (time.time() - startTime1)
            print('Runtime (size experiment) in minutes: {}'.format(executionTime/60),file=f)  
            print('Runtime (size experiment) in minutes: {}\n'.format(executionTime/60))       
            df = pd.DataFrame.from_records(results)
            df.to_csv(save_path+graphStr+'-size-deg'+str(beta)+'-SNIPE'+'.csv')

        ################################################
        # Run Experiment: Varying Treatment Probability 
        ################################################
        if False:    
            startTime2 = time.time()
            if graphStr == "sw":
                n = 96
            else:
                n = 1   # number of nodes in network, default 500
            diag = 10       # maximum norm of direct effect
            r = 2           # ratio between indirect and direct effects
            sigma = 0.02
            pct = 1

        
            results = []
            p_treatments = np.array([0.1, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]) # treatment probabilities [0.1, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]

        
            for p in p_treatments:
                print("Treatment Probability: {}".format(p))
                startTime3 = time.time()

                results.extend(run_experiment(G,T,n,p,r,sigma,pct,graphStr,diag,beta))

                executionTime = (time.time() - startTime3)
                print('Runtime (in seconds) for p = {} step: {}'.format(p,executionTime),file=f)
                print('Runtime (in seconds) for p = {} step: {}'.format(p,executionTime))

            executionTime = (time.time() - startTime2)
            print('Runtime (tp experiment) in minutes: {}'.format(executionTime/60),file=f)  
            print('Runtime (tp experiment) in minutes: {}\n'.format(executionTime/60))           
            df = pd.DataFrame.from_records(results)
            df.to_csv(save_path+graphStr+'-tp-deg'+str(beta)+'-SNIPE-'+str(seed)+'.csv')

        ###########################################################
        # Run Experiment: Varying Ratio of Indirect & Direct Effects 
        ###########################################################
        if False:
            startTime2 = time.time()
            if graphStr == "sw":
                n = 96
            else:
                n = 1     # number of nodes in network, default 500
            p = 0.2             # treatment probability
            diag = 10           # maximum norm of direct effect
            sigma = 0.02
            pct = 1

            results = []
            ratio = [0.01, 0.1, 0.25,0.5,0.75,1,1/0.75,1/0.5,3,1/0.25] #[0.01, 0.1, 0.25,0.5,0.75,1,1/0.75,1/0.5,3,1/0.25]

        
            for r in ratio:
                print('ratio: {}'.format(r))
                startTime3 = time.time()

                results.extend(run_experiment(G,T,n,p,r,sigma,pct,graphStr,diag,beta))

                executionTime = (time.time() - startTime3)
                print('Runtime (in seconds) for r = {} step: {}'.format(r,executionTime),file=f)
                print('Runtime (in seconds) for r = {} step: {}'.format(r,executionTime))

            executionTime = (time.time() - startTime2)
            print('Runtime (ratio experiment) in minutes: {}'.format(executionTime/60),file=f)   
            print('Runtime (ratio experiment) in minutes: {}\n'.format(executionTime/60))           
            df = pd.DataFrame.from_records(results)
            df.to_csv(save_path+graphStr+'-ratio-deg'+str(beta)+'-SNIPE-'+str(seed)+'.csv')

        ###########################################################
        # Run Experiment: Varying Percent of Covariate Explanation  
        ###########################################################
        if False:
            startTime2 = time.time()
            if graphStr == "sw":
                n = 96
            else:
                n = 5000     # number of nodes in network, default 500
            p = 0.2             # treatment probability
            diag = 10           # maximum norm of direct effect
            sigma = 0.02
            r = 2           # ratio between indirect and direct effects

            results = []
            pct_X = [0, 0.25,0.5,0.75,1,1/0.75,1/0.5,1/0.25,10]

            for pct in pct_X:
                print('percent: {}'.format(pct))
                startTime3 = time.time()

                results.extend(run_experiment(G,T,n,p,r,sigma,pct,graphStr,diag,beta))

                executionTime = (time.time() - startTime3)
                print('Runtime (in seconds) for pct = {} step: {}'.format(r,executionTime),file=f)
                print('Runtime (in seconds) for pct = {} step: {}'.format(r,executionTime))

            executionTime = (time.time() - startTime2)
            print('Runtime (percent experiment) in minutes: {}'.format(executionTime/60),file=f)   
            print('Runtime (percent experiment) in minutes: {}\n'.format(executionTime/60))           
            df = pd.DataFrame.from_records(results)
            df.to_csv(save_path+graphStr+'-percent-deg'+str(beta)+'-SNIPE-'+str(seed)+'.csv')

        executionTime = (time.time() - startTime1)
        print('Runtime (whole script) in minutes: {}'.format(executionTime/60),file=f)
        print('Runtime (whole script) in minutes: {}'.format(executionTime/60))

        f.close()

def run_experiment(G,T,n,p,r,sigma,pct,graphStr,diag=1,beta=2):
    
    offdiag = r*diag   # maximum norm of indirect effect

    results = []
    dict_base = {'p': p, 'ratio': r, 'n': n, 'beta': beta, 'pct': pct}

    sz = str(n) + '-'
    
    np.random.seed(2025)
    #############################################
    # generation setting for covariate adjustment
    #############################################
    b = 5*pct*np.ones((covariate_dim,))
    X_uncentered = np.random.randn(n, covariate_dim) 
    X_c_coeff = np.ones((covariate_dim, n))+np.random.randn(covariate_dim, n)
    X = X_uncentered - np.mean(X_uncentered, axis=0)

    for g in range(G):
        graph_rep = str(g)
        dict_base.update({'Graph':sz+graph_rep})

        if graphStr == "er":
            deg = 10 # default 10
            A = ncls.erdos_renyi(n,deg/n)
        elif graphStr == "srgg":
            A = ncls.soft_RGG(X,n,sigma)

        rand_wts = np.random.rand(n,3)
        
    
        alpha = rand_wts[:,0].flatten()
        
        if beta == 1: 
            C = ncls.simpleXWeights_old(A, X, X_c_coeff, diag, offdiag, rand_wts[:,1].flatten(), rand_wts[:,2].flatten())
        if beta == 2:
            C_linear = ncls.simpleXWeights_old(A, X, X_c_coeff, diag, offdiag, rand_wts[:,1].flatten(), rand_wts[:,2].flatten())
            C = ncps.simpleWeights_deg2(X, A, diag, offdiag, rand_wts[:,1].flatten(), rand_wts[:,2].flatten())

        if beta == 1:
            fy = lambda z: ncls.linear_cov_adj(b,C,alpha,z,X)
        else:
            fy = ncps.ppom_cov_adj(beta, C, alpha, X, b, C_linear)

        # compute and print true TTE
        TTE = 1/n * np.sum((fy(np.ones(n)) - fy(np.zeros(n))))
        dict_base.update({'TTE':TTE})


        ####### Estimate ########
        estimators = []
        if beta == 1:
            estimators.append(lambda y,z,w: ncps.Reg_beta(n, y, X, w))
            estimators.append(lambda y,z,w: ncls.VIM_beta(n, y, X, w, components))
            estimators.append(lambda y,z,w: ncls.SNIPE_deg1(n,y,w))
            estimators.append(lambda y,z,w: ncls.SNIPE_beta_Lin(n, y, X, w, z, p))
        else:
            estimators.append(lambda y,z,w: ncps.Reg_beta(n, y, X, w))
            estimators.append(lambda y,z,w: ncls.VIM_beta(n, y, X, w, components))
            estimators.append(lambda y,z,w: ncps.SNIPE_beta(n,y,w))
            estimators.append(lambda y,z,w: ncls.SNIPE_beta_Lin(n, y, X, w, z, p))
            
        estimators.append(lambda y,z,w: ncls.diff_in_means_naive(y,z))
        alg_names = ['Reg', 'VIM', 'SNIPE('+str(beta)+')', 'Lin\'s', 'DM']

        
        if beta == 1:
            component_B = ncls.compute_component_B_deg1(X, A, p)
        else:
            component_B = ncps.compute_component_B_deg2(X, A, p)


        for i in range(T):
            #np.random.seed(seed)
            dict_base.update({'rep':i, 'Rand': 'Bernoulli'})
            z = ncls.bernoulli(n,p)
            y = fy(z)
            if beta == 1:
                treatment_vec = p * np.ones((n,))
                c_est_list = [ncls.get_c_est(A, z, y, treatment_vec, i) for i in range(n)]
                component_D = ncls.compute_component_D_deg1(X, A, p, c_est_list)
                zz = z/p - (1-z)/(1-p)
                w = A.dot(zz)
            else:
                treatment_vec = p * np.ones((n,))
                w = ncps.SNIPE_weights(n, p, A, z, beta)
                c_est_list = [ncps.get_c_est_deg2(A, z, y, treatment_vec, p, i) for i in range(n)]
                component_D = ncps.compute_component_D_deg2(X, A, p, c_est_list)
                
            components = (component_B, component_D)            

            for ind in range(len(estimators)):
                est = estimators[ind](y,z,w)
                dict_base.update({'Estimator': alg_names[ind], 
                                  'TTE_Estimate': est,
                                  'Absolute_Bias': est-TTE,
                                  'Bias_squared': (est-TTE)**2,
                                  'Relative_Bias': (est-TTE)/TTE})
                results.append(dict_base.copy())

    return results


if __name__ == "__main__":
    main(sys.argv[1:])
