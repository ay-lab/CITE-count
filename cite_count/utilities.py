import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import stats
from Bio import SeqIO
from Bio import pairwise2
from scipy.stats import norm,rayleigh,gamma,chi2
from scipy.optimize import least_squares
from sklearn.mixture import BayesianGaussianMixture
from math import ceil
from tqdm import tqdm

# From http://cf.10xgenomics.com/supp/cell-exp/notebook_tutorial-3.0.0.html
import collections
import scipy.sparse as sp_sparse
import h5py

# Functions in AB_count
def Threshold_finder(data,max_population = 3, min_population_size = 0.2, confidence_interval = 0.90, verbose = False):
    import warnings
    warnings.filterwarnings("ignore")
    '''
    data: 1D data array with count numbers
    max_population: Define the maximal number of populations exist in sample datasets
    min_population_size: The smallest population should have at least 20% population
    confidence_interval: if unimodel was used, select the confidence interval for lower bound; 0.90 = 5% confidence one tail test
    '''

    best_population = np.inf
    best_loglike = -np.inf
    best_mdoel = None
    model_kind = 'Gaussian' # Set Gaussian to be the default model type

    for n_components in [n+1 for n in list(reversed(np.arange(max_population)))]:
        BGM = BayesianGaussianMixture(n_components=n_components,verbose=0).fit(data)
        # Proceed only if the model can converge
        if BGM.converged_:
            if verbose:
                print('%s populations converged'%str(n_components))
            dict_wp = dict() # store weighted probability for each population
            for p in np.arange(n_components):
                para = norm.fit(mask_list(data,BGM.predict(data),p)) # fit gaussian model to population p
                dict_wp[p] = norm(para[0],para[1]).pdf(data)*BGM.weights_[p]
            # Compute log likelyhood of prediction
            # wp[0] = norm.pdf(data[i])*weight[0], wp[1] = norm.pdf(data[i])*weight[1] ...
            # log(wp[0]+wp[1]+...) gives total log likelyhood
            loglike = sum([np.log(sum([dict_wp[p][i] for p in np.arange(n_components)])) for i in np.arange(len(data))])[0]
            if loglike > best_loglike and min(BGM.weights_) > min_population_size: # minimal 
                best_loglike = loglike
                best_population = n_components
                best_mdoel = BGM
            if verbose:
                print('%s model with %s population has log likelyhood of %s '%(model_kind,n_components,loglike))    
        else:
            if verbose:
                print('%s populations not converged'%str(n_components))

        if n_components == 1: # A gaussian model may not best fit one distribution; Other models should also being tested to decide if better 1 model fit exist 

            para = rayleigh.fit(data)
            loglike = sum(np.log(rayleigh(para[0],para[1]).pdf(data)))[0]
            if loglike > best_loglike:
                best_loglike = loglike
                best_population = 1
                best_mdoel = rayleigh(para[0],para[1])
                model_kind = 'Rayleigh'
                if verbose:
                    print('%s model with %s population has log likelyhood of %s '%(model_kind,n_components,loglike))
                
    if best_mdoel == None: # nither Gaussian nor Rayleight could fit the data
        para = chi2.fit(data)
        loglike = sum(np.log(chi2(para[0],para[1],para[2]).pdf(data)))[0]
        if loglike > best_loglike:
            best_loglike = loglike
            best_population = 1
            best_mdoel = chi2(para[0],para[1],para[2])
            model_kind = 'Chi-square'
            if verbose:
                print('%s model with %s population has log likelyhood of %s '%(model_kind,n_components,loglike))



    if best_population > 1:
        p = list(best_mdoel.means_).index(min(best_mdoel.means_)) # Get the population id that represent negatives
        threshold = max(mask_list(data,best_mdoel.predict(data),p))[0]
    else:
        if model_kind == 'Rayleigh' or model_kind == 'Chi-square':    
            threshold = min(1,abs(best_mdoel.interval(confidence_interval)[0]))
        else:
            para = norm.fit(data)
            threshold = min(1,abs(norm(data,para[0],para[1]).interval(confidence_interval)[0]))


    print('Best model with %s distribution has %s populations with threshold at %s'%(model_kind,best_population,threshold))

    return threshold,model_kind,best_mdoel,best_population



def Hawk_smash(List):
    '''Flattern lists within a list '''
    return [item for sublist in List for item in sublist]

def mask_list(data, selectors, key_value = -np.inf):
    # Mask a data list based on the boolen selector
    if key_value == -np.inf:
        return [d for d, s in zip(data, selectors) if s]
    else:
        return [d for d, s in zip(data, selectors) if s == key_value]

# Transform vector to matrix and fill 0 if barcode-tag combination not exist
def match_code(cell,tag,df_merge_UMI):
    try:
        return df_merge_UMI.loc[cell+tag]['count']
    except:
        return 0
    
def flip_boolen(val):
    return 1 - val

def platt(coef,X,Y):
    # Modified Platt equation to allow more shape turn at the tipping point
    return coef[0]*(1-10**(-(coef[1]*X)/coef[0]))- Y

def fit_platt(diff):
    # guess init aplha value and set point 1 to zero
    G_alpha = stats.linregress(diff[0][0:5],diff[1][0:5]).slope
    Offset = -diff[1][0]
    G_P = Offset
    coef_0 = np.array([G_P,G_alpha], dtype=float)
    # Start non-linear fitting
    res_lsq = least_squares(platt, coef_0, loss='linear',  args=(diff[0], diff[1]+Offset))

    return res_lsq.x

def get_tag(seq_f,seq_r,tag_list,AB_list,project_name,overwrite = False):
        
    # Check if file has already been processed
    if os.path.isfile('cache/%s/UMI_raw.pkl'%project_name) and overwrite == False:
        df_merge = pd.read_pickle('cache/%s/UMI_raw.pkl'%project_name)

    else:

        # Read in forward and reverse sequences
        print('Read forward sequene')
        df_f = pd.DataFrame([[record.description.split(' ')[0],str(record.seq[:16]),str(record.seq[16:])] for record in tqdm(SeqIO.parse(seq_f, "fastq"))],columns=['ID','cell_code','UMI']).set_index('ID')
        print('Read reverse sequene')
        df_r = pd.DataFrame([[record.description.split(' ')[0],str(record.seq[:15]),str(record.seq[:40])] for record in tqdm(SeqIO.parse(seq_r, "fastq")) if str(record.seq[20:40]).count('A') > 16],columns=['ID','tag_code','ext']).set_index('ID')

        # Remove duplicated UMI
        df_merge_pre_recovery = df_f.merge(df_r,left_index=True, right_index=True).groupby(['cell_code','UMI','tag_code']).first().reset_index()

        # Get length info for stats
        ori_tag = len(df_merge_pre_recovery.set_index('tag_code').loc[AB_list])
        ori_UMI = len(df_merge_pre_recovery)

        # See how the overall threshold for tag recovery was determined - Notebook - Antibody_tag_recovery_benchmark 
        test_list = list(set(df_merge_pre_recovery.tag_code))
        score = [[pairwise2.align.globalms(query,ref,2, -1, -1, -1)[0][2] for ref in AB_list] for query in test_list]
        dict_fix_seq = dict(zip(test_list,[AB_list[int(line.index(max(line)))] if max(line) > 25 else np.nan for line in score]))

        df_merge = df_merge_pre_recovery.copy()
        df_merge['tag_code'] = [dict_fix_seq[x] for x in df_merge_pre_recovery['tag_code']]
        df_merge = df_merge.loc[df_merge['tag_code'].dropna().index]
        df_merge.to_pickle('cache/%s/UMI_raw.pkl'%project_name)
        
        print('Total reads per library = ' "{:,}".format(len(df_f)))
        print('Sequences removed by poly A filter = ' "{0:.3%}".format(1-len(df_r)/len(df_f)))
        print('Reads removed by UMI filter = ' "{0:.3%}".format(len(df_r)/len(df_f) - len(df_merge)/len(df_f)))
        print('Total filtered UMI reads per library = ' "{:,}".format(len(df_merge)))
        print('Reads with no ref tag in UMI counts = '+"{0:.3%}".format(1-len(df_merge_pre_recovery.set_index('tag_code').loc[AB_list])/len(df_merge_pre_recovery)))
        print('Reads recovered in UMI counts = '+"{0:.3%}".format((len(df_merge)-ori_tag)/ori_UMI))
    return df_merge

def find_min_count(count_data):
        return fit_platt(count_data)[0]/fit_platt(count_data)[1]
    
# Functions in CITE_count  

def dist_exclusion(A,B):
    return 4*len(set(A)-set(B))*len(set(B)-set(A))/len(set(A).union(set(B)))**2

def get_pair(df_dist,rank):
    list_tag = list(np.sort([x for x in Hawk_smash(df_dist.values) if x > 0]))
    list_tag.reverse()
    target_value = list_tag[rank]
    df_dist = df_dist[df_dist == target_value].dropna(how = 'all').dropna(axis = 1)
    return df_dist.index[0],df_dist.columns[0],target_value
    
def group_cell(df_CITE):
    dict_label = dict(zip(df_CITE.columns,[list(df_CITE[df_CITE[col] == 1].index) for col in df_CITE.columns]))
    df_dist = pd.DataFrame([[dist_exclusion(dict_label[A],dict_label[B]) for A in dict_label] for B in dict_label],columns=dict_label.keys(),index=dict_label.keys())
    return df_dist.where(np.triu(np.ones(df_dist.shape)).astype(np.bool))

def find_node(name,root_tree):
    return [obj for obj in Hawk_smash(root_tree) if obj.name == name][0]

def get_branch(obj,root_tree):
    if obj.Node.depth == 0:
        return []
    elif obj.Node.depth > 0:
        pair = [obj.pair]
        for depth in np.arange(obj.Node.depth-1):
            parent = find_node(obj.Node.parent.name.split('(')[0][:-1],root_tree)
            obj = parent
            if len([obj.pair])>0:
                pair.append(obj.pair)
    else: 
        pair = [obj.pair]
    return pair

def get_key_by_name(dict_item,target):
    return [A for A, B in dict_item.items() if B == target]

def get_geneid_by_name(list_gene,dict_item):
    return Hawk_smash([get_key_by_name(dict_item,name) for name in list_gene])

# From http://cf.10xgenomics.com/supp/cell-exp/notebook_tutorial-3.0.0.html
import collections
import scipy.sparse as sp_sparse
import h5py


FeatureBCMatrix = collections.namedtuple('FeatureBCMatrix', ['feature_ids', 'feature_names', 'barcodes', 'matrix'])

def get_matrix_from_h5(filename):
    with h5py.File(filename,'r') as f:
        if u'version' in f.attrs:
            if f.attrs['version'] > 2:
                raise ValueError('Matrix HDF5 file format version (%d) is an newer version that is not supported by this function.' % version)
        else:
            raise ValueError('Matrix HDF5 file format version (%d) is an older version that is not supported by this function.' % version)
        
        feature_ids = [x.decode('ascii', 'ignore') for x in f['matrix']['features']['id']]
        feature_names = [x.decode('ascii', 'ignore') for x in f['matrix']['features']['name']]        
        barcodes = list(f['matrix']['barcodes'][:])
        matrix = sp_sparse.csc_matrix((f['matrix']['data'], f['matrix']['indices'], f['matrix']['indptr']), shape=f['matrix']['shape'])
        return FeatureBCMatrix(feature_ids, feature_names, barcodes, matrix)

def get_expression(fbm, gene_name):
    try:
        gene_index = feature_bc_matrix.feature_names.index(gene_name)
    except ValueError:
        raise Exception("%s was not found in list of gene names." % gene_name)
    return fbm.matrix[gene_index, :].toarray().squeeze()