from utilities import *

import scanpy as sc
import numpy as np
import pandas as pd
from anytree import Node,RenderTree
from anytree.exporter import DotExporter
from graphviz import render
from tqdm import tqdm
from sklearn.decomposition import PCA,TruncatedSVD
import umap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib_venn import venn2
import seaborn as sns
import pickle

import warnings


class AB_Node:
    def __init__(self,name,parent = None,list_cell = [],pair = [],group_id = None):
        self.name = name
        if parent == None:
            self.Node = Node(name +' ('+str(len(list_cell))+')')
        else:
            self.Node = Node(name +' ('+str(len(list_cell))+')',parent = parent)
        self.list_cell = list_cell
        self.pair = pair
        self.group_id = group_id

class RNA_count:
    
    'Read in dataframe from Seurat output or count table'
    
    def __init__(self,PROJ_NAME):
        self.PROJ_NAME = PROJ_NAME
        os.system('mkdir -p ./cache/%s'%PROJ_NAME)
        self.data_type = 'RNA'
    
    def import_Seurat(self,df_count):
        'count table for checking RNA data and scale data for PCA and clustering'
        self.df_RNA_count = df_count.T
        df_count = None
        self.list_barcodes = list(self.df_RNA_count.index)
        
    def import_count_table(self,count_table,overwrite = False):
        '''
        Import count table with cells in index and genes in columns (flip if needed)
        '''
        symbol = lambda x: {'csv':',','tsv':'\t'}[x]
        self.df_RNA_count = pd.read_csv(count_table,sep = symbol(count_table.split('.')[-1]),index_col = 0)
        self.df_RNA_count = self.df_RNA_count[self.df_RNA_count.sum(axis = 1)>0]
        self.list_barcodes = list(self.df_RNA_count.index)
        
    def import_h5(self,path_h5,overwrite = False):
        self.feature_bc_matrix = get_matrix_from_h5(path_h5)
        self.list_barcodes = [str(x).split('-')[0][2:] for x in self.feature_bc_matrix.barcodes]
        
        if os.path.isfile('cache/%s/df_RNA_count.pkl'%self.PROJ_NAME) and overwrite == False:
            self.df_RNA_count = pd.read_pickle('cache/%s/df_RNA_count.pkl'%self.PROJ_NAME)
            
        else:
            self.df_RNA_count = pd.DataFrame([list(self.feature_bc_matrix.matrix[:,cell].T.toarray()[0]) for cell in tqdm(np.arange(len(self.list_barcodes)))],index=self.list_barcodes,columns=self.feature_bc_matrix.feature_names)
            self.df_RNA_count = self.df_RNA_count.T[self.df_RNA_count.T.sum(axis = 1)>0].T
            self.df_RNA_count.to_pickle('cache/%s/df_RNA_count.pkl'%self.PROJ_NAME)
    
    def save_cache(self):
        with open('cache/%s/cache_RNA.pkl'%self.PROJ_NAME, 'wb') as f:
            pickle.dump([self.data_type,self.list_barcodes,self.df_RNA_count],f)

    def load_cache(self):
        try:
            with open('cache/%s/cache_RNA.pkl'%self.PROJ_NAME, 'rb') as f:
                self.data_type,self.list_barcodes,self.df_RNA_count = pickle.load(f)
        except OSError:
            print('Loading error, cache file may not exist.')
        
    def qc_results(self,bins = 100,dpi = 100):
        umis_per_cell = self.df_RNA_count
        f,ax = plt.subplots(1,2,figsize = (10,3),dpi = 80)
        ax[0].hist(np.log10(umis_per_cell.sum(axis = 1)), bins=bins)
        ax[0].set_xlabel('UMI counts per cell (log10)')
        ax[0].set_ylabel('Frequency')
        ax[0].set_title('UMI Distribution')
        ax[1].hist(np.log10((umis_per_cell>0).sum(axis = 1)), bins=bins)
        ax[1].set_xlabel('gene counts per cell (log10)')
        ax[1].set_title('gene counts Distribution')
        ax[1].set_ylabel('Frequency')
        
class AB_count:

    "Process raw sequence files and predict labels for each barcoded cell"

    def __init__(self,PROJ_NAME):
        "Check if file has already been processed and init run"
        # direct input
        self.PROJ_NAME = PROJ_NAME
        self.data_type = 'AB'
        self.high_signal_list = []
        os.system('mkdir -p ./cache/%s'%PROJ_NAME)
        
        self.dict_threshold = dict() # for +/- singals
        
    def process_fastq(self,seq_f,seq_r,tag_list,AB_list,overwrite = False):    
        self.AB_list = AB_list
        
        # get tag annotation
        df_tag = pd.read_csv(tag_list).set_index('Barcode Sequence')
        df_tag['Barcode'] = df_tag['Barcode'].astype(int).apply(lambda x: str(f"{x:04d}"))
        df_tag['annotation'] = df_tag['Specificity']
        self.dict_tag = dict(zip(list(df_tag.index),df_tag['annotation']))
        
        # data import and pre-process 
        self.df_merge = get_tag(seq_f,seq_r,tag_list,self.AB_list,self.PROJ_NAME,overwrite = overwrite)
        self.df_tag_count = self.df_merge.groupby('cell_code').size().reset_index(name = 'count').set_index('cell_code')
        self.log_tag_data = [np.arange(0,100),np.log10([len(self.df_tag_count[self.df_tag_count['count']>i]) for i in np.arange(0,100)])]
        
    def import_counts(self,count_table):
        symbol = lambda x: {'csv':',','tsv':'\t'}[x]
        self.df_merge_UMI_mx = pd.read_csv(count_table,sep = symbol(count_table.split('.')[-1]),index_col = 0)
        self.AB_list = list(self.df_merge_UMI_mx.columns)
        self.dict_tag = dict(zip(self.AB_list,self.AB_list))
        self.df_tag_count = pd.DataFrame(self.df_merge_UMI_mx.sum(axis = 1),columns=['count'])
        self.log_tag_data = [np.arange(0,max(self.df_tag_count.values)[0],round(max(self.df_tag_count.values)[0]/100)),np.log10([len(self.df_tag_count[self.df_tag_count['count']>i]) for i in np.arange(0,max(self.df_tag_count.values)[0],round(max(self.df_tag_count.values)[0]/100))])]
    
    def qc_tag_count(self,figsize = (5,3),dpi = 100, THRESHOLD_UMI = None,saturation = 'half'):
        "Plot minimal count threshold for QC; choose half or full saturation setting"
        
        self.THRESHOLD_UMI = THRESHOLD_UMI
        if saturation == 'half':
            S_factor = 1
        elif saturation == 'full':
            S_factor = 2
        plt.subplots(1,1,figsize = figsize,dpi = dpi)
        plt.plot(self.log_tag_data[0],self.log_tag_data[1],'.')
        if self.THRESHOLD_UMI == None:
            self.THRESHOLD_UMI = S_factor*find_min_count(self.log_tag_data)
            plt.axvline(x = self.THRESHOLD_UMI,ymax=1,linestyle = '--')
            plt.plot(self.log_tag_data[0],platt(fit_platt(self.log_tag_data),self.log_tag_data[0],-self.log_tag_data[1][0]))
            plt.text(self.THRESHOLD_UMI+1,np.mean(self.log_tag_data[1])*1.2,'<- Auto threshold = %.1f'%self.THRESHOLD_UMI)
        else:
            plt.axvline(x = self.THRESHOLD_UMI,ymax=1,linestyle = '--')
            plt.text(self.THRESHOLD_UMI+1,np.mean(self.log_tag_data[1])*1.2,'<- Input threshold = %.1f'%self.THRESHOLD_UMI)
                
        # plt.xticks(np.arange(0,max(self.df_tag_count.values)[0],round(max(self.df_tag_count.values)[0]/1000)))
        plt.xlabel('minimal tag UMI counts per barcode')
        plt.ylabel('Total barcodes (log10)');
    
    def filter_tag_count(self, P_OUTLIER = 1e-6, figsize = (10,3),dpi = 100,view = True):
        # Pre-filter low count sequence based on low count QC result (auto filter if QC not done)
        
        if self.THRESHOLD_UMI == None:
            self.THRESHOLD_UMI = find_min_count(self.log_tag_data)
        
        self.df_tag_count_filter = self.df_tag_count[self.df_tag_count['count']>self.THRESHOLD_UMI]
        
        # Filter out super high count reads if their probability within the norm distribution are less than P_OUTLIER
        model = norm.fit(self.df_tag_count_filter['count'])
        THRESHOLD_OUTLIER = norm.interval(1-2*P_OUTLIER,model[0],model[1])[1]
        print('select tag counts between %.1f and %.1f (P = %s)'%(self.THRESHOLD_UMI,THRESHOLD_OUTLIER,str(P_OUTLIER)))
        if view == True:
            f,ax = plt.subplots(1,2,figsize = figsize,dpi = dpi,sharex = True,sharey = True)
            hist = ax[0].hist(self.df_tag_count_filter['count'],density = 0,bins = int(max(self.df_tag_count_filter['count'])/50),color = 'red',label = 'Pre-filter',alpha = 0.5)
            ax[0].axvline(x = THRESHOLD_OUTLIER,ymax=1,linestyle = '--')
            self.df_tag_count_filter = self.df_tag_count_filter[self.df_tag_count_filter['count'] < THRESHOLD_OUTLIER]
            ax[1].hist(self.df_tag_count_filter['count'],density = 0,bins = int(max(self.df_tag_count_filter['count'])/50),color = 'green',label = 'Outlier removed',alpha = 0.5)
            ax[0].set_yscale('log')
            ax[0].set_xlabel('UMI counts per barcode')
            ax[1].set_xlabel('UMI counts per barcode')
            ax[0].set_ylabel('Frequency')
            ax[0].legend()
            ax[1].legend()
        else:
            self.df_tag_count_filter = self.df_tag_count_filter[self.df_tag_count_filter['count'] < THRESHOLD_OUTLIER]
               
    def get_tag_counts_table(self,filter_barcodes = True, overwrite = False):
        "Get UMI counts for each antibody tag and transform vector to m by n matrix"
        "This step takes a lot of time if counting all barcodes 'filter_barcodes = False' "
        
        # Check if file has already been processed
        if os.path.isfile('cache/%s/UMI_Matrix.pkl'%self.PROJ_NAME) and overwrite == False:
            self.df_merge_UMI_mx = pd.read_pickle('cache/%s/UMI_Matrix.pkl'%self.PROJ_NAME)    
        else:
            df_merge_UMI = self.df_merge.groupby(['cell_code','tag_code']).size().reset_index(name = 'count')
            df_merge_UMI.index = df_merge_UMI['cell_code']+df_merge_UMI['tag_code']
            self.df_merge_UMI_mx = pd.DataFrame(0,columns=self.AB_list,index = set(df_merge_UMI['cell_code']))
            if filter_barcodes:
                self.df_merge_UMI_mx = self.df_merge_UMI_mx.loc[self.df_tag_count_filter.index]
            for cell in tqdm(list(self.df_merge_UMI_mx.index)):
                for tag in self.AB_list:
                    self.df_merge_UMI_mx[tag].loc[cell] = match_code(cell,tag,df_merge_UMI)
            # Temp file to save time for future use
            self.df_merge_UMI_mx.columns = [self.dict_tag[col] for col in self.df_merge_UMI_mx.columns]
            self.df_merge_UMI_mx.to_pickle('cache/%s/UMI_Matrix.pkl'%self.PROJ_NAME)
    
    def qc_tag_label(self,max_population = 3, min_population_size = 0.2, confidence_interval = 0.90, verbose = False, ncol =3 ,dpi = 100,view = True):
        '''
        Convert tag UMI counts to labels based on stats results
        max_population: Define the maximal number of populations exist in sample datasets
        min_population_size: The smallest population should have at least 20% population
        confidence_interval: if unimodel was used, select the confidence interval for lower bound; 0.90 = 5% confidence one tail test
        '''
        
        if view:
            f,ax = plt.subplots(ceil(len(self.AB_list)/ncol),ncol,figsize = (min(16,ncol*5),3*ceil(len(self.AB_list)/ncol)),dpi = dpi)
            ax = ax.ravel()
            plt.tight_layout(pad = 1.2)
            
        for i,id_AB in enumerate(self.AB_list):
            data = self.df_merge_UMI_mx[self.dict_tag[id_AB]]
            data = data[data!=0]
            data = data.values.reshape(-1, 1)
            threshold,model_kind,best_mdoel,best_population = Threshold_finder(data, max_population = max_population, min_population_size = min_population_size, confidence_interval = confidence_interval, verbose = verbose)
            self.dict_threshold[self.dict_tag[id_AB]] = threshold
            
            if view:    
                hist = ax[i].hist(data, bins=min(300,len(data)),density = 1,alpha = 0.6,color = 'grey')
                ymax_plot = 0 # ymax - Get max y value for plotting
                
                if model_kind == 'Rayleigh' or model_kind == 'Chi-square':
                    ax[i].plot(hist[1],best_mdoel.pdf(hist[1]), label = model_kind);
                    ymax_plot = best_mdoel.interval(0.99)[1]
                else:             
                    for p in np.arange(best_population):
                        para = norm.fit(mask_list(data,best_mdoel.predict(data),p))
                        sub_group = Hawk_smash(mask_list(data,best_mdoel.predict(data),p))
                        ax[i].plot(hist[1],norm.pdf(hist[1],para[0],para[1])*best_mdoel.weights_[p], label='Gaussian %s'%(p+1));
                        ymax_plot = max(ymax_plot,norm.interval(0.99,para[0],para[1])[1])         
                ax[i].axvline(x = threshold,ymax=1,color = 'b',linestyle = '--') 
                ax[i].text(0.6,0.5,'Threshold = '+str(threshold),transform=ax[i].transAxes)
                ax[i].legend();
                if best_population > 1:
                    ax[i].set_ylim([0,max(hist[0][5:]*3)]);
                ax[i].set_xlim([0,ymax_plot*1.2]);
                ax[i].set_title(self.dict_tag[id_AB])
                       
    def get_tag_label(self):
        self.df_tag_label = self.df_merge_UMI_mx[[]].copy()
        for id_AB in self.dict_threshold.keys():
            # Only if the minimal threshold is > 1, the AB signal should be used
            if self.dict_threshold[id_AB] > 1:
                self.df_tag_label[id_AB] = self.df_merge_UMI_mx[id_AB].apply(lambda x: 1 if x>= self.dict_threshold[id_AB] else 0)     
        self.df_tag_label.to_pickle('./cache/%s/df_tag_label.pkl'%self.PROJ_NAME)
        self.list_barcodes = list(self.df_tag_label.index)

class CITE_count:
    '''
    Merge pre-processed AB and RNA data and conduct down stream processes
    '''
    def __init__(self,obj_AB,obj_RNA):
        if obj_AB.data_type == 'AB' and obj_RNA.data_type == 'RNA':
            self.obj_AB = obj_AB
            self.obj_RNA = obj_RNA
        else:
            print('Warning: Wrong data type, please check your preprocessed dataset!')
        
    def merge_data(self,view = True,dpi = 100):    
        self.cell_index = set(self.obj_AB.list_barcodes).intersection(self.obj_RNA.list_barcodes)
        self.df_AB_tag = self.obj_AB.df_tag_label.loc[self.cell_index]
        self.df_AB_count = self.obj_AB.df_merge_UMI_mx.loc[self.cell_index]
        self.df_RNA_count = self.obj_RNA.df_RNA_count.loc[self.cell_index]
        
            
        if view == True:
            plt.subplots(1,1,dpi = dpi)
            venn2([set(self.obj_AB.list_barcodes), set(self.obj_RNA.list_barcodes)], set_labels = ('AB_barcodes', 'RNA_barcodes'))
    
    def group_aliquot(self,min_cell_number = 100,min_portion = 0.5,max_depth = 5, min_dist = 0.3):
    
        df_CITE_All = self.df_AB_tag
        id_node = 0
        # Store all AB_Node information within this list
        self.Obj_tree = []
        for depth in np.arange(max_depth):
            if depth == 0:
                self.Obj_tree.append([AB_Node('labeled_cells ID[%s]'%str(id_node),list_cell = list(df_CITE_All.index),group_id = id_node)])
                id_node+=1
                branch = True
            
            elif branch == True: # group can still be divided
                branch = False
                list_obj_temp = []
                for group in self.Obj_tree[depth-1]:
                    if len(group.list_cell) >= min_cell_number:
                        df_CITE = df_CITE_All.drop(columns=list(set(Hawk_smash(get_branch(group,self.Obj_tree)))))
                        try:
                            # Select only cells within this subgroup
                            df_CITE = df_CITE.loc[group.list_cell]
                            df_dist = group_cell(df_CITE)
                            
                            for rank in np.arange(int((len(df_dist)**2-len(df_dist))/2)-1):
                                p1,p2,score = get_pair(df_dist,rank)
                                A = set(df_CITE[df_CITE[p1]==1].index)
                                B = set(df_CITE[df_CITE[p2]==1].index)
                                portion = len(set(A).union(set(B)))/len(df_CITE)
                                if score > min_dist and p1!=p2 and portion >= min_portion:
                                    # Get p1 positive only
                                    p1_pos = list(df_CITE.loc[set(df_CITE[df_CITE[p1]==1].index)-set(df_CITE[df_CITE[p2]==1].index)].index)
                                    # Get p2 positive only
                                    p2_pos = list(df_CITE.loc[set(df_CITE[df_CITE[p2]==1].index)-set(df_CITE[df_CITE[p1]==1].index)].index)
                                    # Get dual positive
                                    dual_pos = list(df_CITE.loc[set(df_CITE[df_CITE[p2]==1].index).intersection(df_CITE[df_CITE[p1]==1].index)].index)
                                    # Get dual negative
                                    dual_neg = list(df_CITE.loc[set(df_CITE[df_CITE[p2]==0].index).intersection(df_CITE[df_CITE[p1]==0].index)].index)

                                    if p1 in self.obj_AB.high_signal_list:
                                        symbol1 = ['high','low']
                                    else:
                                        symbol1 = ['+','-']
                                        
                                    if p2 in self.obj_AB.high_signal_list:
                                        symbol2 = ['high','low']
                                    else:
                                        symbol2 = ['+','-']
                                    
                                    # Add all Nodes to list
                                    list_obj_temp.append(AB_Node(p1+symbol1[0]+'/'+p2+symbol2[1]+' ID[%s]'%str(id_node),parent = group.Node,list_cell = p1_pos,pair = [p1,p2],group_id = id_node))
                                    id_node+=1
                                    list_obj_temp.append(AB_Node(p1+symbol1[1]+'/'+p2+symbol2[0]+' ID[%s]'%str(id_node),parent = group.Node,list_cell = p2_pos,pair = [p1,p2],group_id = id_node))
                                    id_node+=1
                                    list_obj_temp.append(AB_Node(p1+symbol1[0]+'/'+p2+symbol2[0]+' ID[%s]'%str(id_node),parent = group.Node,list_cell = dual_pos,pair = [p1,p2],group_id = id_node))
                                    id_node+=1
                                    list_obj_temp.append(AB_Node(p1+symbol1[1]+'/'+p2+symbol2[1]+' ID[%s]'%str(id_node),parent = group.Node,list_cell = dual_neg,pair = [p1,p2],group_id = id_node))
                                    id_node+=1
                                    break
                        except (ZeroDivisionError,KeyError,IndexError) as e:
                            print('Warning: Error detected for %s %s at rank %s'%(p1,p2,rank))
                            continue
                            
                if len(list_obj_temp)> 0:
                    branch = True 
                    self.Obj_tree.append(list_obj_temp)
                else:
                    brach = False
                    pass


    def plot_tree(self,export = None):
        for pre, fill, node in RenderTree(self.Obj_tree[0][0].Node):
            print("%s%s" % (pre, node.name))
        if export != None:
            DotExporter(self.Obj_tree[0][0].Node).to_dotfile("%s.dot"%export)
            render('dot', 'png', '%s.dot'%export)
    
    def get_all_groups(self):
        return [obj.name for obj in Hawk_smash(self.Obj_tree)]
    
    def get_all_leaves(self):
        return [obj.name for obj in Hawk_smash(self.Obj_tree) if obj.Node.is_leaf == True]
    
    def get_group_by_depth(self,depth = 0):
        return [obj.name for obj in Hawk_smash(self.Obj_tree) if obj.Node.depth == depth]
    
    def get_group_by_ID(self,list_ID = []):
        return [obj.name for obj in Hawk_smash(self.Obj_tree) if obj.group_id in list_ID]
            
    def get_group_barcode(self,name):
        for obj in Hawk_smash(self.Obj_tree):
            if obj.name == name:
                return obj.list_cell
    
    def qc_results(self,bins = 100,dpi = 100):
        umis_per_cell = self.df_RNA_count
        f,ax = plt.subplots(1,2,figsize = (10,3),dpi = 80)
        ax[0].hist(np.log10(umis_per_cell.sum(axis = 1)), bins=bins)
        ax[0].set_xlabel('UMI counts per cell (log10)')
        ax[0].set_ylabel('Frequency')
        ax[0].set_title('UMI Distribution')
        ax[1].hist(np.log10((umis_per_cell>0).sum(axis = 1)), bins=bins)
        ax[1].set_xlabel('gene counts per cell (log10)')
        ax[1].set_title('gene counts Distribution')
        ax[1].set_ylabel('Frequency')
        
    
    def export_to_scanpy(self,tag_export = 'count'):
        '''
        Export to scanpy object
        tag_export = 'count' - export tag counts
        tag_export = 'label' - export tag label, 0 for low and 1 for high
        
        '''
        sc_input = self.df_RNA_count    
        sc_export = sc.AnnData(X = sc_input.values)
        sc_export.obs.index = sc_input.index
        if tag_export == 'count':
            sc_export.obs = self.df_AB_count.loc[sc_input.index]
        elif tag_export == 'label':
            sc_export.obs = self.df_AB_tag.loc[sc_input.index]
        sc_export.var.index = sc_input.columns
        return sc_export
        
    def PCA(self,group = None,n_components=200, view = False):
        sc_input = self.df_RNA_count   
        self.sc = sc.AnnData(X = sc_input.values)
        self.sc.obs.index = sc_input.index
        self.sc.var.index = sc_input.columns
        sc.tl.pca(self.sc, n_comps=n_components)
        self.PCA_result = self.sc.obsm['X_pca']
        if view == True:
            plt.plot(self.sc.uns['pca']['variance_ratio']*100)
            plt.xlabel('Principle components')
            plt.ylabel('Variance %')
            plt.title('PCA for group %s'%str(group))
    
    def umap(self,n_components = 10,n_neighbors=5, min_dist=0.3, metric='correlation'):
        np.random.seed(0)
        
        self.embedding = umap.UMAP(n_neighbors=n_neighbors,
                      min_dist=min_dist,
                      metric=metric).fit_transform(self.PCA_result[:,:n_components])
    
    def umap_plot(self,label = None,kind = 'tree',color_scheme = 'Paired',marker_size = 2,dpi = 100):
        self.color = []
        self.legend_elements = []
        
        self.group_norm = self.get_group_by_ID([0]) # This could be edited in the future for subsetting cells
        
        embedding = self.embedding
        
        if kind == 'tree':
            if label == None:
                group_label = self.get_all_leaves()
            else:
                group_label = label
            # Color code all labels, no color for groups not selected for labeling
            dict_label_color = dict(zip(group_label,sns.color_palette(color_scheme, len(group_label))))
            for group in self.get_all_groups():
                if group not in group_label:
                    dict_label_color[group] = (1.0,1.0,1.0)


            for pre, fill, node in RenderTree(self.Obj_tree[0][0].Node):
                self.legend_elements.append(Line2D([0], [0], marker='o',color = 'w', markerfacecolor = dict_label_color[node.name.split('(')[0][:-1]], label="%s%s" % (pre, node.name)))        

            # Color code cells        
            dict_cell_color = dict()
            for group in group_label:
                c = dict_label_color[group]
                dict_cell_color.update(dict(zip(self.get_group_barcode(group),[c for xx in np.arange(len(self.get_group_barcode(group)))])))

            for i,group in enumerate(self.group_norm):
                df_group = self.df_RNA_count.loc[self.get_group_barcode(group)]
                for cell in self.get_group_barcode(group):
                    try:
                        self.color.append(dict_cell_color[cell])
                    except KeyError:    
                        self.color.append((0.9,0.9,0.9))        
                        
            plt.subplots(1,1,figsize = (8,8),dpi = dpi)
            plt.scatter(embedding[:, 0], embedding[:, 1],marker = '.',s = marker_size,c = self.color)
            plt.gca().set_aspect('equal', 'datalim')
            plt.title('UMAP projection with antibody labels', fontsize=18);
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,handles=self.legend_elements,fontsize = min(16,250/len(self.get_all_groups())))
        
        elif kind == 'gene':
            dict_RNA_selection_log = dict(self.df_RNA_count[label].apply(lambda x: np.log(x+1)))
            for group in self.group_norm:
                df_group = self.df_RNA_count.loc[self.get_group_barcode(group)]
                for cell in self.get_group_barcode(group):
                    try:
                        self.color.append(dict_RNA_selection_log[cell])
                    except KeyError:    
                        self.color.append(0)

            plt.subplots(1,1,figsize = (10,8),dpi = dpi)
            plt.scatter(embedding[:, 0], embedding[:, 1],marker = '.',s = marker_size,c = self.color, cmap= color_scheme)
            plt.gca().set_aspect('equal', 'datalim')
            plt.title('UMAP projection with antibody labels', fontsize=18);
            plt.colorbar()
            
        elif kind == 'antibody':
            dict_AB_selection = dict(self.df_AB_count[label].apply(lambda x: np.log(x+1)))
            
            for group in self.group_norm:
                df_group = self.df_RNA_count.loc[self.get_group_barcode(group)]
                for cell in self.get_group_barcode(group):
                    try:
                        self.color.append(dict_AB_selection[cell])
                    except KeyError:    
                        self.color.append(0)

            plt.subplots(1,1,figsize = (10,8),dpi = dpi)
            plt.scatter(embedding[:, 0], embedding[:, 1],marker = '.',s = marker_size,c = self.color, cmap= color_scheme)
            plt.gca().set_aspect('equal', 'datalim')
            plt.title('UMAP projection with antibody labels', fontsize=18);
            plt.colorbar()