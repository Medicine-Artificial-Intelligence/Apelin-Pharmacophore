import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from rdkit.Chem import AllChem, Descriptors
from rdkit import DataStructs

class cluster_scatter_plot:
    """""
    Input:
    - df: dataframe
        must have Molecules column and cluster index
    - no_cls: int
        number of cluster in scatter plot
    - mol_col: string
        name of column containg Molecules
    - cluster_col: string
        name of column containg cluster index
    - algo: string
        name of clustering algorithm
    - radius: int (2)
        ECFP value
    - nBits: int (2048)
        ECFP value
    Return:
    - scatter plot
    """""
    def __init__(self,data, no_cls, mol_col, cluster_col, algo, radius = 2, nBits = 2048,palette=sns.color_palette("husl", 8)):
        self.data = data
        self.no_cls = no_cls
        self.mol_col = mol_col
        self.cluster_col = cluster_col
        self.algo = algo
        self.radius=radius
        self.nBits = nBits
        self.palette = palette
    
    def mol2ecfp(self, mol):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=self.radius, nBits = self.nBits)
        ar = np.zeros((1,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, ar)
        return ar
    
    def processing(self):
        self.df_cls = self.data[[self.mol_col, self.cluster_col]]
        idx = self.df_cls[self.df_cls[self.cluster_col] > self.no_cls].index
        self.df_cls.drop(idx, axis = 0, inplace = True)
        self.df_cls.reset_index(drop = True, inplace = True)
        
        # fp
        self.df_cls["FPs"] = self.df_cls[self.mol_col].apply(self.mol2ecfp)
        X = np.stack(self.df_cls.FPs.values)
        X_df = pd.DataFrame(X)
        self.df_cls= pd.concat([self.df_cls, X_df], axis = 1).drop(["FPs", self.mol_col], axis =1)
        
        # tSNE
        pca  = PCA(n_components=50)
        pca_components = pca.fit_transform(self.df_cls.drop(self.cluster_col, axis = 1))
        tSNE=TSNE(n_components=2, random_state = 42)
        tSNE_result= tSNE.fit_transform(pca_components)
        self.x=tSNE_result[:,0]
        self.y=tSNE_result[:,1]
        
    def visualize(self):
        self.processing()
        sns.set_theme(font_scale=1, context ='notebook',style='darkgrid',)
        
        self.df_cls['x']=self.x
        self.df_cls['y']=self.y       

        fig =plt.figure(figsize=(10,8))
        fig =sns.scatterplot(x='x',y='y',hue = 'Cluster',palette=self.palette,
                             data=self.df_cls,legend="full")
        fig.set_title(f'{self.algo}', fontsize = 24, weight = 'semibold')
        fig.set_xlabel("tSNE-1", fontsize=16)
        fig.set_ylabel("tSNE-2", fontsize=16)
        fig.legend(loc='upper left', bbox_to_anchor=(1.00, 0.75), ncol=1);
        plt.savefig(f"{self.algo}.png", dpi = 600)
        
        
class cluster_heat_map:
    """""
    Input:
    - cls_cps: list
        list of molecules, must have name (GetProp('_Name'))
    - radius: int (2)
        ECFP value
    - nBits: int (2048)
        ECFP value
    Return:
    - 
    """""
    def __init__(self,cls_cps, radius = 2, nBits = 2048):
        self.cls_cps = cls_cps
        self.radius=radius
        self.nBits = nBits
        sns.set(font_scale=2)
        self.fig = plt.figure(figsize = (20,20))
        self.processing()
        
    def processing(self):
        fps= [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 2048) for mol in self.cls_cps]
        size=len(self.cls_cps)
        hmap=np.empty(shape=(size,size))
        self.table=pd.DataFrame()
        for index, i in enumerate(fps):
            for jndex, j in enumerate(fps):
                similarity=DataStructs.FingerprintSimilarity(i,j, metric=DataStructs.TanimotoSimilarity)
                hmap[index,jndex]=similarity
                self.table.loc[self.cls_cps[index].GetProp('_Name'),self.cls_cps[jndex].GetProp('_Name')]=similarity
    
    def visualize_triangle(self): 
        
        corr = self.table.corr()
        mask = np.tril(np.ones_like(corr, dtype=bool))
            # generating the plot
        
        self.fig = sns.heatmap(self.table, annot = True, annot_kws={"fontsize":15}, center=0,
                    square=True,  linewidths=.1, cbar_kws={"shrink": .1}, mask = mask)

        plt.title('Heatmap of Tanimoto Similarities', fontsize = 25 )
        
        
    def visualize_square(self):
            # generating the plot
         
        self.fig = sns.heatmap(self.table, annot = True, annot_kws={"fontsize":10}, center=0,
                    square=True,  linewidths=.7, cbar_kws={"shrink": .5}, vmin = 0, vmax = 1)

        plt.title('Heatmap of Tanimoto Similarities', fontsize = 24) # title with fontsize 20
    
