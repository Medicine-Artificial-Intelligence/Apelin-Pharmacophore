import pandas as pd
import glob
import os
import numpy as np
def data_process(path):
    os.chdir(path)
    data_dir = str(os.getcwd())
    data_name = []
    for i in sorted(glob.glob(f"{path}/*.txt")):
        data_name.append(i[len(data_dir)+1:-4])
    print(data_name)
    
    
    table = pd.DataFrame()
    for i in data_name:
        data = pd.read_csv(f"{i}.txt")
        post =Postprocessing.pharmacophore_postprocess(data = data, model =i, ID = 'ID', ref=dude, 
                                 scores = 'rmsd', rescore = 'minmax')
        post.fit()
        #display(post.ref.head())
        post.ref[f'{i}_predict'] = post.ref['predict']
        #post.ref[f'{i}_scores'] = post.ref['scores']
        #display(post.data)
        post.ref = post.ref.drop(['ID', 'Canomicalsmiles','predict', 'scores'], axis =1)
        table = pd.concat([table, post.ref], axis =1)
    table = table.T.drop_duplicates().T
    
   
    display(table)
    return table