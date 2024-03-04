import numpy
from element_focused_CS_analyzer_part1 import filter
import main
import pandas as pd
import itertools
import downstream_filters
from pathlib import Path
import glob
import data_analysis
import os
from htmlmerger import HtmlMerger
from pathlib import Path

A = ['Cs', 'K', 'Na', 'Rb']
B = ['In', 'Sn', 'Sb', 'Pb', 'Bi']
#C = ['Ag', 'In', 'Sn', 'Sb', 'Pb', 'Bi']
X = ['Cl','I','Br']

charge_mixing = False
charge_frequency = 0.3
stoichimetric_spread = 20 #percentage
dirname = 'Results/Ternary_statespace_spread_%s' % stoichimetric_spread
cur_dir = Path(os. getcwd())
path = Path(cur_dir/dirname)


if not os.path.exists(path):
    os.makedirs(path)

# Existing Stoichioemtry
stoichiometry = [[3, 1, 6],[3, 2, 9],[1, 1, 4],[2, 1, 6],[1, 1, 3],[1, 2, 5],[4, 1, 6],[1, 2, 7],[2, 1, 5],[3, 1, 5]]

combinations = pd.Series(itertools.product(A, B, X))
filtered_combinations = ([combination for combination in combinations if combination[:][1] != combination[:][2]])
idx = len(filtered_combinations)

df_combined = pd.DataFrame(['composition','Atoms','reduced_formula','mp-id', 'icsd_id', 'Make', 'Matches'])
df_mpid_all = []
for i in range(idx):
    elems = (filtered_combinations[i])
    num = len(elems)        # number of elements in the compound
    df = filter(elems, charge_mixing, num, charge_frequency, path)
    df, df_mpid = downstream_filters.stoichiometry_main(df, path, elems,stoichimetric_spread,stoichiometry)
    df_subset = df[['composition','Atoms','reduced_formula','mp-id', 'icsd_id', 'Make', 'Matches']].copy()
    df_combined = pd.concat([df_subset, df_combined], ignore_index=True)
    data_analysis.main(df,elems, path)
    #df_mpid_all.append(df_mpid) 
#df_mpid_all = df_mpid_all.drop_duplicates()
#df_mpid_all = df_mpid_all.reset_index(drop=True)
#df_mpid_all.to_csv(path/'df_mpid_all.csv')
df_combined.to_csv(path/'df_all.csv')


merger = HtmlMerger(files=path.glob("*.html"), output_file=path/'merged.html')  # result will be in ./merged.html
merger.merge()
