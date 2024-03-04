# functions to plot the generated data
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os
import numpy 
import plotly.express as px
import pymatgen.core as mg
import plotly.graph_objects as go


def data(df):
    formula = df['reduced_formula']
    stoichiometry = df['Atoms']
    matches = df['Matches']
    make = df['Make']
    mpid = df['mp-id']
    return formula, stoichiometry, matches, make, mpid



def plot(df, elems, path):
    formula, stoichiometry, matches, make, mpid = data(df)
    num_mp=[]
    num_match = []
    num_make = []
    num_rest = []
    for j in range(len(formula)):
        mpid_temp = mpid[j]
        df_temp = pd.DataFrame.from_dict(mg.Composition.get_el_amt_dict(mg.Composition(formula[j])), orient="index")
        df_temp = (df_temp.T).reset_index()
        del df_temp['index']
        corners = df_temp.columns
      
        if j ==0:
            df_mp = pd.DataFrame(columns=corners)
            df_match = pd.DataFrame(columns=corners)
            df_make = pd.DataFrame(columns=corners) 
            df_rest = pd.DataFrame(columns=corners)
            

        if 'mp-' in mpid_temp:
            df_mp = pd.concat((df_mp,df_temp), ignore_index=True, axis = 0)
            num_mp = numpy.append(num_mp, j)

        elif matches[j]==True and 'mp-' not in mpid_temp:
            df_match = pd.concat((df_match,df_temp), ignore_index=True, axis = 0)
            num_match = numpy.append(num_match, j)

        
        elif make[j] == True and  matches[j]==False:
            df_make = pd.concat((df_make,df_temp), ignore_index=True, axis = 0)
            num_make = numpy.append(num_make, j)

        else : 
            df_rest = pd.concat((df_rest,df_temp), ignore_index=True, axis = 0)
            num_rest = numpy.append(num_rest, j) 
       
    num = numpy.concatenate((num_match, num_make, num_rest, num_mp))
    hover_name = formula[num]
        
    hover_name = hover_name.reset_index()
    del hover_name['index']

    df_mp['Materials'] = 'Pre-existing'
    df_match['Materials'] = 'Missing Known Stoichiometry'
    df_make['Materials'] = 'Materials within range'
    df_rest['Materials'] = 'New Materials which cannot be made'
    

    
    corners = df_rest.columns
    df_all = pd.concat((df_match, df_make, df_rest, df_mp), ignore_index=True)
    df_all['hover_name'] = hover_name
    symbols = ['star', 'star', 'circle-open' , 'square-open']
    colors = ['blue', 'red', 'black', 'green']


    def makeAxis(title, tickangle):
        return {
        'title': title,
        'titlefont': { 'size': 20 },
        'tickangle': tickangle,
        'tickfont': { 'size': 15 },
        'tickcolor': 'rgba(0,0,0,0)',
        'ticklen': 5,
        'showline': True,
        'showgrid': True
        }
    
    fig = px.scatter_ternary(df_all, a = corners[0], b = corners[1], c = corners[2],
                            symbol = df_all['Materials'],
                            color = df_all['Materials'],
                            symbol_sequence= symbols,
                            color_discrete_sequence = colors,
                            width=600, height=400,
                            hover_name=df_all['hover_name'],
                            hover_data={corners[0]:False,
                                        corners[1]:False,
                                        corners[2]:False,
                                        },
                            title=("Ternary compositions of : %s" % corners[0]+corners[1]+corners[2])
                            )
    

    
    fig.update_layout({
    'ternary': {
        'sum': 100,
        'aaxis': makeAxis(corners[0], 0),
        'baxis': makeAxis(corners[1], 45),
        'caxis': makeAxis(corners[2], -45)
    },
    
    })
    
    #fig.show()
       

    
    file = str(elems)+'.html'  
    print(file)
    fig.write_html(path/file,full_html=False)


def main(df, elems, path):
    elems = ''.join(elems)
    path = path
    plot(df, elems, path)


 