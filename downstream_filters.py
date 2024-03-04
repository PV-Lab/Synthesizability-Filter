# filter generated materials for their presence in the ICSD and materials project databases.
# written by Basita Das
import numpy
import pandas as pd
import pymatgen.core as mg
import itertools



def mpid(df : pd.DataFrame):
    """
    Filters the input DataFrame based on the 'mp-id' column.

    Args:
    df (pandas.DataFrame): The input DataFrame to be filtered.

    Returns:
    pandas.DataFrame: A new DataFrame that only includes rows where 'mp-id' is not equal to '[]' or that mp-id is available.
    """
    df_mpid = df[df['mp-id'] != '[]'].copy().reset_index(drop=True)
    return df_mpid


def mpid_range(df : pd.DataFrame,stoichimetric_spread: float):
    """
    Calculates the maximum and minimum of the range of stoichiometric ratios and
    then adds the stoichimetric spread to create the new range where we 
    will look for new material.

    Args:
    df (pandas.DataFrame): The input DataFrame, expected to have a column named "Atoms".
    stoichimetric_spread (float): The stoichiometric spread to be used in the calculation.

    Returns:
    tuple: A tuple containing two pandas.Series, max_mpid_range and min_mpid_range.
    """
    df = df["Atoms"]
    num_of_rows = df.size
    mpid_range = pd.DataFrame([])
  
    for i in range(num_of_rows):
        row = df[i]
        combinations = list(itertools.combinations(row, 2))
        mpid_range[i] = ([(a / b) for a, b in combinations])
    mpid_range = mpid_range.round(3)
    if num_of_rows == 1:
        max_mpid_range = mpid_range * (1 + stoichimetric_spread/100)
        min_mpid_range = mpid_range * (1 - stoichimetric_spread/100)
        max_mpid_range = max_mpid_range.iloc[:,0] # convert them into series for  data compatibility
        min_mpid_range = min_mpid_range.iloc[:,0] # convert them into series for  data compatibility
    elif num_of_rows > 1:
        max_mpid_range = mpid_range.max(axis = 1) * (1 + stoichimetric_spread/100)
        min_mpid_range = mpid_range.min(axis = 1) * (1 - stoichimetric_spread/100)
    elif num_of_rows == 0:
        max_mpid_range = None
        min_mpid_range = None
    return max_mpid_range, min_mpid_range
    
   
def stoich_ratio(df : pd.DataFrame):
    """
    Calculates the ratio between the stoichiometry of every element present in a compound.
    For example, for a compound with the stoichiometry [1, 2, 3], the ratios would be [1/2, 1/3, 2/3].

    Args:
    df (pandas.DataFrame): The input DataFrame, expected to have a column named "Atoms".

    Returns:
    pandas.DataFrame: A new DataFrame with the calculated stoichiometric ratios.
    """
    df_atoms = df["Atoms"]
    num_of_rows = df_atoms.size
    stoich_range = pd.DataFrame([])
    for i in range(num_of_rows):
        row = df_atoms[i]
        combinations = list(itertools.combinations(row, 2))
        stoich_range[i] = ([(a / b) for a, b in combinations])
    stoich_range = stoich_range.round(3)
    return stoich_range

def in_range(max : pd. Series, min : pd.Series, all_range :pd.DataFrame):
    """
    Checks if the ratio between the stoichimetry of different elements in a compound are within the specified range.
    
    Args:
    max (pandas.Series): The maximum value of the range.
    min (pandas.Series): The minimum value of the range.
    all_range (pandas.DataFrame): The DataFrame containing the stoichiometric ratios.

    Returns:
    tuple: A tuple containing three pandas.DataFrames - max_range, min_range, and in_range. 
           max_range and min_range are the reshaped input max and min, and in_range is a DataFrame 
           indicating whether each stoichiometric ratio is within the specified range.
    """
    cols  = all_range.shape[1]
    max_range = pd.DataFrame(max.repeat(cols)).values.reshape(all_range.shape)
    min_range = pd.DataFrame(min.repeat(cols)).values.reshape(all_range.shape)
    in_range = pd.DataFrame([])
    in_range = all_range.le(max_range) & all_range.ge(min_range)
    return max_range, min_range, in_range

def addtodf(df : pd.DataFrame, new_column :str, new_column_data : pd.Series):
    """
    Adds a new column to the input DataFrame.

    Args:
    df (pandas.DataFrame): The input DataFrame to which the new column will be added.
    new_column (str): The name of the new column.
    new_column_data (pandas.Series, list, or array-like): The data for the new column.

    Returns:
    pandas.DataFrame: The updated DataFrame with the new column added.
    """
    df[new_column] = new_column_data
    return df

def truefalse(in_range : pd.DataFrame):
    """
    Checks if all values in the input DataFrame are True.

    Args:
    in_range (pandas.DataFrame): The input DataFrame to be checked.

    Returns:
    pandas.Series: A Series indicating whether all values in each row of the input DataFrame are True.
    """
    result = in_range.all()
    return result

def match_stoichimetric_combinations(stoichiometry : list ,df: pd.DataFrame):
    """
    Checks which stoichiometries have been seen before in other material systems.

    Args:
    stoichiometry (list): The stoichiometry to be matched.
    df (pandas.DataFrame): The input DataFrame, expected to have a column named "Atoms".

    Returns:
    pandas.Series: A Series indicating whether the atoms in each row of the input DataFrame match the given stoichiometry.
    """
    atoms = pd.Series(df["Atoms"])
    matches = atoms.isin(stoichiometry)
    return matches
    


def savetocsv(df : pd.DataFrame ,path : str ,elems : list):
    """
    Saves the input DataFrame to a CSV file.

    Args:
    df (pandas.DataFrame): The DataFrame to be saved.
    path (str or pathlib.Path): The directory where the CSV file will be saved.
    elems (list): The elements to be included in the filename.

    Returns:
    None
    """
    elems = ''.join(elems)
    filename = 'statespace'+ '_'+ elems + '.csv'
    df.to_csv(path/filename)

def create_df(df : pd.DataFrame ,elems : list, path : str, max_mpid_range : pd. Series, min_mpid_range : pd.Series,
               stoichimetric_ratios : pd.DataFrame, matches : pd.Series):
    """
    Creates a DataFrame with additional columns based on the input parameters.

    Args:
    df (pandas.DataFrame): The input DataFrame to which the new columns will be added.
    elems (list): The elements to be included in the filename.
    path (str or pathlib.Path): The directory where the CSV file will be saved.
    max_mpid_range (pandas.Series or None): The maximum value of the range.
    min_mpid_range (pandas.Series or None): The minimum value of the range.
    stoichimetric_ratios (pandas.DataFrame): The DataFrame containing the stoichiometric ratios.
    matches (pandas.Series or None): A Series indicating whether the atoms in each row of the input DataFrame match the given stoichiometry.

    Returns:
    None
    """
    if max_mpid_range is None:
        df = addtodf(df, 'Stoich x/y x/z x/z ', stoichimetric_ratios.transpose().values.tolist())
        df = addtodf(df, 'Max mpid range', None)
        df = addtodf(df, 'Min mpid range', None)
        df = addtodf(df, 'In range', None)
        df = addtodf(df, 'Make', None)
        df = addtodf(df, 'Matches', None)
    else:
        max_range, min_range, candidates_in_range = in_range(max_mpid_range, min_mpid_range, stoichimetric_ratios)
        make = truefalse(candidates_in_range)
        df = addtodf(df, 'Stoich x/y x/z x/z ', stoichimetric_ratios.transpose().values.tolist())
        df = addtodf(df, 'Max mpid range', max_range.transpose().tolist())
        df = addtodf(df, 'Min mpid range', min_range.transpose().tolist())
        df = addtodf(df, 'In range', candidates_in_range.transpose().values.tolist())
        df = addtodf(df, 'Make', make)
        df = addtodf(df, 'Matches', matches)
    savetocsv(df, path, elems)

def stoichiometry_main(df : pd.DataFrame, path : str, elems : list ,stoichimetric_spread : float,stoichiometry : list):
    """
    Main function to filter the DataFrame based on stoichiometry and save the result to a CSV file.

    Args:
    df (pandas.DataFrame): The input DataFrame, expected to have a column named "Atoms".
    path (str or pathlib.Path): The directory where the CSV file will be saved.
    elems (list): The elements to be included in the filename.
    stoichimetric_spread (float): The spread for the stoichiometric ratios.
    stoichiometry (list): The stoichiometry to be matched.

    Returns:
    tuple: A tuple containing two pandas.DataFrames - the updated DataFrame and the 'Atoms' column of df_mpid.
    """
    matches = match_stoichimetric_combinations(stoichiometry,df)
    df_mpid = mpid(df)
    max_mpid_range, min_mpid_range = mpid_range(df_mpid,stoichimetric_spread)
    stoichimetric_ratios = stoich_ratio(df)
    create_df(df,elems, path, max_mpid_range, min_mpid_range, stoichimetric_ratios,matches)
    return df, df_mpid['Atoms']
    


    