import numpy as np
import pandas as pd

def find_elements(elements, df_charge_probability_database):
    recorded_ele_list = pd.Series(df_charge_probability_database['elements'])
    elements_2find = elements[~elements.isin(recorded_ele_list)]
    return elements_2find




def main(elements, charge_mixing):
    charge_probability_database = 'charge_probability_database.csv'
    df_charge_probability_database = pd.read_csv(charge_probability_database)   
    elements_2find = find_elements(elements, df_charge_probability_database)
    if len(elements_2find) > 0:
        print('The following elements are not in the database:',elements_2find)
        

    

if __name__ == '__main__':
    main()