# Purpose: This code is used to generate charge state probability for a given ternary system.

import pandas as pd
import pymatgen.core as mg
import numpy as np
import re
import math
import pickle
import os
import json
from itertools import combinations_with_replacement as CWR
from itertools import product as pdt
from statistics import mean
import torch

def save_to_csv(attached_data2save, elements_2find, path):
    # filename = 'Summary_ternary'
    # for i in elements_2find:
    #     filename = filename + '_' + i
    # file_name = filename + '.csv'
    # attached_data2save.to_csv("./Results/Results_Eg_sites=20/%s" % file_name)
    # print('Saved ', file_name)
    elems = ''.join(elements_2find)
    filename = 'statespace'+ '_'+ elems + '.csv'
    attached_data2save.to_csv(path/filename)
    print('-'*30)

def create_final_dataframe(elements_2find, cp, charges_splits, formula_list, pretty_formula, mp_ids, icsd_ids, mp_FEPA, mp_eabovehull, d2s):
    charge_breakdown_ind_ele = {}
    for ind, ele in enumerate(elements_2find):
        for chg in cp[ind]:
            for ind2, chg2comp in enumerate(charges_splits):
                if ind2 == 0:
                    charge_breakdown_ind_ele[ele+str(chg)] = []
                if chg in chg2comp[ind]:
                    charge_breakdown_ind_ele[ele+str(chg)].append(True)
                else:
                    charge_breakdown_ind_ele[ele+str(chg)].append(False)

    data_to_add = { 'formula':formula_list,
                    'reduced_formula': pretty_formula,
                    'mp-id': mp_ids,
                    'icsd_id': icsd_ids,
                    'mp_FEPA': mp_FEPA,
                    'mp_eabovehull': mp_eabovehull}

    charge_breakdown_ind_ele = pd.DataFrame(charge_breakdown_ind_ele)
    data_to_add = pd.DataFrame(data_to_add)

    frames = [d2s, charge_breakdown_ind_ele, data_to_add ]
    attached_data2save = pd.concat(frames, axis=1)

    return attached_data2save

def match_formulas_with_database(formula_list, df_filtered):
    icsd_ids = ['[]' for _ in range(len(formula_list))]
    mp_ids = ['[]' for _ in range(len(formula_list))]
    mp_FEPA = ['' for _ in range(len(formula_list))]
    mp_eabovehull = ['NA' for _ in range(len(formula_list))]
    pretty_formula = []

    for ind, comp in enumerate(formula_list):
        comp = mg.Composition(comp).reduced_formula
        for ind2, comp2comp in enumerate(df_filtered['pretty_formula']):
            comp2comp = mg.Composition(comp2comp).reduced_formula
            if comp == comp2comp:
                if df_filtered['icsd_ids'][ind2]:
                    icsd_ids[ind]= df_filtered['icsd_ids'][ind2]
                    #icsd_ids[ind].append(df_filtered['icsd_ids'][ind2])
                #mp_ids[ind].append(df_filtered.index[ind2])
                #mp_FEPA[ind] = df_filtered['formation_energy_per_atom'][ind2]
                #mp_eabovehull[ind] = df_filtered['e_above_hull'][ind2]
                mp_ids[ind] = df_filtered.index[ind2]
                mp_FEPA[ind] = df_filtered['formation_energy_per_atom'][ind2]
                mp_eabovehull[ind] = df_filtered['e_above_hull'][ind2]
                    

        pretty_formula.append(comp)

    return icsd_ids, mp_ids, mp_FEPA, mp_eabovehull, pretty_formula

def filter_dataframe(df_all, elements_2find, d2s):
    formula_list = list(d2s["composition"])

    df_filtered = df_all[['formation_energy_per_atom', 
            'unit_cell_formula', 'pretty_formula', 
            'elements', 'nelements', 'e_above_hull', 
            'icsd_ids' ]]
    df_filtered = df_filtered[df_filtered['nelements']==3]

    x1 = [elements_2find[0] in x for x in df_filtered['elements']]
    x2 = [elements_2find[1] in x for x in df_filtered['elements']]
    x3 = [elements_2find[2] in x for x in df_filtered['elements']]
    #x4 = [elements_2find[3] in x for x in df_filtered['elements']]
    df_filtered = df_filtered[[a and b and c for a,b,c  in zip(x1,x2,x3)]]
    
    return formula_list, df_filtered

def create_dataframe(elements_2find, filtered_atom_counts, filtered_charge_comb, charges_splits, ele_charge_counts, ischargemixed, raw_perc_for_each_charg, perc_calculated, one_number_probability):
    name_comp = []
    for i in filtered_atom_counts:
        name_ind = ''
        for ind, j in enumerate(i):
            name_ind = name_ind + elements_2find[ind] + str(j)
        name_comp.append(name_ind)

    d2s = pd.DataFrame ({'composition' : name_comp, 'Atoms' : filtered_atom_counts, 
                        'All Charges': filtered_charge_comb, 'Ind Charges': charges_splits,
                        elements_2find[0] + '(Uniq Charges)': [x[0] for x in ele_charge_counts],
                        elements_2find[1] + '(Uniq Charges)': [x[1] for x in ele_charge_counts],
                        elements_2find[2] + '(Uniq Charges)': [x[2] for x in ele_charge_counts],
                        #elements_2find[3] + '(Uniq Charges)': [x[3] for x in ele_charge_counts],
                        'Is c_mixed': ischargemixed,
                        'Raw_percentage' : raw_perc_for_each_charg, 
                        'Charge_prob_per_ele' : perc_calculated , 'Charge_prob': one_number_probability,
                        'normalized_charge_prob' : one_number_probability/max(one_number_probability)})
    return d2s

def unique(list):
    unique_list = []
    for x in list:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

def calculate_charge_probabilities(filtered_atom_counts, filtered_charge_comb, cp_perc, cp):
    one_number_probability = []
    perc_calculated = []
    raw_perc_for_each_charg = []
    charges_splits = []
    ischargemixed = []
    ele_charge_counts = []

    for ind, i in enumerate(filtered_atom_counts):
        one_formula = []
        one_formula_raw = []
        one_formula_charge = []
        one_ele_charge_counts = []
        c = 0
        charge_mix = False
        for indc,j in enumerate(i):
            one_formula_average = []
            one_ele_charges = []
            for k in range(j):
                temp_c = filtered_charge_comb[ind][c]
                ind_perc = cp_perc[indc][cp[indc].index(int(temp_c))]
                one_formula_average.append(ind_perc)
                one_formula_raw.append(ind_perc)
                one_ele_charges.append(temp_c)
                c+=1      
            unique_one_ele_charges = []
            for chg in one_ele_charges:
                if not chg in unique_one_ele_charges:
                    unique_one_ele_charges.append(chg)
            one_formula.append(mean(one_formula_average))     
            one_formula_charge.append(unique_one_ele_charges)  
            one_ele_charge_counts.append(len(unique_one_ele_charges))
            if len(unique_one_ele_charges) > 1:
                charge_mix = True
        perc_calculated.append(one_formula)
        raw_perc_for_each_charg.append(one_formula_raw)
        one_number_probability.append(np.product(one_formula))
        charges_splits.append(one_formula_charge)
        ele_charge_counts.append(one_ele_charge_counts)
        ischargemixed.append(charge_mix)

    return one_number_probability, perc_calculated, raw_perc_for_each_charg, charges_splits, ischargemixed, ele_charge_counts

def filter_charge_combinations(valid_comb, cp, charge_mixing):
    filtered_charge_comb = []
    filtered_atom_counts = []

    for indc, Nall in enumerate(valid_comb):
        previous_cbn = []
        Nall = list(Nall)
        
        for ind, acharges_of_ele in enumerate(cp):
            #if charge mixing is allowed, compute combination with replacement 
            if charge_mixing:
                combination = list(CWR(acharges_of_ele, int(Nall[ind])))
            else:
                temp_list_c = []
                for k in acharges_of_ele:
                    temp_list_c.append(list([k])*int(Nall[ind]))
                combination = temp_list_c
            
            if ind == 0:
                previous_cbn = combination
            else:
                t_comb_tmp = []
                for j in previous_cbn:
                    for k in combination:
                        temp = j + k
                        t_comb_tmp.append(temp)
                previous_cbn = t_comb_tmp       
        
        #charge neutrality filter
        comb2remove = [(sum(i) == 0) for i in previous_cbn]
        previous_cbn = [i for (i,remove) in zip(previous_cbn, comb2remove) if remove]

        filtered_charge_comb = filtered_charge_comb + previous_cbn
        for j in previous_cbn:
            filtered_atom_counts.append(Nall)

    return filtered_charge_comb, filtered_atom_counts

def filter_charge_polarity(elen, cp, cp_perc):
    if not 0 in elen:
        mini = elen.index(min(elen))
        index = [i > 0 for i in cp[mini]]
        cp[mini] = [i for (i, remove) in zip(cp[mini], index) if remove]
        cp_perc[mini] = [i for (i, remove) in zip(cp_perc[mini], index) if remove]
    
        maxi = elen.index(max(elen))
        index = [i < 0 for i in cp[maxi]]
        cp[maxi] = [i for (i, remove) in zip(cp[maxi], index) if remove]
        cp_perc[maxi] = [i for (i, remove) in zip(cp_perc[maxi], index) if remove]
    return cp, cp_perc

def filter_charges_by_cutoff(elements_2find, ele_symbol_charge, ele_symbol_percentage, cut_off_freq_percentage):
    cp = []
    cp_perc = []
    for ind, acharges_of_ele in enumerate(elements_2find):
        chargs = []
        chargs_per = []
        for indc, j in enumerate(ele_symbol_percentage[ind]):
            if j >= cut_off_freq_percentage:
                chargs.append(ele_symbol_charge[ind][indc])
                chargs_per.append(ele_symbol_percentage[ind][indc])
        cp.append(chargs)
        cp_perc.append(chargs_per)
    return cp, cp_perc

def filter_formulas_by_element_charge(pele, pcharge, round1_ele_charg_dict_list, pformula_r1, r1_o1index, round1_oindex_list):
    temp_formula = []
    for ind, i in enumerate(round1_ele_charg_dict_list):
        if (pele in i) and (i[pele] == pcharge):
            temp_formula.append(pformula_r1[r1_o1index.index(round1_oindex_list[ind])])
            print(pformula_r1[r1_o1index.index(round1_oindex_list[ind])], i)
    return temp_formula

def update_charge_probability_database(elements_2find, recorded_ele_list, ele_symbol_charge, ele_symbol_percentage, elements_2find_original):
    recorded_charge_states = []
    recorded_charge_probabilities = []

    if len(elements_2find):
        for ind, acharges_of_ele in enumerate(elements_2find):
            if not acharges_of_ele in recorded_ele_list:
                recorded_ele_list.append(acharges_of_ele)
                recorded_charge_states.append(str(ele_symbol_charge[ind]))
                recorded_charge_probabilities.append(str(ele_symbol_percentage[ind]))

        data_to_add = {'elements':recorded_ele_list,
                        'charge_states': recorded_charge_states,
                        'probabilities': recorded_charge_probabilities}

        data2save = pd.DataFrame(data_to_add)

        file_name = 'charge_probability_database.csv'
        data2save.to_csv(file_name)
        print('updated ', file_name)

    elements_2find = elements_2find_original
    ele_symbol_charge = []
    ele_symbol_percentage = []

    return elements_2find, ele_symbol_charge, ele_symbol_percentage

def process_result(result, elements_2find):
    tot_a = [0]*len(result)
    for ind, acharges_of_ele in enumerate(result['elements']):
        sub_res = (result[result['elements']==acharges_of_ele])
        total = 0
        for j in sub_res['count']:
            total = total + j
        tot_a[ind] = (total)
        
    result['total_count'] = tot_a
    result['percentage_a'] = round(result['count'] / result['total_count'], 3)

    data = result
    ele_symbol_charge = [[] for _ in range(len(elements_2find))]
    ele_symbol_percentage = [[] for _ in range(len(elements_2find))]
    for index, acharges_of_ele in enumerate(elements_2find):
        charge_state_templist = []
        percentage_templist = []
        for ind, j in enumerate(data['elements']):
            if acharges_of_ele == j:
                charge_state_templist.append(data['charges'].iloc[ind])
                percentage_templist.append(data['percentage_a'].iloc[ind])
        ele_symbol_charge[index] = charge_state_templist
        ele_symbol_percentage[index] = percentage_templist

    return result, ele_symbol_charge, ele_symbol_percentage

def process_element_charges(round1_ele_charg_dict_list):
    r1_ele_c = []
    r1_chag_c = []
    for ind, acharges_of_ele in enumerate(round1_ele_charg_dict_list):
        for key in acharges_of_ele: 
            r1_ele_c.append(key)
            r1_chag_c.append(acharges_of_ele[key])

    if len(round1_ele_charg_dict_list) == 0:
        ele_list = []
        charge_list = []
    else:
        ele_list, charge_list = zip(*sorted(zip(r1_ele_c, r1_chag_c)))

    res = {'elements' : ele_list,
        'charges' : charge_list}

    res = pd.DataFrame(res)
    result = res.groupby(['elements', 'charges']).size().reset_index(name='count')
    tot_a = [0]*len(result)

    return result, tot_a

def save_charge_dict(elements_2find, round1_ele_charg_dict_list, temp_formula):
    filename = 'charge_dict'
    for i in elements_2find:
        filename = filename + '_' + i
    file_name = filename + '.pkl'
    with open(file_name,'wb') as output:
        pickle.dump(round1_ele_charg_dict_list, output)

    print()
    print(len(temp_formula), 'compounds were solved and saved as', file_name)

def get_temp_formula(round1_ele_charg_dict_list, pformula_r1, r1_o1index, round1_oindex_list):
    temp_formula = []
    for ind, i in enumerate(round1_ele_charg_dict_list):
        temp_formula.append(pformula_r1[r1_o1index.index(round1_oindex_list[ind])])
        if ind < 5:
            print(pformula_r1[r1_o1index.index(round1_oindex_list[ind])], i)
    return temp_formula

def process_formula(pformula_r1, r1_int_ele, r1_int_charg, ele_neg, r1_o1index):
    round1_ele_list = []
    round1_charge_list = []
    round1_oindex_list = []
    round1_ele_charg_dict_list = []
    count = 0

    for ind, acharges_of_ele in enumerate(pformula_r1):
        compound = mg.Composition(str(acharges_of_ele))

        elements = []
        na_list = []
        ele_flag = []
        for ele in compound.elements:
            elements.append(str(ele))
            na_list.append(compound.get_el_amt_dict()[str(ele)])
            ele_flag.append(str(ele) in r1_int_ele)
        
        if sum(ele_flag) == len(elements):
            elements, zeropairs = solve_charge_round1(elements, na_list, r1_int_ele, r1_int_charg, ele_neg)
        
            if sum([j==1 for j in [len(i) for i in zeropairs]]) == len(elements):
                for acharges_of_ele in elements:
                    if not acharges_of_ele in round1_ele_list:
                        round1_ele_list.append(acharges_of_ele)
                        round1_charge_list.append(zeropairs[elements.index(acharges_of_ele)][0])
                    elif acharges_of_ele in round1_ele_list:
                        ind1 = round1_ele_list.index(acharges_of_ele)
                        if isinstance(round1_charge_list[ind1], int):
                            c1 = [round1_charge_list[ind1]]   
                        else:
                            c1 = list(round1_charge_list[ind1])

                        if not zeropairs[elements.index(acharges_of_ele)][0] in c1:
                            c1.append(zeropairs[elements.index(acharges_of_ele)][0])
                            round1_charge_list[ind1] = c1
                round1_oindex_list.append(r1_o1index[ind])
                round1_ele_charg_dict = {}
                for indi, elei in enumerate(elements):
                    round1_ele_charg_dict[elei] = zeropairs[indi][0]
                round1_ele_charg_dict_list.append(round1_ele_charg_dict)
                count +=1

    return round1_oindex_list, round1_ele_charg_dict_list

def solve_charge_round1 (Eall, Nall, ele_lvl, charg_lvl, ele_neg):
    zeropairs = [[] for _ in range(len(Nall))]
    
    cp = polish_cpr1 (Eall, ele_lvl, charg_lvl, ele_neg)
    
    flag = 1
    for i in cp:
        if not i:
            flag = 0
            #print('There is at least one zero pool. Skip...')
    
    #print(cp)
    total_cstate = np.prod([len(i) for i in cp])
    #print(total_cstate, 'combinations should be calculated')
    
    #excluding the single elements
    if flag == 1 and (not len(cp) == 1):
        
        initial = cp[-1]
        count = -2
        new_frame = []
        for ind, i in enumerate(initial*len(cp[count])):
            new_frame.append([cp[count][int(ind/len(initial))], i])
        count -= 1
        previous_frame = new_frame
        
        for indc, c in enumerate(range(len(Eall)-2)):
            new_frame1 = []
            for ind, i in enumerate(previous_frame*len(cp[count])):
                new_frame1.append([cp[count][int(ind/len(previous_frame))]]+ i)

            previous_frame = new_frame1
            count -=1
        #print(len(previous_frame),'combinations found.')
    
        #counter check
        if not total_cstate == len(previous_frame):
            print("Check for Errors!")
        else:
            for i in previous_frame:
                total = 0
                for ind, j in enumerate(i):
                    total = total + (j*Nall[ind])
                if total == 0:
                    for ind, j in enumerate(i):
                        zeropairs[ind].append(j)
        #print('Possible solutions :' , zeropairs)
    return Eall, zeropairs


def polish_cpr1 (Eall, ele_lvl, charg_lvl, ele_neg):
    cp = []
    for i in Eall:
        if isinstance(charg_lvl[ele_lvl.index(i)], int):
            cp.append([charg_lvl[ele_lvl.index(i)]])
        else:
            cp.append(charg_lvl[ele_lvl.index(i)])

    #print(elements, 'CP', cp)
    elen = find_electron_negativity_single(Eall, ele_neg)
    
    if not 0 in elen:
        mini = elen.index(min(elen))
        #print(elen, 'min is', elements[mini], 'with',  cp[mini])
        #print("Removing negative charges from", elements[mini])
        temp_charge = cp[mini]
        temp_charge = [i for i in temp_charge if i > 0]
        cp[mini] = temp_charge
        #print(elements[mini], 'in now with',  cp[mini])

        maxi = elen.index(max(elen))
        #print(elen, 'max is', elements[maxi], 'with',  cp[maxi])
        #print("Removing positive charges from", elements[maxi])
        temp_charge = cp[maxi]
        temp_charge = [i for i in temp_charge if i < 0]
        cp[maxi] = temp_charge
        #print(elements[maxi], 'in now with', cp[maxi])
        #print('CP Final : ', Eall, cp)
    #else:
        #print('skip elen discrimination. 0 in elen')
        
    return cp

def find_electron_negativity_single (ele, ele_neg):
    elen = []
    for e in ele:
        ele = 0
        temp_frame = ele_neg.loc[ele_neg['symbol'] == e]
        if not temp_frame.empty:
            ele = float(temp_frame['en_pauling'])
        elen.append(ele)
    return elen

def extract_wiki_data(file_name):
    wiki_table = pd.read_excel(file_name, header=None)
    wiki_ele = []
    wiki_charges = []
    for ind, acharges_of_ele in enumerate(wiki_table[2]):
        wiki_ele.append(acharges_of_ele)
        temp = []
        for j in range(3,18,1):
            if not (pd.isnull(wiki_table[j][ind]) or wiki_table[j][ind]==0):
                if isinstance(wiki_table[j][ind], float):
                    temp.append(int(wiki_table[j][ind]))
                else:
                    temp.append(int(re.sub(u"\u2212", "-", wiki_table[j][ind])))
        wiki_charges.append(temp)
    return wiki_ele, wiki_charges

def find_elements(elements_2find_original, recorded_ele_list):
    elements_2find = []
    for ele in elements_2find_original:
        if not ele in recorded_ele_list:
            elements_2find.append(ele)
    return elements_2find

def filter_dataframe_by_elements(df, elements_2find):
    ind_comp = [False]*len(df['elements'])
    for acharges_of_ele in elements_2find:
        indexes = [str(acharges_of_ele) in s for s in df['elements']]
        ind_comp = [a or b for a, b in zip(indexes, ind_comp)] 
        
    filtered_df = df[ind_comp]
    return filtered_df

def load_charge_probability_database(filename):
    if os.path.isfile(filename):
        cb_data_base = pd.read_csv(filename)
        recorded_ele_list = list(cb_data_base['elements'])
        recorded_charge_states = list(cb_data_base['charge_states'])
        recorded_charge_probabilities = list(cb_data_base['probabilities'])
    else:
        recorded_ele_list = []
        recorded_charge_states = []
        recorded_charge_probabilities = []

    return recorded_ele_list, recorded_charge_states, recorded_charge_probabilities

def filter(elements_2find_original, charge_mixing,num, freq, path):
    # Section 1 (define variables)
    df = pd.read_hdf("all_materials_9June2022.h5")

    #maximum number of atoms in a compound to be generated
    max_atom_in_compound = 20

    #to be used to ignore low frequency charge states
    cut_off_freq_percentage = freq

    # Section 1a (get data from database)
    filename = 'charge_probability_database.csv'    
    recorded_ele_list, recorded_charge_states, recorded_charge_probabilities = load_charge_probability_database(filename)
    print('Recorded elements :', recorded_ele_list)
    print('Ternary system :', elements_2find_original)
    elements_2find = find_elements(elements_2find_original, recorded_ele_list)

    print('Elements to find CS probability :', elements_2find)

    # Section 2 (load data)

    df_all = df #save it for later to find icsd_ids
    df = df[[bool(s) for s in df['icsd_ids']]]
    #df = df[[s==0 for s in df['e_above_hull']]]

    
    wiki_ele, wiki_charges = extract_wiki_data('Wiki_charge_list.xlsx')

    ##adding +3 to 'Sr'
    ind1 = wiki_ele.index('Sr')
    wiki_charges[ind1] = [1,2,3]
        
    '''
    for ind, i in enumerate(wiki_ele):
        print(i, wiki_charges[ind])
    '''

    ele_neg = pd.read_excel('Pauling_ENs.xlsx')
    ele_neg = ele_neg[[not math.isnan(x) for x in ele_neg['en_pauling']]]

    # Section 3 (Filter dataset by elements)   
    filtered_df = filter_dataframe_by_elements(df, elements_2find)
    print(filtered_df)

    # Section 4 (Assign variables)

    pformula_r1 = list(filtered_df['pretty_formula'])
    r1_o1index = list(filtered_df.index)
    #print('Total compounds =' , len(pformula_r1))


    r1_int_ele = list(wiki_ele)
    r1_int_charg = list(wiki_charges)

    round1_oindex_list, round1_ele_charg_dict_list = process_formula(pformula_r1,
                                                                    r1_int_ele,
                                                                    r1_int_charg,
                                                                    ele_neg, 
                                                                    r1_o1index)



    # Section 6a (printing 5 compuonds that had been solved)
    temp_formula = get_temp_formula(round1_ele_charg_dict_list, pformula_r1, r1_o1index, round1_oindex_list)
    print(temp_formula)

    # Section 6b (Save CS probability temporarly)
    
    if False:
        save_charge_dict(elements_2find, round1_ele_charg_dict_list, temp_formula)


    # Section 7 (NO_CHARGE MIXING METHOD - results)
    result, tot_a = process_element_charges(round1_ele_charg_dict_list)

    result, ele_symbol_charge, ele_symbol_percentage = process_result(result, elements_2find)

    
    # Section 7s (Update charge probability Database)
    elements_2find, ele_symbol_charge, ele_symbol_percentage = update_charge_probability_database(elements_2find,
                                                                                                recorded_ele_list,
                                                                                                ele_symbol_charge,
                                                                                                ele_symbol_percentage, 
                                                                                                elements_2find_original)

    for acharges_of_ele in elements_2find:
        ele_symbol_charge.append(json.loads(recorded_charge_states[recorded_ele_list.index(acharges_of_ele)]))
        ele_symbol_percentage.append(json.loads(recorded_charge_probabilities[recorded_ele_list.index(acharges_of_ele)]))

    for ind, acharges_of_ele in enumerate(elements_2find):
        print(acharges_of_ele, ele_symbol_charge[ind], ele_symbol_percentage[ind])

    
    # Section 7t (NO_CHARGE MIXING METHOD - analyze result details)

    if False:
        temp_formula = filter_formulas_by_element_charge('Ge', 2, round1_ele_charg_dict_list,
                                                        pformula_r1, r1_o1index, round1_oindex_list)



    # Section 8 (Charge preparation - elenge rule - ignore very few occurances)


    cp, cp_perc = filter_charges_by_cutoff(elements_2find, ele_symbol_charge, 
                                           ele_symbol_percentage, cut_off_freq_percentage)

    ### electronegativity rule out
    elen = find_electron_negativity_single(elements_2find, ele_neg)

    cp, cp_perc = filter_charge_polarity(elen, cp, cp_perc)
    
    #for ind, acharges_of_ele in enumerate(elements_2find):
    #    print(acharges_of_ele, cp[ind], cp_perc[ind])

    # Section 9 (Generate valid compounds - charge mixing)

    #charge_mixing = True

    range_max = max_atom_in_compound - len(cp) + 1
    number_of_atoms_per_ele = list(range(1,range_max+1))

    #valid comb are the list of atoms with less than the max
    all_possible_comb = list(pdt(number_of_atoms_per_ele, repeat = len(cp)))
    comb2remove = [(sum(i) <= max_atom_in_compound) for i in all_possible_comb]
    valid_comb = [i for (i,remove) in zip(all_possible_comb, comb2remove) if remove]

    #print('list of availabe charges', cp)
    #print('first 10 valid combo pairs', valid_comb[:10])
    #print('valid number of atoms per elements pairs = ', len(valid_comb))


    # Section 9a (filter charge neutral compounds)
    filtered_charge_comb, filtered_atom_counts = filter_charge_combinations(valid_comb, cp, charge_mixing)

    #print('first five charge_combo', filtered_charge_comb[:5])
    #print('first five atom counts', filtered_atom_counts[:5])
    #print('check sum : valid combo & atom counts :', len(filtered_atom_counts),'&', len(filtered_charge_comb))



    # Section 10 (Assign probability to the charge calculated)
    one_number_probability, perc_calculated, raw_perc_for_each_charg, charges_splits, ischargemixed, ele_charge_counts = calculate_charge_probabilities(filtered_atom_counts, filtered_charge_comb, cp_perc, cp)

    #print(filtered_atom_counts)
    #print(filtered_charge_comb)
    #print(raw_perc_for_each_charg)
    #print(perc_calculated)
    #print(len(one_number_probability), 'percentages calculated')

    # Section 10e (Check how many unique formula were generated)
    


    unique_formula = unique(filtered_atom_counts)
    #print('number of unique formulas = ', len(unique_formula))

    #for ind, i in enumerate(filtered_charge_comb[:5]):
        #print(i, charges_splits[ind], ele_charge_counts[ind], ischargemixed[ind])

    # Section 11 (prepare Results)
    d2s = create_dataframe(elements_2find, filtered_atom_counts, filtered_charge_comb,
                        charges_splits,
                        ele_charge_counts, 
                        ischargemixed,
                        raw_perc_for_each_charg,
                        perc_calculated,
                        one_number_probability)


    #d2s.head(2)
    # formula_list = list(d2s["composition"])
    # reduced_formula = [mg.Composition(i).reduced_formula for i in formula_list]
    # d2s["formula"] = formula_list
    # d2s["reduced_formula"] = reduced_formula
    # #remove duplicates
    # d2s = d2s.drop_duplicates(subset=['reduced_formula'])
    # save_to_csv(d2s, elements_2find)


    # Section 15 (Definitions to find icsd_ids)

    formula_list, df_filtered = filter_dataframe(df_all, elements_2find, d2s)

    # Section 16 (find icsd_id and others paramters)

    icsd_ids, mp_ids, mp_FEPA, mp_eabovehull, pretty_formula = match_formulas_with_database(formula_list, df_filtered)


    #Secion 17 (put everything together into dataframe format)
    attached_data2save = create_final_dataframe(elements_2find, cp, charges_splits,
                                                formula_list, 
                                                pretty_formula, 
                                                mp_ids, icsd_ids, 
                                                mp_FEPA, 
                                                mp_eabovehull, 
                                                d2s)
    #attached_data2save.head(2)

    attached_data2save = attached_data2save.drop_duplicates(subset=['reduced_formula']).reset_index(drop=True)
    
    #Section 18 (Save Data)
    save_to_csv(attached_data2save, elements_2find, path)
    return attached_data2save

def main():
    elems = ['Cs','Sn','I']
    charge_mixing = False
    filter(elems, charge_mixing)


if __name__ == "__main__":
    main()