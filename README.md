# Filters
Code developed by Basita Das to screen for materials. This code is based on the origanl code developed by Muang Thway, and can be foudn at https://github.com/maungthway/Oxidation-State-Probability-Filter.git 

The dataset used to the analysis can be download rom here (https://www.dropbox.com/scl/fi/yc7pgthfllij9ne4icw1d/all_materials_9June2022.h5?rlkey=72wl4m6a61dtvsvtqqxpfykm5&dl=0). This data was downloaded from the Materials Project database on June 2022.

## Files: 
1.	**run.py:** main file from where the code is run
2.	**charge_neutrality_and_oxidation_state_filter.py:** file where the charge neutrality filter, electronegativity balance filter, and the two-oxidation state filter is implemented.\
  **File dependencies:**\
  a)	Pualing_ENs.xlsx\
  b)	Wiki_charge_list.xlsx\
  c)	charge_probabiloty_database.csv
3.	**stoichiometry_filters.py:** file where the intra phase diagram filter, and cross phase diagram filter is implemented.
4.	**Data_analysis.py:** file used to generate the ternary plots shown in Figure 3 in the paper.

## Running instructions : 
1. Create python environment according to **requirements.txt**
2.  The code is run from the **run.py** file.
3. Specify the results directory by updating the 'dir_name'
4. Hit run, the results will be generated in the results directory in the sub_directory specified by the user.

