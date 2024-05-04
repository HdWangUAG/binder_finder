import requests
import re
import pandas  as pd
import ast
import os

###################
def classify_ligand(name,formula):

    # Define a list of amino acids
    amino_acids = ['glycine', 'alanine', 'valine', 'leucine', 'isoleucine', 'proline', 'phenylalanine', 'tryptophan', 'methionine', 'cysteine', 'serine', 'threonine', 'asparagine', 'glutamine', 'tyrosine', 'histidine', 'lysine', 'arginine', 'aspartic acid', 'glutamic acid']
    for amino_acid in amino_acids:
        if amino_acid in name.lower():
            return 'Amino Acid'
    
    if 'ion' in name.lower():
        return 'Ion'
    if 'phosphate' in name.lower():
        return 'Phosphate'
    if 'acid' in name.lower():
        return 'Acid'

    ########### Extract the counts of C, H, and O from the formula
    c_count = int(re.search('C(\\d+)', formula).group(1)) if re.search('C(\\d+)', formula) else 0
    h_count = int(re.search('H(\\d+)', formula).group(1)) if re.search('H(\\d+)', formula) else 0
    o_count = int(re.search('O(\\d+)', formula).group(1)) if re.search('O(\\d+)', formula) else 0
    ########### classify the formula
    if h_count == 2* c_count and  o_count ==2:
        return'Lipid'
    elif h_count == 2*c_count and o_count == c_count:
        return 'Sugar'
    elif h_count == 2*c_count + 2 and o_count == c_count:
        return 'Alcohols'
    else:
        return 'Complex'
    

#####################
def get_chem_comp_details(ccd_id):
    base_url = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/"
    response = requests.get(base_url + ccd_id)
    data = response.json()
    return data


#####################
def get_ligand_chem_name(ccd_id):
    try:
        if ccd_id =='8NK':
            ccd_id = 'M7G'

        chem_comp_details = get_chem_comp_details(ccd_id)
        name = chem_comp_details[ccd_id][0]['name']
        formula = chem_comp_details[ccd_id][0]['formula']
        ligand_type = classify_ligand(name, formula)

    except KeyError:
        ligand_type = []
    return ligand_type


# Define the function to get chem_name
def get_chem_name(row):
    list_elements = ast.literal_eval(row)
    row_chemname_list = []
    for element in list_elements:
        ccd_id = element
        ligand_type = get_ligand_chem_name(ccd_id)
        row_chemname_list.append(ligand_type)
        print(row_chemname_list)

    return row_chemname_list


#############################
base_input = '/home/hdwang/Ampal_learning/3.Find_ligand_@homopolymer/filtered_ligand_info/'
for i in range(1,55):
    file = 'filter_ligand_group' + str(i)
    input_pathway = os.path.join(base_input,file)

    df = pd.read_csv(input_pathway)

    #df = df.drop(columns=[ 'Unnamed: 0'])
    df = df[df['ligand_id'] != "['NAG']"] 
    df = df[df['ligand_id'] != "['GOL']"] 
    df.to_csv(input_pathway,index=False)




''' concact all files into one '''
import pandas as pd
import glob


def concatenate_csv_files(path_to_files, output_file):
    # Get a list of all CSV files in the target directory
    csv_files = glob.glob(path_to_files + '/*')

    # Create a list to hold dataframes
    dfs = []

    # Loop through the list of file names and read each file into a dataframe, then append it to the list
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        dfs.append(df)

    # Concatenate all dataframes in the list
    combined_df = pd.concat(dfs, ignore_index=True)

    # Write the concatenated dataframe to a new CSV file
    combined_df.to_csv(output_file, index=False)
path_to_files = '/home/hdwang/Ampal_learning/3.Find_ligand_@homopolymer/filtered_ligand_info'

csv_files = glob.glob(path_to_files+ '/*')

concatenate_csv_files('/home/hdwang/Ampal_learning/3.Find_ligand_@homopolymer/filtered_ligand_info','combined.csv')






