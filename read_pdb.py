import pandas as pd
import os
import ampal


################################################## tool function
def pdb2df(protPdb):
    columns = ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']
    data = []
    with open(protPdb, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_type = line[0:6].strip()
                try:
                    atom_id = int(line[6:11].strip())
                except:
                    atom_id = str(line[6:11].strip())
 
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                if chain_id == '':
                    chain_id = "A"
                try:
                    res_id = int(line[22:26].strip())
                except:
                    res_id = str(line[22:26].strip())
 
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                try:
                    z = float(line[46:54].strip())
                except:
                    z = str(line[46:54].strip())
                try:
                    occupancy = float(line[54:60].strip())
                except:
                    occupancy = str(line[54:60].strip())
          
                try:
                    temp_factor = float(line[60:66].strip())
                except:
                    temp_factor = " "*6
                element = line[76:78].strip()
 
                data.append([atom_type, atom_id, atom_name, res_name, chain_id, res_id, x, y, z, occupancy, temp_factor, element])
 
    return pd.DataFrame(data, columns=columns)

'''using ampal we determine if the ligand's atom could interact other atom from different chains, edit by handing 23/04/24'''

def deter_if_interact(input_pathway):
#for index,row in iter
    df = pdb2df(input_pathway)
    filtered_df_ligand = df[(df['ATOM'] == 'HETATM') & (df['RES_NAME']  != 'HOH')]
   # print(filtered_df_ligand)
    try:
        myprotein = ampal.load_pdb(input_pathway)
       
    except ValueError as e:
        print(f"An error occurred with file {input_pathway}: {e}")
        return None  # return None or an appropriate value to indicate that the file was skipped

    result = {}
    for index,row in filtered_df_ligand.iterrows():
        coordinate_x = row['X']
        coordinate_y = row['Y']
        coordinate_z = row['Z']
        parent_chain = row['CHAIN_ID']
        ligand_name = row['RES_NAME']   
    
        each_atom_around = myprotein[0].is_within(3,(coordinate_x,coordinate_y, coordinate_z))
        
    
        interact_atom_onchain =[]
        interact_atom_onchain_id = []
        interact_atom_onchain_AA = []
        for i in each_atom_around:
        #    print(i.parent.parent.id)
            if parent_chain != i.parent.parent.id:
                interact_atom_onchain.append(i)
                interact_atom_onchain_id.append(i.id)
                interact_atom_onchain_AA.append(i.parent)
        
        if interact_atom_onchain and len(interact_atom_onchain) >=2 :
            result[ligand_name] = interact_atom_onchain_id,interact_atom_onchain_AA,interact_atom_onchain

    return result

def find_res_name_with_multiple_chains(df):
    res_list =[]

    chain_counts = df.groupby('RES_NAME')['CHAIN_ID'].nunique() 
    ####### 1.1 The dataframe is group by the res_name, and for each group we counts unique chains_ids,return a series

    # 1.2 Find 'RES_NAME's with more than 2 unique 'CHAIN_ID's
    res_names_with_multiple_chains = chain_counts[chain_counts >= 2].index.tolist()
    res_list.append(res_names_with_multiple_chains)

    return res_list

def read_pdb_files_in_batch(directory):

    structures = []
    pdb_id = []
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):
            pdb_id.append(filename) 
            pdb_path = os.path.join(directory, filename)
            structures.append(pdb_path)        
    return structures,pdb_id

######################################################################### write and call above
def wirte_to_file(folder_in,folder_out):
    for input,output in zip(folder_in,folder_out):

        pdb_batch_pathway = input
        output_file = output
        dir_path_list = []

        num = 0
        # step1:Loop over all directories
        for dir_name in os.listdir(pdb_batch_pathway):
           
            dir_path = os.path.join(pdb_batch_pathway, dir_name) 
            dir_path_list.append(dir_path) # Loop over all directories, get folder_pathway list
            num += 1
            batch_list = dir_path_list
        #    print(batch_list)
           
            # step2: Process the file one by one
            for dir_batch in batch_list:
                pdb_batch_list,pdb_id = read_pdb_files_in_batch(dir_batch)
             #    print(pdb_id)
                
                pdb_id_dict = {'pro_id':pdb_id} # make a dict to write
                df_pdb_id = pd.DataFrame(pdb_id_dict)
                            
            ligand_list =[]
            ligand_interact_atom_list = []

            for input_pathway in pdb_batch_list:

                df = pdb2df(input_pathway)
                
                resiDf = df[(df["ATOM"] == "HETATM") & (df["RES_NAME"] != "HOH")]
                ligand_dict_iter = find_res_name_with_multiple_chains(resiDf)[0]
                
                if ligand_dict_iter != []:
                    ligand_interact_atom = deter_if_interact(input_pathway) ## already a dict
                    ligand_interact_atom_list.append(ligand_interact_atom)
                else:
                    ligand_interact_atom_list.append(None)

                ligand_list.append(ligand_dict_iter)

                ligand_dict = {'ligand_id':ligand_list } 
                df_ligand_id = pd.DataFrame(ligand_dict)

                ##### thrid column
                ligand_interact_atom_dict = {'ligand_interact_atom':ligand_interact_atom_list}
                df_ligand_interact_atom = pd.DataFrame(ligand_interact_atom_dict)

                       
            # Write result_ligand to CSV in sequence, if the file does not exist, mode='a' will create it.
            result_ligand = pd.concat([df_pdb_id,df_ligand_id,df_ligand_interact_atom],axis=1)
            result_ligand.to_csv(output_file, mode = 'a',index= False)

##########################################################

base_input = '/mnt/e/'
base_oupt = '/home/hdwang/Ampal_learning/3.Find_ligand_@homopolymer/test_class'

folder_in = [] # get a whole list of folder pathway
folder_out = []

for i in range(1,2):
    
    # Construct the folder name
    folder_name_in = 'unzip_group' + str(i)
    folder_name_out = 'ligand_info_group' + str(i)
    
    # Construct the full path
    full_input_path = os.path.join(base_input, folder_name_in)
    full_output_path = os.path.join(base_oupt, folder_name_out)
    if not os.path.exists(full_output_path):
        df = pd.DataFrame()
        df.to_csv(full_output_path, index=False) 
    # Get the list of files in the folder
    folder_in = full_input_path
    folder_out = full_output_path
   # print(folder_in,folder_out)
 
    wirte_to_file([folder_in],[folder_out])

################################################################








