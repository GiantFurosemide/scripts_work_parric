"""
if you want a .py command, use pymol API

like
> cmd.save(filename[, selection[, state[, format]]])

more information can be found at:
https://pymolwiki.org/index.php/Category:Commands

this script creates .pml file with cmd style command


目的：做成一个可以 在指定数值范围（contact score），连续显示结构上contact 的分布的 pdb.file 就好像 movie一样

面对：trajectory 的file也可以做出同样的东西

distance trace
"""


########################################################################################################################
########################################################################################################################
########################################################################################################################
#
#   BEFORE USING, PLEASE CHANGE CONFIGURATION IN CLASS 'ConfigData'
#
########################################################################################################################
########################################################################################################################
########################################################################################################################
# 1, configure file
class ConfigData:
    PROTEIN_PDB_CODE = '6ibb'  # should always be lower case
    OUT_PATH = './'
    PYMOL_OUT_PATH = f'{OUT_PATH}/{PROTEIN_PDB_CODE}_show_contact.pml'
    # change this
    JSON_FILE_PATH = '/Users/muwang/Documents/scripts_work_parric/RRCS_based/pdb/6ibb_clean.cif.pdb.cscore.json'

    PDB_OUT_PATH = f'{PROTEIN_PDB_CODE}_contact_marked.pse'  # this what we finally want,pse will save color change.


########################################################################################################################
########################################################################################################################
########################################################################################################################
# 2, make command
# include save and open


def load_json_dict(fp: str) -> dict:
    import json
    with open(str(fp).strip()) as inf:
        temp_dict = json.load(inf)
    return temp_dict


def brief_report_contact_score(contact_dict: dict) -> float:
    import pandas as pd

    contact_df = pd.DataFrame.from_dict(contact_dict)
    contact_df = contact_df.fillna(0)

    length = contact_df.max().max()
    print('max of score is: ', length)
    return float(length)


def contact_dict_parser(contact_dict: dict, limit_range: tuple = (0, 99999)) -> list:
    """
    # to parse the resi and resj name:str to 3 part:
    # 1. chain ID
    # 2. res number
    # 3. res name
    # here we only use 1 and 2 for pml command
    # only the resi-j's score in limit_range will be appended to list.
    :param limit_range: res pair which score in this range will be recorded
    :param contact_dict: dict stored [resi][resj][score]
    :return: [(chain ID,res number), ....]


    """
    res_list = []  # list to store output
    for resi in contact_dict:
        resi_flag = 0
        resj_list = contact_dict[resi].keys()
        for resj in resj_list:
            score = contact_dict[resi][resj]
            if limit_range[0] <= score < limit_range[1]:
                resi_chain_id = str(resi).split(':')[0].strip()
                resi_number = str(resi).split(':')[1].split('_')[0].strip()
                resj_chain_id = str(resj).split(':')[0].strip()
                resj_number = str(resj).split(':')[1].split('_')[0].strip()
                res_list.append((resj_chain_id, resj_number))
                if resi_flag == 0:
                    resi_flag = 1
                    res_list.append((resi_chain_id, resi_number))
    return res_list


def format_cmd2pml(res_list: list,pdb_code:str) -> list:
    """
    make pml cmd sentence from res list and store in a list.
    :param res_list: [(chain ID,res number), ....]
    :return: list of command
    """
    cmd_list = []
    for chain_id, res_number in res_list:
        chain_id = str(chain_id).strip()
        res_number = str(res_number).strip()
        cmd_list.append(f'color red , /{pdb_code.lower()}/{chain_id}/{chain_id}/{res_number}\n')
    return cmd_list


def make_pml_cmd(fp: "json file path",pdb_code,limit_range=(0, 9999)) -> list:
    contact_dict = load_json_dict(fp)
    res_list = contact_dict_parser(contact_dict, limit_range)
    cmd_list = format_cmd2pml(res_list,pdb_code)
    return cmd_list


# 3, make pml
# read from json


json_file_path = '/Users/muwang/Documents/scripts_work_parric/RRCS_based/pdb/6ibb_clean.cif.pdb.cscore.json'


def make_pml_file(json_path: str, pml_out_path: str, pdb_code, pdb_out, limit_range=(0, 9999)):
    header = [
        f'fetch {pdb_code}\n',
        'as cartoon\n',
        'show sticks,all\n',
        'color green ,all\n'
    ]
    ending = [
        f'save {pdb_out}'
    ]

    command_load = header + make_pml_cmd(json_path, pdb_code,limit_range) + ending
    with open(pml_out_path, 'w+') as fi:
        fi.writelines(command_load)


# 4, run pml
# Pymol -c command.pml
# Pymol -d command
# we can do this a cmd-line style
def run_pymol(json_path, pymol_out_path, pdb_code, pdb_out_path, limit_range: tuple = (0, 9999)):
    make_pml_file(json_path, pymol_out_path, pdb_code, pdb_out_path, limit_range)
    import os
    os.system(f'pymol -cQ {pymol_out_path}')
    os.system(f'pymol {pdb_out_path}')


if __name__ == '__main__':
    print("""
    BEFORE USING, PLEASE CHANGE CONFIGURATION IN CLASS 'ConfigData'.
    """)

    cf = ConfigData()
    cf.OUT_PATH = './'
    cf.PROTEIN_PDB_CODE = '6rnk'  # should always be lower case
    cf.SCORE_LIMIT = (7, 10)
    cf.PYMOL_OUT_PATH = f'{cf.OUT_PATH}/{cf.PROTEIN_PDB_CODE}_{cf.SCORE_LIMIT[0]}_{cf.SCORE_LIMIT[1]}_show_contact.pml'
    # json file can be calculated by RRCS_change.py
    cf.JSON_FILE_PATH = '/Users/muwang/Documents/scripts_work_parric/RRCS_based/pdb/6rnk_clean.cif.pdb.cscore.json'
    # this what we finally want,pse will save color change.
    cf.PDB_OUT_PATH = f'{cf.PROTEIN_PDB_CODE}_{cf.SCORE_LIMIT[0]}_{cf.SCORE_LIMIT[1]}_contact_marked.pse'

    # get max score from contact map as reference for cf.SCORE_LIMIT
    max_score = brief_report_contact_score(load_json_dict(cf.JSON_FILE_PATH))
    max_score_int = int(max_score) + 1

    import numpy as np

    limit_range = np.arange(0, max_score_int, 0.5)
    limit_range = np.array_split(limit_range, 5)

    #run_pymol(cf.JSON_FILE_PATH, cf.PYMOL_OUT_PATH, cf.PROTEIN_PDB_CODE, cf.PDB_OUT_PATH, limit_range=cf.SCORE_LIMIT)
    for i in range(len(limit_range)):
        min_s = 0
        max_s = 0
        if i < len(limit_range)-1:
            min_s = limit_range[i][0]
            max_s = limit_range[i + 1][0]
        elif i == len(limit_range) - 1:
            min_s = limit_range[i][0]
            max_s = limit_range[i][-1]

        cf.SCORE_LIMIT = (min_s, max_s)
        cf.PYMOL_OUT_PATH = f'{cf.OUT_PATH}/{cf.PROTEIN_PDB_CODE}_{cf.SCORE_LIMIT[0]}_{cf.SCORE_LIMIT[1]}_show_contact.pml'
        cf.PDB_OUT_PATH = f'{cf.PROTEIN_PDB_CODE}_{cf.SCORE_LIMIT[0]}_{cf.SCORE_LIMIT[1]}_contact_marked.pse'
        run_pymol(cf.JSON_FILE_PATH, cf.PYMOL_OUT_PATH, cf.PROTEIN_PDB_CODE, cf.PDB_OUT_PATH,
                  limit_range=cf.SCORE_LIMIT)

        if i == 3:
            min_s = limit_range[i+1][0]
            max_s = limit_range[i+1][-1]

