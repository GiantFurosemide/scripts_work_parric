# -*- coding: utf-8 -*-

'''
process pdb file to calculate contact score and write to csv.
    main function is prosess_save2csv()

usage: python3 RRCS_change.py file.pdb
input .pdb
output .csv
'''
import sys


def calc_contact(pdbbase: "pdb file") -> (dict, dict):
    fi = open(pdbbase, 'r')
    all_lines = fi.readlines()
    fi.close()
    atom_lines = [l for l in all_lines if l[0:6] == 'ATOM  ']
    dict_coord = {}  # dict to store coordinates. dict_coord[res][atom] = (x,y,z,occ)
    atomnum_2_name = {}  # map atom number to atom name, in order to find N, CA, C, O
    contact_score = {}  # dict to store final results. contact_score[ires][jres] = contact_score.
    for line in atom_lines:
        # retrive info from each atom line
        atom_num = int(line[6:11].strip())
        atom_name = line[12:16].replace(' ', '_')
        res_name = line[17:20]
        res_num = int(line[22:26].strip())
        chain_id = line[21:22]
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occ = float(line[54:60].strip())
        res = chain_id + ':' + str(res_num) + '_' + res_name  # res
        atomnum_2_name[atom_num] = atom_name
        if res_num <= 0:
            continue
        if res not in dict_coord:
            dict_coord[res] = {}
        dict_coord[res][atom_num] = (x, y, z, occ)  # store res:atom:(atom coord x,y,z,occ) in dict_coord

    ####
    # above we read coord from .pdb file and store to a dict in hierarchy res::atom::(x,y,z,occ)
    ####

    # list shortening
    num_of_aa = len(dict_coord)
    res_list = dict_coord.keys()
    ires_contact = {}  ## to save interact res' contact score
    for ires in res_list:  # ires is res:str , loop over every res in protein loop1
        ires_contact[ires] = []  # ires_contact:dict {res:str=[]:list}
        ires_num = int(ires.split(':')[1].split('_')[0].strip())  # got res number from ires
        for jres in res_list:  # jres is res:str. for res in loop1, loop over every res in protein.
            jres_num = int(jres.split(':')[1].split('_')[0].strip())  # got res number from jres
            if jres_num <= ires_num:  # move to next step until jres_number > ires_number to avoid repeat
                continue
            jres_flag = 0
            atom_in_ires = dict_coord[ires].keys()
            atom_in_jres = dict_coord[jres].keys()

            # 之后对每个res-res loop over atom
            for iatom in atom_in_ires:
                (ix, iy, iz, iocc) = dict_coord[ires][iatom]
                for jatom in atom_in_jres:
                    (jx, jy, jz, jocc) = dict_coord[jres][jatom]
                    dx = abs(ix - jx)
                    dy = abs(iy - jy)
                    dz = abs(iz - jz)
                    if dx < 4.63 and dy < 4.63 and dz < 4.63:
                        jres_flag = 1
                        ires_contact[ires].append(jres)
                        break
                    # 4.63 来源未知
                if jres_flag:
                    break

                # 如果对于一个res-res对中的至少一个atom-atom对 的 xyz 都小于 4.63, 那么将在 dir ires_contact 存储 这一对 res
        # 以上保存了对于 res i 的所有距离小于4.63的res-res对 于ires_contact

        # loop over the shortened list
        contact_score[ires] = {}
        for kres in ires_contact[ires]:  # 取出刚才找到的相近的jres 重命名为kres
            atom_in_kres = dict_coord[kres].keys()  # return a list of natoms
            kres_num = int(kres.split(':')[1].split('_')[0].strip())
            contact_score[ires][kres] = 0
            total_score = 0
            if abs(ires_num - kres_num) < 5:
                for iatom in atom_in_ires:
                    iatom_name = atomnum_2_name[iatom]
                    if iatom_name in ['_N__', '_CA_', '_C__', '_O__']:
                        continue
                    # 距离小于4 跳过主链的atom
                    (ix, iy, iz, iocc) = dict_coord[ires][iatom]
                    for katom in atom_in_kres:
                        katom_name = atomnum_2_name[katom]
                        if katom_name in ['_N__', '_CA_', '_C__', '_O__']:
                            continue
                        (kx, ky, kz, kocc) = dict_coord[kres][katom]
                        # 距离小于4 跳过主链的atom，只比较侧链
                        d2 = (ix - kx) ** 2 + (iy - ky) ** 2 + (iz - kz) ** 2
                        if d2 >= 21.4369:  # 4.63*4.63 = 21.4369
                            score = 0
                        elif d2 <= 10.4329:  # 3.23*3.23 = 10.4329
                            score = 1.0 * iocc * kocc
                        else:
                            score = (1 - (d2 ** 0.5 - 3.23) / 1.4) * iocc * kocc
                        total_score = total_score + score

                        # 侧链穷尽atom对，距离在3.23:3.23~4.63:4.63打分，计入总分

            elif abs(ires_num - kres_num) > 4:
                for iatom in atom_in_ires:
                    (ix, iy, iz, iocc) = dict_coord[ires][iatom]
                    for katom in atom_in_kres:
                        (kx, ky, kz, kocc) = dict_coord[kres][katom]
                        d2 = (ix - kx) ** 2 + (iy - ky) ** 2 + (iz - kz) ** 2
                        if d2 >= 21.4369:  # 4.63*4.63 = 21.4369
                            score = 0
                        elif d2 <= 10.4329:  # 3.23*3.23 = 10.4329
                            score = 1.0 * iocc * kocc
                        else:
                            score = (1 - (d2 ** 0.5 - 3.23) / 1.4) * iocc * kocc
                        total_score = total_score + score
            contact_score[ires][kres] = total_score
    return contact_score, dict_coord


# 当所有atom距离都大于 4.63 是不被记录的，应该计为0

# a function to eliminate 0 score from calc_contact()
def calc_contact_no_zero(pdbbase: "pdb file") -> dict:
    contact, coord_dict = calc_contact(pdbbase)
    contact_no_zero = {}
    for resi in contact:
        j_res_list = contact[resi].keys()
        for resj in j_res_list:
            score = contact[resi][resj]
            if score > 0:
                if resi not in contact_no_zero:
                    contact_no_zero[resi] = {}
                contact_no_zero[resi][resj] = score
    return contact_no_zero, coord_dict


# pack contact score into csv file
def calc_contact_save2csv(pdbbase: 'pdb file'):
    """
    process pdb file to calculate contact score and write to csv
    """
    contact, _ = calc_contact(pdbbase)
    outf = pdbbase + '.cscore.csv'
    fout = open(outf, 'w')
    for a_res in contact:
        b_res_list = contact[a_res].keys()
        for b_res in b_res_list:
            score = contact[a_res][b_res]
            if score > 0:
                # fout.write('%-12s\t%-12s%10.6f\n' %(a_res, b_res, score))
                fout.write('%-12s,%s,%10.6f\n' % (a_res, b_res, score))
    fout.close()


def calc_contact_save2json(pdbbase: 'pdb file'):
    """
    process pdb file to calculate contact score and write to csv
    """
    contact, _ = calc_contact_no_zero(pdbbase)
    outf = pdbbase + '.cscore.json'
    dict2json(contact, outf)


def delta_rrcs(pdbbase1: 'pdb file', pdbbase2: 'pdb file'):
    """
    delta RRCS

    calculate delta RRCS between two states of receptor, i.e. active inactive.

    only resi exist in the intersection set of pdb_base1 and pdb_base2 have been seen as calculable.

    only calculable resi exist in the union set of pdb_base1 and pdb_base2 contact score,

    have been calculated delta RRCS

    output: (delta_RRCS:dict,(set_not_calculated_a,set_not_calculated_b))

    """
    contact_a, dict_coord_a = calc_contact(pdbbase1)
    contact_b, dict_coord_b = calc_contact(pdbbase2)

    # make intersection of ires for checking
    list_ires_all_a = dict_coord_a.keys()
    list_ires_all_b = dict_coord_b.keys()
    list_ires_contact_a = contact_a.keys()
    list_ires_contact_b = contact_b.keys()

    # this is a set in which elements exist both in protein A and B
    set_ires_inter = set(list_ires_all_a) & set(list_ires_all_b)

    # contact res intersection set
    # set_ires_contact_inter = set(list_ires_contact_a)&set(list_ires_contact_b)
    # contact res union set
    # set_ires_contact_union = set(list_ires_contact_a)|set(list_ires_contact_b)

    # NOT CALCULATED RES FORM A AND B. for report.
    set_not_calculated_a = set(list_ires_all_a) - set_ires_inter
    set_not_calculated_b = set(list_ires_all_b) - set_ires_inter

    # delta_RRCS
    # list [((resi,resj),score),...]
    contact_list_a = []
    contact_list_b = []
    ## shorten list
    for ires in list_ires_contact_a:
        for jres, score in contact_a[ires].items():
            if score > 0:
                contact_list_a.append(((ires, jres), score))

    for ires in list_ires_contact_b:
        for jres, score in contact_b[ires].items():
            if score > 0:
                contact_list_b.append(((ires, jres), score))
    ## shorten list end
    # contact list done
    # print(len(contact_list_a))
    # print(len(contact_list_b))

    delta_RRCS = []
    # filter (resi,resj) in A,B or only A
    for resij_score_a in contact_list_a:
        equal_flag = 0
        for resij_score_b in contact_list_b:
            if resij_score_a[0] == resij_score_b[0]:  # (resi,resj) true means resii in both protein A and B contact
                delta_RRCS.append((resij_score_a[0], resij_score_b[1] - resij_score_a[1]))
                equal_flag = 1
                break
        if equal_flag == 0:
            delta_RRCS.append((resij_score_a[0], 0 - resij_score_a[1]))
    # filter (resi,resj) in only B
    for resij_score_b in contact_list_b:
        equal_flag = 0
        for resij_score_a in contact_list_a:
            if resij_score_a[0] == resij_score_b[0]:  # (resi,resj)
                equal_flag = 1
                break
        if equal_flag == 0:
            delta_RRCS.append((resij_score_b[0], resij_score_b[1] - 0))
    # make delta list
    delta_rrcs_dict = {}
    for (resi, resj), score in delta_RRCS:
        if resi in (set_not_calculated_a | set_not_calculated_b):  # make sure resi is 'calculable'
            continue
        if resi not in delta_rrcs_dict:
            delta_rrcs_dict[resi] = {}
        delta_rrcs_dict[resi][resj] = score
    # delta_RRCS.sort(key=lambda elem:elem[0][0].split(':')[1].split('_')[0])
    # delta_RRCS end

    return delta_rrcs_dict, (set_not_calculated_a, set_not_calculated_b)


# # calculate delta RRCS save to delta_RRCS
# delta_RRCS={}
# for ires_a in list_ires_all_a: # loop over every res in protein A
#     if ires_a in set_ires_inter: # if ture means this res is calculable, exist both in A and B
#         delta_RRCS[ires_a]={}
#         if ires_a in set_ires_contact_union: # ture means this res have at least one RRCS is not 0 :
#             if ires_a in set_ires_contact_inter: # means res probably a "repacking"
#                 list_a_jres_contact = contact_a[ires_a].keys()
#                 list_b_jres_contact = contact_b[ires_a].keys()
#                 for jres_a in list_a_jres_contact:
#                     if jres_a in list_b_jres_contact:# if true means same res pair
#                         delta_RRCS[ires_a][jres_a] = contact_b[ires_a][jres_a]-contact_a[ires_a][jres_a]
#                     elif jres_a not in list_b_jres_contact:# jres's RRCS is 0 in b
#                         delta_RRCS[ires_a][jres_a] = 0-contact_a[ires_a][jres_a]
#                 for jres_b in list_b_jres_contact:
#                     if jres_b not in list_a_jres_contact:# jres's RRCS is 0 in a
#                         delta_RRCS[ires_a][jres_b]=contact_b[ires_a][jres_b]-0
#
#             elif ires_a not in set_ires_contact_inter:
#                 if ires_a in set(list_ires_contact_a): # ires_a's RRCS is 0 in b
#                     list_a_jres_contact = contact_a[ires_a].keys()
#                     for jres_a in list_a_jres_contact:
#                         delta_RRCS[ires_a][jres_a] = 0 - contact_a[ires_a][jres_a]
#                 elif ires_a in set(list_ires_contact_b): # ires_a's RRCS is 0 in a
#                     list_b_jres_contact = contact_b[ires_a].keys()
#                     for jres_b in list_b_jres_contact:
#                         delta_RRCS[ires_a][jres_b]=contact_b[ires_a][jres_b]-0
#         elif ires_a not in set_ires_contact_union:# ture means this A portein res's RRCS is 0
#             continue
#     else:
#         continue
# return delta_RRCS,(set_not_calculated_a,set_not_calculated_b)


def delta_rrcs_save2csv(pdbbase1: 'pdb file', pdbbase2: 'pdb file', outf_path='./delta.csv'):
    """
    input pdbbase1,pdbbase2
    output csv
    """

    deltaRRCS_dict, not_calculated_res_list = delta_rrcs(pdbbase1, pdbbase2)
    outf_path = str(outf_path)

    with open(outf_path, 'w') as outf:
        for ires in deltaRRCS_dict.keys():
            for jres in deltaRRCS_dict[ires].keys():
                outf.write('%-12s,%s,%10.6f\n' % (ires, jres, deltaRRCS_dict[ires][jres]))
    print('Brief report:')
    print('>>>>>>>>>>>>')
    print(not_calculated_res_list[0])
    print('>>>>>>>>>>>>')
    print(not_calculated_res_list[1])
    print('>>>>>>>>>>>>')


def delta_rrcs_save2json(pdbbase1: 'pdb file', pdbbase2: 'pdb file', outf_path='./delta.json'):
    """
    input pdbbase1,pdbbase2
    output python dict [resi][resj][delta_score] as json file
    """
    import json
    dict_rrcs, not_calculated_res_list = delta_rrcs(pdbbase1, pdbbase2)
    outf_path = str(outf_path)

    with open(outf_path, 'w') as outf:
        json.dump(dict_rrcs, outf)
    print('Brief report from delta_rrcs_save2json:')
    print('output file is ./delta.json')
    print('>>>>>>>>>>>>')
    print(not_calculated_res_list[0])
    print('>>>>>>>>>>>>')
    print(not_calculated_res_list[1])
    print('>>>>>>>>>>>>')


def dict2json(python_dict: dict, outf: str):
    """
    save python dict to json,
    input: dict_object, out file path
    output: json file
    """
    import json
    if isinstance(python_dict, dict):
        with open(str(outf), 'w') as ofi:
            json.dump(python_dict, ofi)


def contact_map_from_dict(pdbbase: "pdb file") -> "np.array":
    import numpy as np

    contact_dict_a, coord_dict_a = calc_contact_no_zero(pdbbase)
    data_array = np.zeros((len(coord_dict_a), len(coord_dict_a)), dtype=np.float32)

    resi_list = coord_dict_a.keys()
    for index_i, resi in zip(range(len(resi_list)), resi_list):
        for index_j, resj in zip(range(len(resi_list)), resi_list):
            if resi in contact_dict_a:
                if resj in contact_dict_a[resi]:
                    data_array[index_i][index_j] = contact_dict_a[resi][resj]
    return data_array


if __name__ == '__main__':
    fin = sys.argv[1]
    calc_contact_save2json(fin)

    print('''
    Before using, if there is more than one chain in pdb file, 
    please make sure what you want to calculate is in right order of chain number.
    ''')
