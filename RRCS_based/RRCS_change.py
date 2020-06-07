# -*- coding: utf-8 -*-
'''
process pdb file to calculate contact score and write to csv.
    main function is prosess_save2csv()

usage: python3 RRCS_change.py file.pdb
# input .pdb
# output .csv
'''
import sys


def calc_contact(pdbbase:"pdb file")->dict:
    fi = open(pdbbase, 'r')
    all_lines = fi.readlines()
    fi.close()
    atom_lines = [l for l in all_lines if l[0:6] == 'ATOM  ']
    dict_coord = {} # dict to store coordinates. dict_coord[res][atom] = (x,y,z,occ)
    atomnum_2_name = {} # map atom number to atom name, in order to find N, CA, C, O
    contact_score = {} # dict to store final results. contact_score[ires][jres] = contact_score.
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
        res = chain_id + ':' + str(res_num) + '_' + res_name # res 
        atomnum_2_name[atom_num] = atom_name
        if res_num <= 0:
            continue
        if res not in dict_coord:
            dict_coord[res] = {}
        dict_coord[res][atom_num] = (x, y, z, occ)# store res:atom:(atom coord x,y,z,occ) in dict_coord
    
    
    #### 
    # above we read coord from .pdb file and store to a dict in hierarchy res::atom::(x,y,z,occ)
    ####
    
    #list shortening  
    num_of_aa = len(dict_coord)
    res_list = dict_coord.keys()
    ires_contact = {} ## to save interact res' contact score
    for ires in res_list: # ires is res:str , loop over every res in protein loop1
        ires_contact[ires] = [] # ires_contact:dict {res:str=[]:list}
        ires_num = int(ires.split(':')[1].split('_')[0].strip()) # got res number from ires
        for jres in res_list: # jres is res:str. for res in loop1, loop over every res in protein. 
            jres_num = int(jres.split(':')[1].split('_')[0].strip()) # got res number from jres
            if jres_num <= ires_num: # move to next step until jres_number > ires_number to avoid repeat
                continue
            jres_flag = 0
            atom_in_ires = dict_coord[ires].keys()
            atom_in_jres = dict_coord[jres].keys() 

            # 之后对每个res-res loop over atom
            for iatom in atom_in_ires:
                (ix, iy, iz, iocc) = dict_coord[ires][iatom]
                for jatom in atom_in_jres:                  
                    (jx, jy, jz, jocc) = dict_coord[jres][jatom]
                    dx = abs(ix-jx)
                    dy = abs(iy-jy)
                    dz = abs(iz-jz)
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
        for kres in ires_contact[ires]: # 取出刚才找到的相近的jres 重命名为kres
            atom_in_kres = dict_coord[kres].keys() # return a list of natoms
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
                        d2 = (ix-kx)**2 + (iy-ky)**2 + (iz-kz)**2
                        if d2 >= 21.4369:  # 4.63*4.63 = 21.4369
                            score = 0
                        elif d2 <= 10.4329:  # 3.23*3.23 = 10.4329
                            score = 1.0*iocc*kocc
                        else:
                            score = (1-(d2**0.5 - 3.23)/1.4)*iocc*kocc
                        total_score = total_score + score
                        
                        # 侧链穷尽atom对，距离在3.23:3.23~4.63:4.63打分，计入总分

            elif abs(ires_num - kres_num) > 4:
                for iatom in atom_in_ires:
                    (ix, iy, iz, iocc) = dict_coord[ires][iatom]
                    for katom in atom_in_kres:
                        (kx, ky, kz, kocc) = dict_coord[kres][katom]
                        d2 = (ix-kx)**2 + (iy-ky)**2 + (iz-kz)**2
                        if d2 >= 21.4369:  # 4.63*4.63 = 21.4369
                            score = 0
                        elif d2 <= 10.4329:  # 3.23*3.23 = 10.4329
                            score = 1.0*iocc*kocc
                        else:
                            score = (1-(d2**0.5 - 3.23)/1.4)*iocc*kocc
                        total_score = total_score + score
            contact_score[ires][kres] = total_score
    return contact_score

## 当所有atom距离都大于 4.63 是不被记录的，应该计为0

# pack contact score into csv file                                 
def prosess_save2csv(pdbbase:'pdb file'):
    '''
    process pdb file to calculate contact score and write to csv
    '''
    contact = calc_contact(pdbbase)
    outf = pdbbase + '.cscore.csv'
    fout = open(outf, 'w')
    for a_res in contact:
        b_res_list = contact[a_res].keys()
        for b_res in b_res_list:
            score = contact[a_res][b_res]
            if score > 0:
                #fout.write('%-12s\t%-12s%10.6f\n' %(a_res, b_res, score))
                fout.write('%-12s,%s,%10.6f\n' %(a_res, b_res, score))
    fout.close()

def deltaRRCS_save2csv(pdbbase1:'pdb file',pdbbase2:'pdb file'):
    '''
    delta RRCS

    calculate delta RRCS between two states of receptor, i.e. active inactive.

    only calculate

    save intersection of two protein 

    residue number is same, but residue name(because mutation) atom number and chain number maybe different.
    '''
    contact_a = calc_contact(pdbbase1)
    contact_b = calc_contact(pdbbase2)

    # make intersection of ires for checking
    list_ires_a = contact_a.keys()
    list_ires_b = contact_b.keys()
    
    set_ires = set(list_ires_a)&set(list_ires_b)

    #
    for ires_a in list_ires_a:
        if ires_a in set_ires:
            list_jres_a = contact_a[ires_a].keys()
            list_jres_b = contact_b[ires_a].keys()









if __name__=='__main__':
    

    fin = sys.argv[1]
    prosess_save2csv(fin)

    print('''
    Before using, if there is more than one chain in pdb file, 
    please make sure what you want to calculate is in right order of chain number.
    ''')