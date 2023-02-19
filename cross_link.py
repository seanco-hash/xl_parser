import argparse
import multiprocessing
import subprocess
import os
import numpy as np
import general_utils
from os import listdir
import alpha_fold_files
import pdb_files_manager
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB.DSSP import make_dssp_dict
import copy
import re
import data_table_parser
import general_xl_parser
from os.path import isfile
import csv
import pandas as pd

SS_DICT = {'H': 0, 'E': 1, 'B': 2, 'G': 3, 'I': 4, 'T': 5, 'S': 6, '-': 7}
SS_HARD_DICT = {'H': 0, 'E': 1, 'B': 2, 'G': 2, 'I': 2, 'T': 2, 'S': 2, '-': 2}

DSSO_SPACER = 10.3
DSS_SPACER = 11.4
BDP_NHP_SPACER = 34.4
LEIKER_SPACER = 9.3
QPIR_SPACER = 30
A_DSBSO_SPACER = 14
LINKER_DICT = {'DSSO': 1, 'BDP_NHP': 2, 'DSS': 3, 'UNKNOWN': 0, 'BDP-NHP': 2, '': 0, 'LEIKER': 4, None: 0, 'QPIR': 5,
               'A-DSBSO': 6}
SPACER_DICT = {'DSSO': DSSO_SPACER, 'DSS': DSS_SPACER, 'BDP_NHP': BDP_NHP_SPACER, 'BDP-NHP': BDP_NHP_SPACER,
               'UNKNOWN': 0, '': 0, 'LEIKER': LEIKER_SPACER, None: 0, 'QPIR': QPIR_SPACER, 'A-DSBSO': A_DSBSO_SPACER}
LINKER_AA_DICT = {'DSSO': ['LYS', 'N']}
PEP_A = 0
POS_IN_PEP_A = 1
RES_NUM_A = 2
UNIPORT_A = 3
PEP_B = 4
POS_IN_PEP_B = 5
RES_NUM_B = 6
UNIPORT_B = 7
LINKER = 8
_PDB = 9
ORIGIN_DB = 10
DISTANCE = 11
ERROR = 12
PDB_TYPE = 13
XL_TYPE = 14
RES_A_TYPE = 15
RES_B_TYPE = 16

IN_RES_A = 0
IN_CHAIN_A = 1
IN_RES_B = 2
IN_CHAIN_B = 3
IN_PDB_PATH = 4
IN_DISTANCE = 5
IN_CB_DISTANCE = 6
IN_LINKER = 7
IN_XL_TYPE = 8

OBJ_DIR = '/cs/labs/dina/seanco/xl_parser/obj/xl_objects/'
OBJ_DIR_PREFIX = 'xl_objects/'
PROCESSED_OBJ_DIR = '/cs/labs/dina/seanco/xl_parser/obj/processed/'
PROCESSED_OBJ_DIR_PREFIX = 'processed/'
MISSING_RESIDUES_SCRIPT = '/cs/staff/dina/utils/missing_res'

INTRA_XL = 0
INTER_XL = 1

INVALID_VAL = -1
INVALID_ERROR_VALUE = 100
MEAN_PAE_ERROR = 7


class CrossLink:
    def __init__(self, pep_a, pos_in_pep_a, res_num_a, uniport_a, pep_b, pos_in_pep_b, res_num_b,
                 uniport_b, linker_type, pdb_file=None, db="", pdb_path='', chain_a='', chain_b='', distance=0,
                 cb_distance=0, xl_type=-1, error=INVALID_ERROR_VALUE, res_a_type=None, res_b_type=None):
        self.pep_a = pep_a.upper()
        self.pep_b = pep_b.upper()
        self.pos_in_pep_a = str(pos_in_pep_a)
        self.pos_in_pep_b = str(pos_in_pep_b)
        self.res_num_a = str(res_num_a)
        self.res_num_b = str(res_num_b)
        self.uniport_a = uniport_a
        self.uniport_b = uniport_b
        self.linker_type = linker_type
        self.pdb_file = pdb_file
        self.origin_db = db
        self.distance = distance
        self.error = error
        self.pdb_type = -1
        self.res_a_type = res_a_type
        self.res_b_type = res_b_type
        self.cb_distance = cb_distance
        self.omega = 0  # cos, sin - symmetric, len == 2
        self.theta = 0  # cos(a_to_b), sin(a_to_b), cos(b_to_a), sin(b_to_a), len == 4
        self.phi = 0  # cos(a_to_b), sin(a_to_b), cos(b_to_a), sin(b_to_a), len == 4
        self.chain_a = chain_a
        self.chain_b = chain_b
        self.pdb_path = pdb_path
        self.xl_type = xl_type
        if distance > 0:
            if chain_a != chain_b:
                self.xl_type = INTER_XL
            else:
                self.xl_type = INTRA_XL

    def cross_link_obj_to_list(self):
        lst = [self.pep_a, self.pos_in_pep_a, self.res_num_a, self.uniport_a, self.pep_b,
               self. pos_in_pep_b, self.res_num_b, self.uniport_b, LINKER_DICT[self.linker_type],
               self.pdb_file]
        if hasattr(self, 'origin_db'):
            _origin_db = self.origin_db
        else:
            _origin_db = ""
        if hasattr(self, 'pdb_type'):
            _pdb_type = self.pdb_type
        else:
            _pdb_type = -1
        if hasattr(self, 'xl_type'):
            _xl_type = self.xl_type
        else:
            _xl_type = -1
        if hasattr(self, 'res_a_type'):
            res_a_type = self.res_a_type
        else:
            res_a_type = None
        if hasattr(self, 'res_b_type'):
            res_b_type = self.res_b_type
        else:
            res_b_type = None
        lst += [_origin_db, self.distance, self.error, _pdb_type,
                _xl_type, res_a_type, res_b_type]
        return lst

    def get_keys(self):
        key1 = self.uniport_a + '_' + self.res_num_a + '_' + self.uniport_b + '_' + self.res_num_b + '_' + \
               str(LINKER_DICT[self.linker_type])
        key2 = self.uniport_b + '_' + self.res_num_b + '_' + self.uniport_a + '_' + self.res_num_a + '_' + \
               str(LINKER_DICT[self.linker_type])
        if int(self.res_num_a) < 0:
            print("residue number does not found: ", key1)
        return key1, key2

    @staticmethod
    def get_opposite_key(key):
        k = key.split('_')
        new_key = k[2] + '_' + k[3] + '_' + k[0] + '_' + k[1] + '_' + k[4]
        return new_key

    @staticmethod
    def cross_links_objects_to_list(cross_links):
        xl_list = []
        for xl in cross_links:
            xl_list.append(xl.cross_link_obj_to_list())
        return xl_list

    @staticmethod
    def load_all_xl_objects_as_list(xl_objects=None):
        if xl_objects is None:
            cross_links = CrossLink.load_all_xl_objects()
        else:
            cross_links = xl_objects
        cross_links = CrossLink.cross_links_objects_to_list(cross_links)
        return cross_links

    @staticmethod
    def load_all_xl_objects(dir_=OBJ_DIR, dir_prefix=OBJ_DIR_PREFIX, specific_file=None):
        cross_link_objects = []
        for file in listdir(dir_):
            if specific_file is None or specific_file == file:
                cur_xl = general_utils.load_obj(dir_prefix + file.split('.')[0])
                cross_link_objects += cur_xl
        return cross_link_objects

    @staticmethod
    def load_all_xl_objects_as_dict(dir_=OBJ_DIR, dir_prefix=OBJ_DIR_PREFIX, specific_file=None, xl_list=None):
        xl_dict = {}
        if xl_list is None:
            xl_list = CrossLink.load_all_xl_objects(dir_, dir_prefix, specific_file)
        already_in_dict = set()
        for cl in xl_list:
            key1, key2 = cl.get_keys()
            if not(key1 in already_in_dict or key2 in already_in_dict):
                xl_dict[key1] = cl
                already_in_dict.add(key1)
        return xl_dict

    @staticmethod
    def create_xl_objects_from_csv(file_path):
        """
        Creates xl_objects list from file
        :param file_path: each row represents a object. format: aa_residue_number_a, chain_a, aa_residue_number_b,
         chain_b, pdb_file_path, ca_distance, cb_distance, linker_type
        :return: xl_objects list
        """
        with open(file_path, 'r') as f:
            reader = csv.reader(f, delimiter=",")
            objects = []
            for row in reader:
                tmp_obj = CrossLink('', -1, row[IN_RES_A], '', '', -1, row[IN_RES_B], '', row[IN_LINKER],
                                    '', '', row[IN_PDB_PATH], row[IN_CHAIN_A], row[IN_CHAIN_B],
                                    int(row[IN_DISTANCE]), int(row[IN_CB_DISTANCE]), int(row[IN_XL_TYPE]))
                objects.append(tmp_obj)
            return objects

    def validate_xl(self, polypep_list):
        length = sum(len(p) for p in polypep_list)
        err = alpha_fold_files.get_pred_align_err(self.uniport_a, self.res_num_a, self.res_num_b, length)
        if err == -1:
            err = INVALID_ERROR_VALUE
            self.error = err
            return 1
        self.error = err
        return 0

    def fix_peptides_remove_bracelet(self):
        self.pep_a = re.sub("[\(\[].*?[\)\]]", "", self.pep_a)
        self.pep_b = re.sub("[\(\[].*?[\)\]]", "", self.pep_b)

    @staticmethod
    def search_pep_in_chain(pep, chain):
        i, j = 0, 0
        while i < len(chain):
            while pep[j] == seq1(chain.child_list[i + j].resname):
                j += 1
                if j == len(pep):
                    return i
            i += 1
            j = 0

    @staticmethod
    def search_residue_in_chain(pep, pos_in_pep, chain):
        if pep != '' and int(pos_in_pep) >= 0 and chain != 0:
            i = int(pos_in_pep)
            right_aa = min(i + 5 - min(2, i), len(pep))
            left_aa = max(i - (4 - min(2, len(pep) - 1 - i)), 0)
            short_pep = pep[left_aa: right_aa]
            chain_seq = ''.join([seq1(c.resname) for c in chain.child_list])
            if short_pep.upper() in chain_seq:
                start_idx = CrossLink.search_pep_in_chain(short_pep.upper(), chain)
                res_char = pep[int(pos_in_pep)]
                new_res = chain.child_list[start_idx + (i - left_aa)]
                if res_char == seq1(new_res.resname):
                    return new_res
        return 0

    def check_res_type_by_linker(self, res):
        if res.resname == 'LYS' or res.id[1] == 1:
            return True
        return False

    def check_residues(self, residues, chain_a, chain_b):
        res_a, res_b = residues
        a, b = False, False
        if res_a != 0 and res_b != 0 and res_a.id != res_b.id:
            if self.check_res_type_by_linker(res_a):
                a = True
            elif self.pep_a != '' and 0 <= int(self.pos_in_pep_a) < len(self.pep_a):
                res_a_char = self.pep_a[int(self.pos_in_pep_a)]
                if res_a_char == seq1(res_a.resname):
                    a = True
                else:
                    res_a = CrossLink.search_residue_in_chain(self.pep_a, self.pos_in_pep_a, chain_a)
                    if res_a != 0:
                        residues[0] = res_a
                        self.res_num_a = str(res_a.id[1])
                        a = True
            if self.check_res_type_by_linker(res_b):
                b = True
            elif self.pep_b != '' and 0 <= int(self.pos_in_pep_b) < len(self.pep_b):
                res_b_char = self.pep_b[int(self.pos_in_pep_b)]
                if res_b_char == seq1(res_b.resname):
                    b = True
                elif a:
                    res_b = CrossLink.search_residue_in_chain(self.pep_b, self.pos_in_pep_b, chain_b)
                    if res_b != 0:
                        residues[1] = res_b
                        self.res_num_b = str(res_b.id[1])
                        b = True
        return a and b

    def get_residues_from_single_chain(self, chains):
        self.xl_type = INTRA_XL
        chain = list(chains)[0]
        key = (' ', int(self.res_num_a), ' ')
        if key in chain:
            res_a = chain[key]
        else:
            res_a = CrossLink.search_residue_in_chain(self.pep_a, self.pos_in_pep_a, chain)
            if res_a != 0:
                self.res_num_a = str(res_a.id[1])
                self.chain_a = chain.id
        key = (' ', int(self.res_num_b), ' ')
        if key in chain:
            res_b = chain[key]
        else:
            res_b = CrossLink.search_residue_in_chain(self.pep_b, self.pos_in_pep_b, chain)
            if res_b != 0:
                self.res_num_b = str(res_b.id[1])
                self.chain_b = chain.id
        return res_a, res_b

    @staticmethod
    def get_residue(chain, res_num, pep, pos_in_pep):
        key = (' ', int(res_num), ' ')
        if key in chain:
            res = chain[key]
        else:
            res = CrossLink.search_residue_in_chain(pep, pos_in_pep, chain)
            if res != 0:
                res_num = str(res.id[1])
        return res, res_num

    def process_single_xl(self, polypep_list=None, multichain=False, chains=None, error_objects=None):
        res_a, res_b = 0, 0
        chain_a, chain_b = 0, 0
        if multichain:
            i = 0
            while i < len(chains) and (res_a == 0 or res_b == 0):
                chain = chains[i]
                if chain.id == self.chain_a:
                    res_a, self.res_num_a = CrossLink.get_residue(chain, self.res_num_a, self.pep_a, self.pos_in_pep_a)
                    chain_a = chain
                if chain.id == self.chain_b:
                    res_b, self.res_num_b = CrossLink.get_residue(chain, self.res_num_b, self.pep_b, self.pos_in_pep_b)
                    chain_b = chain
                i += 1
            if chain_a != 0 and chain_b != 0 and chain_a.id == chain_b.id:
                self.xl_type = INTRA_XL
            else:
                self.xl_type = INTER_XL
        else:
            res_a, res_b = self.get_residues_from_single_chain(chains)

        residues = [res_a, res_b]
        check_res = self.check_residues(residues, chain_a, chain_b)
        if not check_res:
            if error_objects is not None:
                error_objects.append(self)
            self.distance = INVALID_VAL
            return 1
        res_a, res_b = residues
        self.distance = CrossLink.get_atom_distance_form_residues(res_a, res_b, 'CA', 'CA')
        self.cb_distance = CrossLink.get_atom_distance_form_residues(res_a, res_b, 'CB', 'CB')
        if self.distance == INVALID_VAL:
            return 1
        self.phi = CrossLink.get_phi(res_a, res_b)
        self.theta = CrossLink.get_theta(res_a, res_b)
        self.omega = CrossLink.get_omega(res_a, res_b)
        self.res_a_type = res_a.resname
        self.res_b_type = res_b.resname
        return 0

    @staticmethod
    def get_closest_residues(chains, res_num, chain_id, k=5):
        res = 0
        chain_res = pdb_files_manager.get_chain_by_chain_id(chains, chain_id)
        key = (' ', int(res_num), ' ')
        res = chain_res[key]
        closest = sorted(chain_res, key=lambda res_i: res['CA'] - res_i['CA'])[1:k + 1]
        closest = [(c.id[1], seq1(c.resname)) for c in closest]
        return closest

    @staticmethod
    def get_closest_res_dict():
        return general_utils.load_obj("closest_residues")

    @staticmethod
    def create_closest_residues_dict(xl_objects, k=10):
        xl_objects = sorted(xl_objects, key=lambda o: o.pdb_path)
        prev_pdb = ''
        pdb_parser = PDBParser(PERMISSIVE=1)
        closest_dict = dict()
        objects_for_ablations = []
        for i, obj in enumerate(xl_objects):
            try:
                if obj.pdb_path != prev_pdb:
                    closest_dict[obj.pdb_path] = dict()
                    if obj.pdb_path[-3:] == 'cif':
                        continue
                    else:
                        structure = pdb_parser.get_structure(obj.pdb_path.split('/')[-1], obj.pdb_path)
                    chains = list(structure.get_chains())
                obj.chain_a = obj.chain_a if obj.chain_a != '' else 'A'
                obj.chain_b = obj.chain_b if obj.chain_b != '' else 'A'
                for r, c in [(obj.res_num_a, obj.chain_a), (obj.res_num_b, obj.chain_b)]:
                    if c not in closest_dict[obj.pdb_path]:
                        closest_dict[obj.pdb_path][c] = dict()
                    if r not in closest_dict[obj.pdb_path][c]:
                        closest_dict[obj.pdb_path][c][r] = CrossLink.get_closest_residues(chains, r, c, k)
            except Exception as e:
                prev_pdb = obj.pdb_path
                print(f"problem with object with pdb_file: {obj.pdb_path}, chains: {obj.chain_a}, {obj.chain_b}\n {e}")
                continue
            prev_pdb = obj.pdb_path
            objects_for_ablations.append(obj)
            if i % 100 == 0:
                print(i)
        general_utils.save_obj(closest_dict, "closest_residues")
        general_utils.save_obj(objects_for_ablations, "ablation_objects")
        print(len(objects_for_ablations))

    @staticmethod
    def create_mutation_pdbs():
        xl_objects = general_utils.load_obj('ablation_objects')
        xl_objects = sorted(xl_objects, key=lambda o: o.pdb_path)
        prev_pdb = ''
        pdb_parser = PDBParser(PERMISSIVE=1)
        closest_dict = general_utils.load_obj('closest_residues')
        processes = []
        for i, obj in enumerate(xl_objects):
            # if obj.pdb_path[-3:] != 'ent':
            #     continue
            if obj.pdb_path != prev_pdb:
                seq_dict = {}
                if obj.pdb_path[-3:] == 'cif':
                    continue
                else:
                    structure = pdb_parser.get_structure(obj.pdb_path.split('/')[-1], obj.pdb_path)
                chains = list(structure.get_chains())
            seq_a, seq_b = pdb_files_manager.get_sequences_for_obj(obj, chains, seq_dict)
            pref, suff = obj.pdb_path.split('.')
            closest_a = closest_dict[obj.pdb_path][obj.chain_a][obj.res_num_a]
            pdb_a = obj.pdb_path if len(chains) == 1 else f"{pref}_{obj.chain_a}.{suff}"
            # p1 = multiprocessing.Process(target=pdb_files_manager.create_scwrl_mutation_files,
            #                             args=(obj.res_num_a, obj.chain_a, closest_a, pdb_a, copy.deepcopy(seq_a),
            #                                               pdb_files_manager.SCWRL_OUT_DIR))
            # p1.start()
            # processes.append(p1)
            pdb_files_manager.create_scwrl_mutation_files(obj.res_num_a, obj.chain_a, closest_a, pdb_a, copy.deepcopy(seq_a),
                                                          pdb_files_manager.SCWRL_OUT_DIR)
            closest_b = closest_dict[obj.pdb_path][obj.chain_b][obj.res_num_b]
            pdb_b = pdb_a if obj.chain_a == obj.chain_b else f"{pref}_{obj.chain_b}.{suff}"
            # p2 = multiprocessing.Process(target=pdb_files_manager.create_scwrl_mutation_files,
            #                              args=(obj.res_num_b, obj.chain_b, closest_b, pdb_b, copy.deepcopy(seq_b),
            #                                    pdb_files_manager.SCWRL_OUT_DIR))
            # p2.start()
            # processes.append(p2)
            pdb_files_manager.create_scwrl_mutation_files(obj.res_num_b, obj.chain_b, closest_b, pdb_b, copy.deepcopy(seq_b),
                                                          pdb_files_manager.SCWRL_OUT_DIR)
            prev_pdb = obj.pdb_path
        # for p in processes:
        #     p.join()


    @staticmethod
    def get_atom_distance_form_residues(res_a, res_b, atom_in_a, atom_in_b):
        if isinstance(res_a, int) or isinstance(res_b, int) or atom_in_a not in res_a or atom_in_b not in res_b:
            return INVALID_VAL
        # if isinstance(res_a[atom_in_a], int) or isinstance(res_b[atom_in_b], int):
        #     return INVALID_VAL
        return res_a[atom_in_a] - res_b[atom_in_b]

    @staticmethod
    def get_theta(res_a, res_b):
        angles = np.zeros(4)
        if not res_a.has_id('CB') or not res_b.has_id('CB'):  # GLY OR UNK OR MISSING CB
            return None

        angle = PDB.calc_dihedral(res_a["N"].get_vector(), res_a["CA"].get_vector(),
                                  res_a["CB"].get_vector(), res_b["CB"].get_vector())
        angles[0] = np.cos(angle)
        angles[1] = np.sin(angle)
        angle = PDB.calc_dihedral(res_b["N"].get_vector(), res_b["CA"].get_vector(),
                                  res_b["CB"].get_vector(), res_a["CB"].get_vector())
        angles[2] = np.cos(angle)
        angles[3] = np.sin(angle)
        # angle_a_to_b = np.degrees(angle) % 360
        return angles

    @staticmethod
    def get_phi(res_a, res_b):
        angles = np.zeros(4)
        if not res_a.has_id('CB') or not res_b.has_id('CB'):  # GLY OR UNK OR MISSING CB
            return None

        angle = PDB.calc_angle(res_a["CA"].get_vector(), res_a["CB"].get_vector(), res_b["CB"].get_vector())
        angles[0] = np.cos(angle)
        angles[1] = np.sin(angle)
        angle = PDB.calc_angle(res_b["CA"].get_vector(), res_b["CB"].get_vector(), res_a["CB"].get_vector())
        angles[2] = np.cos(angle)
        angles[3] = np.sin(angle)
        # angle_a_to_b = np.degrees(angle) % 360
        return angles

    @staticmethod
    def get_omega(res_a, res_b):
        angles = np.zeros(2)
        if not res_a.has_id('CB') or not res_b.has_id('CB'):  # GLY OR UNK OR MISSING CB
            return None

        angle = PDB.calc_dihedral(res_a["CA"].get_vector(), res_a["CB"].get_vector(), res_b["CB"].get_vector(),
                                  res_b["CA"].get_vector())
        angles[0] = np.cos(angle)
        angles[1] = np.sin(angle)
        # angle_a_to_b = np.degrees(angle) % 360
        return angles

    @staticmethod
    def unit_pickle_files():
        xl_objects = CrossLink.load_all_xl_objects(PROCESSED_OBJ_DIR, PROCESSED_OBJ_DIR_PREFIX)
        for file in listdir(PROCESSED_OBJ_DIR):
            os.remove(PROCESSED_OBJ_DIR + file)
        general_xl_parser.XlParser.save_cross_links(xl_objects, None,
                                                    "processed_xl_objects", PROCESSED_OBJ_DIR_PREFIX)

    @staticmethod
    def clear_dup_before_process(xl_objects):
        processed = CrossLink.load_all_xl_objects_as_dict(PROCESSED_OBJ_DIR, PROCESSED_OBJ_DIR_PREFIX)
        filtered = []
        for obj in xl_objects:
            key1, key2 = obj.get_keys()
            if key1 in processed or key2 in processed:
                continue
            filtered.append(obj)
        return filtered

    @staticmethod
    def process_all_xl(xl_objects, uni_chain_dict=None, pid=None, out_name=None, already_processed=False,
                       inter_pdbs=None):
        # xl_objects = CrossLink.clear_dup_before_process(xl_objects)
        print(f"number of xl objects to process: {len(xl_objects)}\n")
        xl_objects.sort(key=lambda x: x.uniport_a)
        counters = np.zeros(5)
        cif_files = general_utils.load_obj('cif_files')
        prev_obj = CrossLink('', 0, 0, '', '', 0, 0, '', '')
        prev_pdb = '-1'
        no_file_error = 0
        no_file_inter_xl = 0
        parser_exceptions = 0
        same_residue = 0
        fasta_dict = general_utils.load_obj('fasta_files_dict_by_uniport')
        distance_exceptions = 0
        error_exceptions = 0
        pdb_parser = PDBParser(PERMISSIVE=1)
        cif_parser = MMCIFParser()
        cif_files_errors = []
        pdb_files_errors = []
        mmcif_dict = None
        chains = None
        # processed_objects = CrossLink.load_all_xl_objects(PROCESSED_OBJ_DIR, PROCESSED_OBJ_DIR_PREFIX)
        processed_objects = []
        error_objects = []
        error_multichain_objects = []
        entity_chain_dict = pdb_files_manager.get_pdb_entity_chain_dict()
        i = 0
        for cur_obj in xl_objects:
            i += 1
            if i > 1000:
                print(f"processed by now: {len(processed_objects)}")
                i = 0
            pdb_file_name = pdb_files_manager.find_pdb_file_from_xl_obj(cur_obj, cif_files, inter_pdbs=inter_pdbs)
            if pdb_file_name != prev_pdb:
                if pdb_file_name is None:
                    if cur_obj.uniport_a == cur_obj.uniport_b:
                        no_file_error += 1
                    else:
                        no_file_inter_xl += 1
                    continue
                try:
                    if pdb_file_name[-3:] == 'cif':
                        structure = cif_parser.get_structure(pdb_file_name.split('/')[-1], pdb_file_name)
                        mmcif_dict = cif_parser.get_mmcif_dict()
                        if mmcif_dict is None or structure is None:
                            cif_files_errors.append(pdb_file_name)
                            continue
                    else:
                        structure = pdb_parser.get_structure(pdb_file_name.split('/')[-1], pdb_file_name)
                        if structure is None:
                            pdb_files_errors.append(pdb_file_name)
                            continue
                        mmcif_dict = None
                except Exception as e:
                    print(e)
                    parser_exceptions += 1
                    continue
            elif cur_obj.res_num_a == prev_obj.res_num_a and cur_obj.res_num_b == prev_obj.res_num_b and \
                    cur_obj.uniport_a == prev_obj.uniport_a and cur_obj.uniport_b == prev_obj.uniport_b:
                continue
            if cur_obj.uniport_a == cur_obj.uniport_b and \
                    (uni_chain_dict is None or (cur_obj.uniport_a, cur_obj.uniport_b) not in uni_chain_dict):
                multichain = False
                if cur_obj.res_num_a == cur_obj.res_num_b:
                    same_residue += 1
                    continue
            else:
                multichain = True
            if pdb_file_name is not None:
                if cur_obj.uniport_a == cur_obj.uniport_b and not multichain:
                    res_new = cur_obj.process_single_xl(None, False, structure.get_chains(), error_objects)
                    if res_new > 0:
                        distance_exceptions += res_new
                        continue
                    error_exceptions += cur_obj.validate_xl(structure.get_chains())
                else:
                    tmp_obj = pdb_files_manager.process_multichain_xl_object(cur_obj, chains, fasta_dict, structure,
                                                                             45, counters, uni_chain_dict,
                                                                             pdb_file_name, mmcif_dict,
                                                                             entity_chain_dict, error_objects, already_processed)
                    if tmp_obj is not None:
                        cur_obj = tmp_obj
                    else:
                        error_multichain_objects.append(cur_obj)
                        continue
                cur_obj.pdb_path = pdb_file_name
                processed_objects.append(cur_obj)
            prev_obj = cur_obj
            prev_pdb = pdb_file_name

        if inter_pdbs is not None:
            print("Retry failed multicahin objects with AF")
            retry_objects = [obj for obj in error_multichain_objects if obj.uniport_a == obj.uniport_b]
            process_failed_with_af = CrossLink.process_all_xl(retry_objects)
            processed_objects = processed_objects + process_failed_with_af
            print("Retry finished")

        print("xl with missing files: ", no_file_error)
        print("inter xl missing file: ", no_file_inter_xl)
        print("parsing errors: ", parser_exceptions)
        print("distance problems single chain", distance_exceptions)
        print("error problems single chain ", error_exceptions)
        print(f"objects with same residues: {same_residue} \n")
        print(f"Successfully processed multichain: {counters[0] + counters[4]}")
        print(f"Multiple options fo xl multichain: {counters[4]}")
        print(f"xl multichain chains found, but no valid xl: {counters[1]}")
        print(f"missing fasta for multichain {counters[2]}")
        print(f"No chain found for multichain: {counters[3]}")
        print("successfully processed: ", len(processed_objects))
        print(cif_files_errors)
        print(pdb_files_errors)
        print(f"files with potential numbering problem: {len(error_objects)}")
        print(f"objects failed with multichain: {len(error_multichain_objects)}")

        if pid is None:
            if out_name is not None:
                general_utils.save_obj(error_multichain_objects, 'error_multichain_objects')
                general_utils.save_obj(error_objects, 'numbering_problem_objects' + out_name)
                general_xl_parser.XlParser.save_cross_links(processed_objects, None,
                                                            out_name, PROCESSED_OBJ_DIR_PREFIX)
        else:
            print("Process number: ", pid)
            general_xl_parser.XlParser.save_cross_links(processed_objects, None,
                                                        out_name + str(pid), PROCESSED_OBJ_DIR_PREFIX)
        return processed_objects

    @staticmethod
    def parallel_process_all_xl(xl_objects, uni_chain_dict=None):
        try:
            with multiprocessing.Manager() as manager:
                available_cpus = len(os.sched_getaffinity(0)) - 1
                print("available cpus: ", available_cpus, flush=True)
                start_idx = 0
                step = int(len(xl_objects) / (available_cpus - 1))
                xl_objects = manager.list(xl_objects)
                processes = []
                for pid in range(available_cpus):
                    end = min(start_idx + step, len(xl_objects))
                    p = multiprocessing.Process(target=CrossLink.process_all_xl,
                                                args=(xl_objects[start_idx: end], uni_chain_dict, pid))
                    p.start()
                    processes.append(p)
                    start_idx += step
                for p in processes:
                    p.join()
        except Exception as e:
            print("Thread_Error")
            print(e)

    @staticmethod
    def unit_processed_objects(list1='processed_xl_objects_inter', list2='processed_xl_objects',
                               new_list='processed_xl_objects_inter_intra'):
        l1 = general_utils.load_obj(list1, PROCESSED_OBJ_DIR)
        l2 = general_utils.load_obj(list2, PROCESSED_OBJ_DIR)
        new_l = l1 + l2
        general_xl_parser.XlParser.save_cross_links(new_l, None, new_list, PROCESSED_OBJ_DIR_PREFIX)

    @staticmethod
    def filter_xl_obj(xl_objects, max_dist=45, pred_err=15, linker_type=None, invalid_error=True):
        filtered = []
        for obj in xl_objects:
            if 0 < obj.distance <= max_dist:
                    # (obj.error <= pred_err or (invalid_error and obj.error == INVALID_ERROR_VALUE)):
                if linker_type == "any but none" and LINKER_DICT[obj.linker_type] != 0:
                    filtered.append(obj)
                elif linker_type is None or linker_type == obj.linker_type:
                    filtered.append(obj)
        print("filtered: ", len(filtered))
        return filtered

    @staticmethod
    def find_unknown_linkers(obj_as_list):
        objects = np.array(obj_as_list)
        u_dict = dict()
        unknown = objects[objects[:, LINKER] == '0']
        for obj in unknown:
            db = obj[ORIGIN_DB]
            dbs = db.split(',')
            for d in dbs:
                if d in u_dict:
                    u_dict[d] += 1
                else:
                    u_dict[d] = 1
        u_dict = sorted(u_dict.items(), key=lambda it: it[1])
        print(u_dict)

    @staticmethod
    def update_linker_by_origin_db(xl_objects):
        upd_objects = []
        db_linker_dict = {'Leiker2016Ecoli_Lei': 'LEIKER', 'HeLa_PTX_titration_Cell_Reports_2019': 'BDP-NHP',
                          'iqPIR_HeLa_17AAG_Chavez2020_Bruce': 'BDP-NHP',
                          'ChavezChemBiol2016_BruceLab': 'BDP-NHP',
                          'KellerJProteomeRes2019_HeLaInVivo_Bruce': 'BDP-NHP',
                          'Schweppe2015ChemBio_BruceLab': 'BDP-NHP',
                          'KellerJProteomeRes2019_HeLaLysate_Bruce': 'BDP-NHP',
                          'HeLaNatureCommunications2015_BruceLab': 'BDP-NHP',
                          'SchweppemousemitoPNAS2017_BruceLab': 'BDP-NHP',
                          'ChavezCellSystems2017_BruceLab': 'BDP-NHP', 'MangoNuclear2018_Bruce': 'BDP-NHP',
                          'XLdb2013_Malmstrom': 'DSS', 'Ecoli2012_BingYang': 'DSS',
                          'Liu2017MCPmousemito_Heck': 'DSSO', 'Chavezetal_MolCellProteomics2013': 'BDP-NHP',
                          'Weisbrodetal_JProteomeRes2013': 'BDP-NHP',
                          'Crcrosslinking_Lesliehicks': 'DSSO', 'ACP_Percoll_Fraction_XL_Bruce': 'BDP-NHP',
                          'Caudal_iqPIR_TACsham_Bruce': 'BDP-NHP', 'MCF7_Hsp90_6plex_Bruce': 'BDP-NHP',
                          '6-plex_iqPIR_HeLa_QE_Bruce': 'BDP-NHP', '6-plex_iqPIR_5prot_mix_QE_Bruce': 'BDP-NHP',
                          'Kaake2014MCP_Huang': 'A-DSBSO'}
        changed_objects = 0
        for obj in xl_objects:
            if LINKER_DICT[obj.linker_type] == 0:
                dbs = obj.origin_db.split(',')
                for i in range(len(dbs)):
                    new_obj = copy.deepcopy(obj)
                    if dbs[i] in db_linker_dict:
                        new_obj.linker_type = db_linker_dict[dbs[i]]
                    new_obj.origin_db = dbs[i]
                    upd_objects.append(new_obj)
                    changed_objects += 1
            else:
                upd_objects.append(obj)

        return upd_objects
        # general_xl_parser.XlParser.save_cross_links(xl_objects, None,
        #                                             "processed_xl_objects", PROCESSED_OBJ_DIR_PREFIX)

    @staticmethod
    def get_only_objects_to_update(xl_objects, n_feat=5):
        none_file_names = []
        feat_dict = pdb_files_manager.get_xl_neighbors_dict()
        cif_files = general_utils.load_obj('cif_files')
        filtered = []
        uniports_included = set()
        for obj in xl_objects:
            file_name = obj.pdb_path
            if file_name is None:
                none_file_names.append(obj)
                continue
            uniport = file_name.split('.')[0]
            uniport = uniport.split('/')[-1]
            if uniport not in feat_dict or uniport in uniports_included:
                filtered.append(obj)
                uniports_included.add(uniport)
                continue
            uniport_feat_dict = feat_dict[uniport]
            if obj.res_num_a not in uniport_feat_dict or obj.res_num_b not in uniport_feat_dict:
                filtered.append(obj)
                uniports_included.add(uniport)
                continue
            res_a_feat = np.stack(uniport_feat_dict[obj.res_num_a][:-1])  # -1 - last is distance
            res_b_feat = np.stack(uniport_feat_dict[obj.res_num_b][:-1])
            if res_a_feat.shape[1] != n_feat or res_b_feat.shape[1] != n_feat:
                filtered.append(obj)
                uniports_included.add(uniport)
        return filtered

    @staticmethod
    def export_objects_to_txt_format(xl_objects, inter_pdbs=None, out_path=None, skip_exist=True, dimer_form=False):
        """
        export all the cross links from same uniport to list in a text file (chain, residue, chain, residue)
        """
        cif_files = general_utils.load_obj('cif_files')
        xl_objects.sort(key=lambda x: x.uniport_a)
        dict_by_uni = {}
        dimer_dict = None
        if dimer_form:
            dimer_dict = general_utils.load_obj('dimers_in_pdb_dict')
        for obj in xl_objects:
            if int(obj.res_num_a) == -1 or int(obj.res_num_b) == -1:
                continue
            # All crosslinks from all chains exported to same file
            pdb_name, _ = pdb_files_manager.get_obj_files_key(obj, cif_files, inter_pdbs, by_chain=False)
            if pdb_name in dict_by_uni:
                dict_by_uni[pdb_name].append(obj)
            else:
                dict_by_uni[pdb_name] = [obj]
        counters = [0, 0, 0, 0, 0]
        exist = [0]
        for key, lst in dict_by_uni.items():
            pdb_files_manager.create_single_pdb_xl_list_txt_file(lst, key, cif_files, counters, 35, inter_pdbs,
                                                                 out_path, skip_exist, exist, dimer_dict)
        print(f"number of created inter xl: {counters[0]}")
        print(f"mismatch in residue number inter_xl / illegal distance: {counters[1]}")
        print(f"missing uniport sequence: {counters[2]}")
        print(f"did'nt found chain: {counters[3]}")
        print(f"multiple options for cross link: {counters[4]}")
        print(f"already exist: {exist[0]}")

    @staticmethod
    def export_objects_to_txt_for_test(xl_objects, inter_pdbs, out_path=None, skip_exist=True):
        entity_chain_dict = pdb_files_manager.get_pdb_entity_chain_dict()
        uni_chain_dict = get_uni_chain_id_dict()
        cif_files = general_utils.load_obj('cif_files')
        xl_objects.sort(key=lambda x: x.uniport_a)
        cif_parser = MMCIFParser()
        dict_by_uni = {}
        for obj in xl_objects:
            if int(obj.res_num_a) == -1 or int(obj.res_num_b) == -1:
                continue
            pdb_name = obj.pdb_path.split('/')[-1].split('.')[0]
            optional_chain_a, optional_chain_b = [], []
            if pdb_name.split('.')[-1] != 'cif':
                pdb_files_manager.find_chain_ids_pdb(obj, entity_chain_dict, optional_chain_a, optional_chain_b,
                                                     uni_chain_dict, pdb_name)
            else:
                structure = cif_parser.get_structure(pdb_name.split('/')[-1], pdb_name)
                mmcif_dict = cif_parser.get_mmcif_dict()
                pdb_files_manager.find_chain_ids_cif(obj, mmcif_dict, optional_chain_a, optional_chain_b)
            obj.chain_a = ''.join(optional_chain_a)
            obj.chain_b = ''.join(optional_chain_b)
            if pdb_name in dict_by_uni:
                dict_by_uni[pdb_name].append(obj)
            else:
                dict_by_uni[pdb_name] = [obj]
        counters = [0, 0, 0, 0, 0]
        exist = [0]
        for key, lst in dict_by_uni.items():
            pdb_files_manager.create_single_pdb_xl_list_txt_file(lst, key, cif_files, counters, 32, inter_pdbs,
                                                                 out_path, skip_exist, exist)
        print(f"number of created inter xl: {counters[0]}")

    @staticmethod
    def create_features_from_xl_objects(xl_objects):
        """
        Currently not in use.
        reads features dict with pdb_files_manager.read_features_to_dict.
        for each xl object, join the features of each residue and the label (distance)
        :return:
        """
        feat_dict = pdb_files_manager.read_features_to_dict()
        all_features = []
        for obj in xl_objects:
            features = []
            if obj.uniport_a == obj.uniport_b and obj.uniport_a in feat_dict:
                uni_key_a = uni_key_b = obj.uniport_a
            elif obj.pdb_file in feat_dict:
                uni_key_a = uni_key_b = obj.pdb_file
            elif obj.uniport_a != obj.uniport_b and obj.uniport_a in feat_dict and obj.uniport_b in feat_dict:
                uni_key_a = obj.uniport_a
                uni_key_b = obj.uniport_b
            else:
                print(f"invalid xl object: {obj.uniport_a}, {obj.res_num_a}, {obj.uniport_b}, {obj.res_num_b}")
                continue
            features.append(feat_dict[uni_key_a][obj.res_num_a])
            features.append(feat_dict[uni_key_b][obj.res_num_b])
            features.append(SPACER_DICT[obj.linker_type])
            features.append(obj.distance)
            all_features.append(features)

        general_utils.save_obj(all_features, "xl_neighbors_features")

    @staticmethod
    def validate_processed_objects(xl_obj):
        filtered = []
        cif_files = general_utils.load_obj('cif_files')
        inter_pdbs = pdb_files_manager.get_inter_pdbs()
        for obj in xl_obj:
            file_name = obj.pdb_path
            if file_name is not None:
                filtered.append(obj)
        return filtered
        # general_xl_parser.XlParser.save_cross_links(filtered, None,
        #                                             "processed_xl_objects", PROCESSED_OBJ_DIR_PREFIX)

    @staticmethod
    def clear_duplicates_of_xl_objects_wo_linker(xl_objects):
        cleared = []
        exists = set()
        for obj in xl_objects:
            key1, key2 = obj.get_keys()
            key1 = key1[:-2]
            key2 = key2[:-2]
            if key1 in exists or key2 in exists:
                continue
            exists.add(key1)
            exists.add(key2)
            cleared.append(obj)
        print(f"No duplicated amount of objects: {len(cleared)}")
        return cleared

    @staticmethod
    def clear_duplicates_of_xl_objects_by_linker(xl_objects):
        print(f"Initial amount of objects: {len(xl_objects)}")
        obj_dict = CrossLink.load_all_xl_objects_as_dict(xl_list=xl_objects)
        to_remove = set()
        for key, obj in obj_dict.items():
            if key[-1] != '0':
                opposite_key = CrossLink.get_opposite_key(key)
                no_linker_key1 = key[:-1] + '0'
                no_linker_key2 = opposite_key[:-1] + '0'
                if no_linker_key1 in obj_dict:
                    to_remove.add(no_linker_key1)
                if no_linker_key2 in obj_dict:
                    to_remove.add(no_linker_key2)

        for key in to_remove:
            del obj_dict[key]

        cleared = [obj for obj in obj_dict.values()]
        print(f"No duplicated amount of objects: {len(cleared)}")
        return cleared

    @staticmethod
    def clear_duplicates_of_dimers(xl_objects, cif_files, inter_pdbs):
        filtered = []
        residue_pairs_by_pdb = dict()
        print(f"before filtering: {len(xl_objects)}")
        dups = 0
        for obj in xl_objects:
            pdb_file = obj.pdb_path
            key = pdb_file.split('/')[-1].split('.')[0]
            if key not in residue_pairs_by_pdb:
                residue_pairs_by_pdb[key] = set()
            if (obj.res_num_a, obj.res_num_b, obj.distance, obj.linker_type) not in residue_pairs_by_pdb[key]:
                residue_pairs_by_pdb[key].add((obj.res_num_a, obj.res_num_b, obj.distance, obj.linker_type))
                residue_pairs_by_pdb[key].add((obj.res_num_b, obj.res_num_a, obj.distance, obj.linker_type))
                filtered.append(obj)
            else:
                dups += 1
        print(f"after filtering {len(filtered)}")
        return filtered

    @staticmethod
    def clear_duplicates_of_xl_objects(xl_objects, linker=False):
        if linker:
            return CrossLink.clear_duplicates_of_xl_objects_by_linker(xl_objects)
        return CrossLink.clear_duplicates_of_xl_objects_wo_linker(xl_objects)

    @staticmethod
    def xl_objects_stats(xl_objects):
        intra_xl = 0
        inter_xl = 0
        for obj in xl_objects:
            if obj.chain_a == obj.chain_b:
                intra_xl += 1
            else:
                inter_xl += 1
        print(f"Intra_xl: {intra_xl}")
        print(f"Inter_xl: {inter_xl}")

    @staticmethod
    def fix_peptides_all_xl_objects(xl_objects):
        for obj in xl_objects:
            obj.fix_peptides_remove_bracelet()

    @staticmethod
    def map_inter_xl_objects_to_pdb(xl_objects):
        pdb_dict = pdb_files_manager.read_pdb_info_jsons()
        uniport_dict = dict()
        inter_pdbs = set()
        pdbs_frequency = dict()
        af_files = 0
        pdb_files = 0
        for pdb, entities in pdb_dict.items():
            if len(entities) > 1:
                ids = list(entities.keys())
                uniports = list(entities.values())
                for i in range(len(uniports)):
                    for j in range(len(uniports)):
                        # if i != j:
                        if uniports[i] not in uniport_dict:
                            uniport_dict[uniports[i]] = dict()
                        if uniports[j] not in uniport_dict[uniports[i]]:
                            uniport_dict[uniports[i]][uniports[j]] = []
                        uniport_dict[uniports[i]][uniports[j]].append((pdb, ids[i], ids[j]))

        for obj in xl_objects:
            if obj.uniport_a in uniport_dict:
                if obj.uniport_b in uniport_dict[obj.uniport_a]:
                    for option in uniport_dict[obj.uniport_a][obj.uniport_b]:
                        if option[0] not in pdbs_frequency:
                            pdbs_frequency[option[0]] = 0
                        pdbs_frequency[option[0]] += 1

        possible_objects = []
        chain_id_dict = dict()
        for obj in xl_objects:
            cur_pdb = None
            if obj.uniport_a in uniport_dict:
                if obj.uniport_b in uniport_dict[obj.uniport_a]:
                    if (obj.uniport_a, obj.uniport_b) not in chain_id_dict:
                        opt_pdb_freqs = [pdbs_frequency[o[0]] for o in uniport_dict[obj.uniport_a][obj.uniport_b]]
                        most_frequent = np.argmax(opt_pdb_freqs)
                        cur_pdb = uniport_dict[obj.uniport_a][obj.uniport_b][most_frequent][0]
                        for option in uniport_dict[obj.uniport_a][obj.uniport_b]:
                            if option[0] == cur_pdb:
                                if (obj.uniport_a, obj.uniport_b) not in chain_id_dict:
                                    chain_id_dict[(obj.uniport_a, obj.uniport_b)] = []
                                if obj.uniport_a != obj.uniport_b or (option[1], option[2]) not in chain_id_dict[(obj.uniport_a, obj.uniport_b)]:
                                    chain_id_dict[(obj.uniport_a, obj.uniport_b)].append((option[1], option[2]))
                                # obj.uniport_a, obj.uniport_b = uniport_dict[obj.uniport_a][obj.uniport_b][1], uniport_dict[obj.uniport_a][obj.uniport_b][2]
                        # if obj.uniport_a == obj.uniport_b:
                        #     idx_tuples = chain_id_dict[(obj.uniport_a, obj.uniport_b)]
                        #     idx_set = set()
                        #     for idx_tuple in idx_tuples:
                        #         idx_set.add(idx_tuple[0])
                        #         idx_set.add(idx_tuple[1])
                        #     for idx in idx_set:
                        #         idx_tuples.append((idx, idx))
                    else:
                        opt_pdb_freqs = [pdbs_frequency[o[0]] for o in uniport_dict[obj.uniport_a][obj.uniport_b]]
                        most_frequent = np.argmax(opt_pdb_freqs)
                        cur_pdb = uniport_dict[obj.uniport_a][obj.uniport_b][most_frequent][0]
                    obj.pdb_file = cur_pdb
                    inter_pdbs.add(cur_pdb)
                    possible_objects.append(obj)
                    pdb_files += 1
            if cur_pdb is None and obj.uniport_a == obj.uniport_b:
                possible_objects.append(obj)
                af_files += 1

        # pdb_files_manager.download_pdbs_from_xl_objects(possible_objects, inter_pdbs)
        # pdb_files_manager.create_cif_set()
        # general_utils.save_obj(chain_id_dict, 'uni_chain_id_dict')
        # general_utils.save_obj(inter_pdbs, 'inter_pdbs')
        # pdb_files_manager.parse_pdb_headers(pdb_dict=pdb_files_manager.get_pdb_entity_chain_dict())
        # pdb_files_manager.parse_pdb_headers()
        # CrossLink.process_all_xl(possible_objects, uni_chain_dict=chain_id_dict)
        return possible_objects, chain_id_dict

    @staticmethod
    def plot_ss_analysis():
        ss_info = general_utils.load_obj("secondary_structure_info")
        labels = ['H', 'B', 'L']
        ss_info[ss_info[:, 0] > 1, 0] = SS_HARD_DICT['-']
        ss_info[ss_info[:, 1] > 1, 1] = SS_HARD_DICT['-']
        labels_to_plot = []
        distances_by_ss = []
        for i in range(3):
            for j in range(i + 1):
                data1 = ss_info[ss_info[:, 0] == i]
                data1 = data1[data1[:, 1] == j, 2]
                if i != j:
                    data2 = ss_info[ss_info[:, 0] == j]
                    data2 = data2[data2[:, 1] == i, 2]
                    data = np.concatenate((data1, data2))
                else:
                    data = data1
                distances_by_ss.append(data)
                labels_to_plot.append(labels[i] + '_' + labels[j])
        lengths = np.array([len(l) / 1000 for l in distances_by_ss])
        general_utils.initialize_plt_params()
        df1 = pd.DataFrame(lengths, columns=['lens'])
        df1['labels'] = labels_to_plot
        # general_utils.box_plot(distances_by_ss, labels_to_plot)
        general_utils.plot_line(df1, None, "", "SS", "Num Samples (K)", bar=True, plot_values=False)

    @staticmethod
    def analyze_secondary_structure():
        xl_objects = general_utils.load_obj('ablation_objects')
        xl_objects = sorted(xl_objects, key=lambda o: o.pdb_path)
        ss_info = np.zeros((len(xl_objects), 5))
        for i, obj in enumerate(xl_objects):
            if not os.path.isfile(obj.pdb_path + '.dssp'):
                dssp_dict_from_pdb_file(obj.pdb_path, DSSP='/cs/staff/dina/software/Staccato/mkdssps')
            dssp_dict = make_dssp_dict(obj.pdb_path + ".dssp")[0]
            if obj.chain_a == '':
                obj.chain_a = obj.chain_b = 'A'
            try:
                aa1, ss1, acc1 = dssp_dict[obj.chain_a, (' ', int(obj.res_num_a), ' ')][:3]
                aa2, ss2, acc2 = dssp_dict[obj.chain_b, (' ', int(obj.res_num_b), ' ')][:3]
            except KeyError as e:
                try:
                    aa1, ss1, acc1 = dssp_dict['A', (' ', int(obj.res_num_a), ' ')][:3]
                    aa2, ss2, acc2 = dssp_dict['A', (' ', int(obj.res_num_b), ' ')][:3]
                except KeyError as e:
                    print(e)
                    print(f"problem with {obj.pdb_path}, {obj.chain_a}, {obj.res_num_a}, {obj.chain_b}, {obj.res_num_b}")
                    continue
            ss_info[i][0] = SS_DICT[ss1]
            ss_info[i][1] = SS_DICT[ss2]
            ss_info[i][2] = obj.distance
            ss_info[i][3] = obj.xl_type
            ss_info[i][4] = LINKER_DICT[obj.linker_type]
            if i % 1000 == 0:
                general_utils.save_obj(ss_info, "secondary_structure_info")
        general_utils.save_obj(ss_info, "secondary_structure_info")


def get_uni_chain_id_dict():
    return general_utils.load_obj('uni_chain_id_dict')


def get_numbering_problem_objects():
    return general_utils.load_obj('numbering_problem_objects')


def linker_type_histogram(xl_objects):
    linker_types = [LINKER_DICT[obj.linker_type] for obj in xl_objects]
    bins = [1, 2, 3, 4, 5, 6]
    xticks_values = [1.5, 2.5, 3.5, 4.5, 5.5]
    xticks = ['DSSO', 'BDP_NHP', 'DSS', 'LEIKER', 'A-DSBSO']
    general_utils.plot_histogram(linker_types, "Crosslinker Type Distribution", x_label="Linker Type", xticks_values=xticks_values,
                                 xticks_names=xticks, bins=bins)


def linker_type_stats(xl_objects):
    counts = dict()
    for obj in xl_objects:
        if LINKER_DICT[obj.linker_type] in counts:
            counts[LINKER_DICT[obj.linker_type]] += 1
        else:
            counts[LINKER_DICT[obj.linker_type]] = 1
    for key, val in counts.items():
        print(f"{key}: {(val / len(xl_objects)) * 100}")
    print(f"total xl: {len(xl_objects)}")


def analyze_predicted_align_error(xl_objects):
    count = 0
    errors = []
    for obj in xl_objects:
        if obj.error == INVALID_ERROR_VALUE or obj.error == 0 or obj.error > 33:
            count += 1
            print(f"{obj.res_num_a}, {obj.res_num_b}")
        else:
            errors.append(obj.error)
    print(count)
    print(np.average(errors))
    general_utils.plot_histogram([errors], 'Predicted align error histogram')


def search_hidden_dimers(xl_objects):
    af_xl_obj = [obj for obj in xl_objects if obj.error != INVALID_ERROR_VALUE]
    af_xl_obj = [obj for obj in af_xl_obj if obj.distance > 35]
    viol_uniport_dict = dict()
    for obj in af_xl_obj:
        if obj.uniport_a not in viol_uniport_dict:
            viol_uniport_dict[obj.uniport_a] = 0
        viol_uniport_dict[obj.uniport_a] += 1
    violated = sorted(list(viol_uniport_dict.items()), key=lambda item: item[1], reverse=True)
    print(violated)


def analyze_structures(xl_objects):
    inter_pdbs = pdb_files_manager.get_inter_pdbs()
    cif_files = general_utils.load_obj('cif_files')
    af_structure = 0
    pdb_structures = 0
    structures = set()
    for obj in xl_objects:
        pdb_path = obj.pdb_path
        structures.add(pdb_path)
    for s in structures:
        if s[-3:] == 'cif' or s[-3:] == 'ent':
            pdb_structures += 1
        else:
            af_structure += 1
    print(len(structures), af_structure, pdb_structures)


def plot_xl_objects_distance_histogram(xl_objects, title='Crosslink Distance Distribution', xl_type='all', bins_lim=46):
    dist_intra_obj = np.array([obj.distance for obj in xl_objects if obj.xl_type == INTRA_XL])
    dist_inter_obj = np.array([obj.distance for obj in xl_objects if obj.xl_type == INTER_XL])
    dist_all = np.array([obj.distance for obj in xl_objects])
    bins = list(range(bins_lim))
    if xl_type == 'intra':
        general_utils.plot_histogram([dist_intra_obj], title, 'Distance', normalize=False, colors=['cadetblue'], bins=bins)
    elif xl_type == 'inter':
        general_utils.plot_histogram([dist_inter_obj], title, 'Distance', normalize=False, colors=['firebrick'], bins=bins)
    elif xl_type == 'all':
        general_utils.plot_histogram([dist_all], title, 'Distance', normalize=False, colors=['olive'], bins=bins)
    elif xl_type == 'union':
        general_utils.plot_histogram([dist_intra_obj, dist_inter_obj], title, 'Distance', data_labels=['Intra', 'Inter'],
                                     normalize=True, colors=['cadetblue', 'firebrick'], bins=bins)
    elif xl_type == 'by_linker':
        dist_dss = np.array([obj.distance for obj in xl_objects if LINKER_DICT[obj.linker_type] == 3])
        dist_dsso = np.array([obj.distance for obj in xl_objects if LINKER_DICT[obj.linker_type] == 1])
        dist_bdp_nhp = np.array([obj.distance for obj in xl_objects if LINKER_DICT[obj.linker_type] == 2])
        general_utils.plot_histogram([dist_dss, dist_dsso, dist_bdp_nhp], title, 'Distance', data_labels=['DSS', 'DSSO', 'BDP-NHP'],
                                     normalize=True, colors=['palegreen', 'firebrick', 'royalblue'], bins=bins)


def missing_res_single_pdb(pdb_name, suff, uni, missing_by_pdbs):
    pdb = pdb_name.split('/')[-1].split('.')[0]
    if pdb == uni:
        pdb_a = pdb_name
    else:
        pdb_a = f"{pdb_name.split('.')[0]}_{uni}.{suff}"
    if isfile(pdb_a):
        key = pdb_a.split('/')[-1].split('.')[0]
        if key[:3] == 'pdb':
            key = key[3:]
        if key not in missing_by_pdbs:
            res = str(subprocess.check_output(f"{MISSING_RESIDUES_SCRIPT} {pdb_a}", shell=True))
            missing_by_pdbs[key] = int(res.split(':')[2].split()[0])
    else:
        print(f"no such file {pdb_a}")


def count_missing_residues(xl_objects, cif_files=None, inter_pdbs=None, upd_exist=True, exclude_cif=False):
    if cif_files is None:
        cif_files = general_utils.load_obj('cif_files')
    if inter_pdbs is None:
        inter_pdbs = pdb_files_manager.get_inter_pdbs()
    missing_by_pdbs = dict()
    if upd_exist:
        missing_by_pdbs = general_utils.load_obj('missing_residues_by_pdb')
    for i, obj in enumerate(xl_objects):
        if obj.uniport_a == obj.uniport_b and obj.pdb_file not in inter_pdbs:
            continue
        pdb_name = obj.pdb_path
        suff = pdb_name.split('.')[-1]  # cif or pdb
        if exclude_cif and suff == 'cif':
            continue
        missing_res_single_pdb(pdb_name, suff, obj.uniport_a.split('-')[0].split('.')[0], missing_by_pdbs)
        if obj.uniport_a != obj.uniport_b:
            missing_res_single_pdb(pdb_name, suff, obj.uniport_b.split('-')[0].split('.')[0], missing_by_pdbs)
        if i % 1000 == 0:
            print(f"processed missing res for {i} objects")
    general_utils.save_obj(missing_by_pdbs, 'missing_residues_by_pdb')
    print("finish count_missing_residues")
    return missing_by_pdbs


def analyze_inter_xl(exclude_cif=False, objects_file_name=None):
    cif_files = general_utils.load_obj('cif_files')
    xl_obj = get_upd_and_filter_processed_objects_pipeline(objects_file_name)
    print("Counting missing residues: ")
    inter_pdbs = pdb_files_manager.get_inter_pdbs()
    missing_res_by_pdbs = count_missing_residues(xl_obj, cif_files, upd_exist=True, exclude_cif=exclude_cif)
    pdb_dict = dict()
    triple_dict = dict()
    pdb_chains = pdb_files_manager.num_of_chains_in_pdb(xl_obj, True)
    # pdb_chains = pdb_files_manager.get_num_of_chains_in_pdb()
    for obj in xl_obj:
        if (exclude_cif and obj.pdb_file in cif_files) or (obj.uniport_a == obj.uniport_b and obj.pdb_file not in inter_pdbs):
            continue
        if obj.pdb_file not in pdb_dict:
            pdb_dict[obj.pdb_file] = 0
        triple_a = (obj.chain_a, obj.chain_b, obj.pdb_file, obj.linker_type)
        triple_b = (obj.chain_b, obj.chain_a, obj.pdb_file, obj.linker_type)
        if triple_a in triple_dict:
            triple_dict[triple_a] += 1
        elif triple_b in triple_dict:
            triple_dict[triple_b] += 1
        else:
            triple_dict[triple_a] = 1
        pdb_dict[obj.pdb_file] += 1
    print(sorted(pdb_dict.items(), key=lambda item: item[1], reverse=True))
    print(sorted(triple_dict.items(), key=lambda item: item[1], reverse=True))
    different_chains_dict = {k: v for k, v in triple_dict.items() if k[0] != k[1]}
    print(sorted(different_chains_dict.items(), key=lambda item: item[1], reverse=True))
    different_chains_pdb_dict = dict()
    manual_hist = dict()
    for k, v in different_chains_dict.items():
        if v not in manual_hist:
            manual_hist[v] = 0
        manual_hist[v] += 1
        missing_res = missing_res_by_pdbs[f"{k[2].lower()}_{k[0]}"] + missing_res_by_pdbs[f"{k[2].lower()}_{k[1]}"]
        if v > 2 and pdb_chains[k[2]] <= 4:
            key = (k[0], k[1], k[2], pdb_chains[k[2]], missing_res, k[3])
            different_chains_pdb_dict[key] = v
    print()
    print("Chain A, Chain B, PDB, # chains, # missing residues: # cross links, linker type")
    different_chains_pdb_dict = sorted(sorted(different_chains_pdb_dict.items(), key=lambda item: item[0][3]), key=lambda item: item[1], reverse=True)
    print(different_chains_pdb_dict)
    print()
    print("Chain A, Chain B, PDB, # chains, # missing residues: # cross links, linker type")
    print([item for item in different_chains_pdb_dict if item[0][4] < 20])
    manual_hist = sorted(manual_hist.items(), key=lambda item: item[0])
    x = [manual_hist[i][0] for i in range(len(manual_hist))]
    y = [manual_hist[i][1] for i in range(len(manual_hist))]
    general_utils.plot_line(x, y, "# Inter-crosslinks between single pair", "# crosslinks", "# pairs")
    plot_xl_objects_distance_histogram(xl_obj)


def filter_pdbs_from_xl_objects(xl_objects):
    filter_pairs = {'3WYO': ('A', 'B')}
    new_xl_objects = []
    for obj in xl_objects:
        if obj.pdb_file in filter_pairs and \
                ((filter_pairs[obj.pdb_file][0] == obj.chain_a and filter_pairs[obj.pdb_file][1] == obj.chain_b)
                 or (filter_pairs[obj.pdb_file][0] == obj.chain_b and filter_pairs[obj.pdb_file][1] == obj.chain_a)):
            continue
        new_xl_objects.append(obj)
    return new_xl_objects


def xlink_db_data_update_pipeline(data_path='/cs/labs/dina/seanco/xl_parser/data/xl_db_data_update.txt',
                                  stop_for_prediction=True, after_prediction=False):
    if not after_prediction:
        # data_table_parser.read_all_clear_dup(data_path)
        # general_xl_parser.XlParser.convert_old_xl_list_to_cross_link_objects()
        data_table_parser.update_uniport_fasta_dict_from_xl_objects()
        data_table_parser.fix_fasta_int_values()

        for i in range(2):
            xl_objects = CrossLink.load_all_xl_objects()
            pdb_files_manager.update_xl_objects_with_obsolete_uniports(xl_objects)
            alpha_fold_files.get_af_of_cross_link_objects(xl_objects)
            pdb_files_manager.download_pdbs_from_xl_objects(xl_objects)
            pdb_files_manager.create_cif_set()
            pdb_files_manager.missing_pdb_sequences_to_fasta()
        # After this, run alphafold on 'fasta_to_predict.fasta'

    if not stop_for_prediction:
        process_objects_pipeline()


def process_objects_pipeline(out_name='new_all_objects', xl_objects=None):
    if xl_objects is None:
        xl_objects = CrossLink.load_all_xl_objects()
    # xl_objects = general_utils.load_obj('error_multichain_objects')
    pdb_files_manager.update_xl_objects_with_obsolete_uniports(xl_objects)
    CrossLink.fix_peptides_all_xl_objects(xl_objects)
    xl_objects, uni_dict = CrossLink.map_inter_xl_objects_to_pdb(xl_objects)
    pdb_files_manager.create_dimers_dict(xl_objects, uni_dict)
    inter_pdbs = pdb_files_manager.get_inter_pdbs()
    download_pdb_cif_af_files(xl_objects, inter_pdbs)
    CrossLink.process_all_xl(xl_objects, uni_chain_dict=uni_dict, out_name=out_name, inter_pdbs=inter_pdbs)
    # CrossLink.process_all_xl(xl_objects, uni_chain_dict=None, out_name='intra_af_objects')
    # CrossLink.process_all_xl(xl_objects)


def get_upd_and_filter_processed_objects_pipeline(specific_files=None):
    if specific_files is None:
        processed_xl_objects = CrossLink.load_all_xl_objects(PROCESSED_OBJ_DIR, PROCESSED_OBJ_DIR_PREFIX)
    else:
        processed_xl_objects = []
        for file in specific_files:
            processed_xl_objects += general_utils.load_obj(file, PROCESSED_OBJ_DIR)
    print(f"initial amount: {len(processed_xl_objects)}")
    processed_xl_objects = CrossLink.update_linker_by_origin_db(processed_xl_objects)
    processed_xl_objects = CrossLink.filter_xl_obj(processed_xl_objects, 45, 100)
    processed_xl_objects = CrossLink.clear_duplicates_of_xl_objects(processed_xl_objects, True)
    inter_pdbs = pdb_files_manager.get_inter_pdbs()
    cif_files = general_utils.load_obj('cif_files')
    processed_xl_objects = CrossLink.clear_duplicates_of_dimers(processed_xl_objects, cif_files, inter_pdbs)
    pdb_files_manager.fix_objects_pdb_path(processed_xl_objects, cif_files, inter_pdbs)
    # processed_xl_objects = pdb_files_manager.filter_objects_from_list_by_pdb(processed_xl_objects, inter_pdbs=inter_pdbs)
    return processed_xl_objects


def download_pdb_cif_af_files(xl_objects, inter_pdbs=None):
    if inter_pdbs is None:
        inter_pdbs = pdb_files_manager.get_inter_pdbs()
    pdb_files_manager.download_pdbs_from_xl_objects(xl_objects, inter_pdbs)
    pdb_files_manager.create_cif_set()
    alpha_fold_files.get_af_of_cross_link_objects(xl_objects)


def extract_features_pipeline(processed_xl_objects=None, objects_file_name=None):
    """
    Extracting features from xk objects - preparation for creating data set for NN. Takes few hours to run.
    :param processed_xl_objects: Objects list or None. if None - name of pickle file must be provided.
    :param objects_file_name: pickle file name of xl objects list.
    :return:
    """
    if processed_xl_objects is None:
        processed_xl_objects = general_utils.load_obj(objects_file_name)
    inter_pdbs = pdb_files_manager.get_inter_pdbs()
    CrossLink.export_objects_to_txt_format(processed_xl_objects, inter_pdbs=inter_pdbs,
                                           out_path=pdb_files_manager.INTER_AND_INTRA_LYS_XL_FILES_PATH, skip_exist=False, dimer_form=False)
    pdb_files_manager.extract_xl_features_from_xl_objects(processed_xl_objects, inter_pdbs=inter_pdbs,
                                                          xl_dir=pdb_files_manager.INTER_AND_INTRA_LYS_XL_FILES_PATH,
                                                          out_dir=pdb_files_manager.INTER_INTRA_LYS_XL_NEIGHBORS_FILES_PATH,
                                                          skip_exist=False)
    pdb_files_manager.read_features_to_dict(False, pdb_files_manager.INTER_INTRA_LYS_XL_NEIGHBORS_FILES_PATH,
                                            pdb_files_manager.XL_NEIGHBORS_FEATURE_DICT_INTER_INTRA_LYS)


def create_objects_from_csv(import_from_csv=True, parse_cfg=None, is_processed=True, out_name='new_all_objects',
                            input_file_path=None, final_save_name=None, download_files=False):
    """
    Creates xl object list from a file.
    :param import_from_csv: if True, creates xl objects from a file in our format
    (identical to the objects exporting format). Else - parse objects from different format.
    :param parse_cfg: Config file with parsing instructions. (only if import_from_csv is False).
    :param is_processed: If False, objects does not have distance data - extracting distances from pdb (very slow)
    :param out_name: If is_processed is False, name for savinf the processed objects
    :param input_file_path: For importing from csv
    :param final_save_name: if not None, saves the object list with the given name (give only name! not path)
    :param download_files: if true - download pdb and af files to directories defined in 'pdb_files_manager.py' and 'alpha_fold_files.py'
    :return:
    """
    if import_from_csv:
        xl_objects = general_xl_parser.XlParser.import_xl_objects_from_csv(input_file_path,
                                                                           pdb_files_manager.PDB_FILES_DIR)
    else:
        xl_objects = general_xl_parser.run_parsing(parse_cfg)
    if not is_processed:
        process_objects_pipeline(out_name, xl_objects)
        xl_objects = get_upd_and_filter_processed_objects_pipeline()
    if download_files and is_processed:
        download_pdb_cif_af_files(xl_objects)
    if final_save_name is not None:
        general_utils.save_obj(xl_objects, final_save_name)
    return xl_objects


class ArgParser(argparse.ArgumentParser):

    def __init__(self, **kwargs):
        super(ArgParser, self).__init__(**kwargs)
        self.add_argument(
            "--run_mode",
            help="Which utility to run: import / parse / extract",
            default="import",
            type=str,
        )
        self.add_argument(
            "--parse_cfg",
            help="config file for parsing instructions (optional)",
            default=None,
            type=str,
        )
        self.add_argument(
            "--is_processed",
            help="Is objects already contain distance info",
            default=True,
            type=bool
        )
        self.add_argument(
            "--out_name",
            help="Name for saving objects if processing is needed (optional)",
            default=None,
            type=str
        )
        self.add_argument(
            "--input_file_path",
            help="Path of file for importing objects",
            default=None,
            type=str
        )
        self.add_argument(
            "--final_save_name",
            help="Name of file for saving objects",
            default=None,
            type=str
        )
        self.add_argument(
            "--download_files",
            help="If true, downloading pdb, cif and alpha_fold pdb files for the objects, necessary for extracting "
                 "distances, features and creating dataset",
            default=True,
            type=bool
        )

    def parse_args(self, args=None, namespace=None):
        """ Parse the input arguments """
        args = super(ArgParser, self).parse_args(args, namespace)
        return args


def main():
    p = ArgParser()
    args = p.parse_args()
    if args.run_mode == 'import':
        create_objects_from_csv(True, args.parse_cfg, args.is_processed, args.out_name, args.input_file_path,
                                args.final_save_name, args.download_files)
    elif args.run_mode == 'parse':
        create_objects_from_csv(False, args.parse_cfg, args.is_processed, args.out_name, args.input_file_path,
                                args.final_save_name, args.download_files)
    elif args.run_mode == 'extract':
        extract_features_pipeline(None, args.final_save_name)
    else:
        raise argparse.ArgumentTypeError(f"Invalid run mode: {args.run_mode}")


if __name__ == "__main__":
    main()
