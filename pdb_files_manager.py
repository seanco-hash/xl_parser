import subprocess
import multiprocessing
import os
import general_utils
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
import torch
from Bio.SeqUtils import seq1
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from os.path import isfile
import shutil
import alpha_fold_files
import numpy as np
import copy
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.firefox import GeckoDriverManager
from selenium.webdriver.firefox.options import Options as FirefoxOptions
from selenium.webdriver.common.by import By
import json

SCWRL_OUT_DIR = '/cs/labs/dina/seanco/xl_parser/scwrl/'
SCWRL_SCRIPT = '~dina/software/progs/scwrl4/Scwrl4 '
DSSP_CIF_SCRIPT = '/cs/labs/dina/seanco/needle/project/venv_needle/bin/biolib run bio_utils/DSSP -i '
DSSP_SCRIPT = '/cs/labs/dina/seanco/xl_mlp_nn/run_dssp.sh '
WORKDIR = "/cs/labs/dina/seanco/xl_parser/"
PDB_FILES_DIR = "/cs/labs/dina/seanco/xl_parser/pdbs/"
FASTA_TO_PREDICT = "/cs/labs/dina/seanco/xl_parser/fasta_to_predict.fasta"
OUTPUT_XL_FILES_PATH = "/sci/labs/dina/seanco/xl_neighbors/xl_files/"
UNPROCESSED_XL_FILES_PATH = "/cs/labs/dina/seanco/xl_neighbors/unprocessed_xl_files/"
INTRA_LYS_XL_FILES_PATH = "/cs/labs/dina/seanco/xl_neighbors/intra_lys_xl_files/"
UNFILTERED_XL_FILES_PATH = "/cs/labs/dina/seanco/xl_neighbors/unfiltered_xl_files/"
INTER_AND_INTRA_LYS_XL_FILES_PATH = "/cs/labs/dina/seanco/xl_neighbors/inter_intra_lys_xl_files/"
XL_NEIGHBORS_EXE = "/cs/labs/dina/seanco/xl_neighbors/xl_neighbors"
XL_NEIGHBORS_FILES_PATH = "/sci/labs/dina/seanco/xl_neighbors/feature_files/"
INTRA_LYS_XL_NEIGHBORS_FILES_PATH = "/cs/labs/dina/seanco/xl_neighbors/intra_lys_feature_files/"
INTER_INTRA_LYS_XL_NEIGHBORS_FILES_PATH = "/cs/labs/dina/seanco/xl_neighbors/inter_intra_lys_feature_files/"
XL_NEIGHBORS_FILES_NO_CIF = "/cs/labs/dina/seanco/xl_neighbors/feature_files_no_cif/"
XL_NEIGHBORS_FEATURE_DICT = "xl_neighbors_feature_dict"
XL_NEIGHBORS_FEATURE_DICT_INTRA_LYS = "xl_neighbors_feature_dict_intra_lys"
XL_NEIGHBORS_FEATURE_DICT_INTER_INTRA_LYS = "xl_neighbors_feature_dict_inter_intra_lys"
XL_NEIGHBORS_FEATURE_DICT_NO_CIF = "xl_neighbors_feature_dict_no_cif"
XL_NEIGHBORS_FEATURE_DICT_ONLY_CIF = "xl_neighbors_feature_dict_only_cif"
UNIPARC_WEB_PATH = 'https://www.uniprot.org/uniparc?query='
UNIPROT_WEB_PATH = 'https://www.uniprot.org/uniprotkb?query='
UNIPORTS_DICT = 'obsolete_new_uniports_dict'
PARSE_HEADER_COMMAND = "grep COMPND pdb5a2q.ent | grep -E 'MOL_ID|CHAIN'"
INVALID_OFFSET = -11111
LINKER_MAX_DIST = {'DSSO': 32, 'BDP_NHP': 35, 'DSS': 32, 'UNKNOWN': 35, 'BDP-NHP': 35, '': 35, 'LEIKER': 35, None: 35,
                   'QPIR': 35, 'A-DSBSO': 35}
EXTENDED_AA_LETTERS = "ACDEFGHIKLMNPQRSTVWYBXZJUO"
AA_TABLE_IDX = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10,
                'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19,
                'B': 20, 'X': 21, 'Z': 22, 'J': 23, 'U': 24, 'O': 25}
AA_LETTERS = "ACDEFGHIKLMNPQRSTVWY"
PDB_2_FASTA = '/cs/staff/dina/utils/pdb2fasta'
FASTA_PATH = '/cs/labs/dina/seanco/xl_parser/fasta_files/'


def copy_wrong_dir_models():
    to_copy = []
    for file_name in os.listdir(WORKDIR + 'best_models'):
        parts = file_name.split('_')
        if len(parts) > 4:
            to_copy.append(file_name)

    for file_name in to_copy:
        shutil.copy(WORKDIR + 'best_models/' + file_name, WORKDIR + '/best_models/inter')
        os.remove(WORKDIR + 'best_models/' + file_name)

    print('finish')


def download_single_pdb(name, downloaded, invalid, parser, file_format):
    """
    Downloads single pdb file with Biopython PDB Parser
    """
    if name not in downloaded:
        if name not in invalid:
            ret_val = parser.retrieve_pdb_file(name, pdir=PDB_FILES_DIR, overwrite=False, file_format=file_format)
            if os.path.exists(ret_val):
                downloaded.add(name)
                return True
            else:
                # invalid.add(name)
                return False
    return True


def download_missing_pdbs():
    pdb1 = PDB.PDBList()
    wanted_format = "mmCif"
    missing_pdbs = general_utils.load_obj("missing_pdbs")
    already_downloaded = set()
    invalid_download = set()
    for pdb in missing_pdbs:
        download_single_pdb(pdb, already_downloaded, invalid_download, pdb1, wanted_format)

    print("Downloaded: " + str(len(already_downloaded)))
    print("Invalid: " + str(len(invalid_download)))


def download_pdbs_from_xl_objects(xl_objects, inter_pdbs=None):
    pdb1 = PDB.PDBList()
    cif_files = general_utils.load_obj('cif_files')
    already_downloaded = set()
    invalid_download = set()
    for obj in xl_objects:
        pdb_file_name = find_pdb_file_from_xl_obj(obj, cif_files, inter_pdbs=inter_pdbs)
        if (pdb_file_name is None or not isfile(pdb_file_name)) and obj.pdb_file != '\n' and \
                obj.pdb_file != '' and obj.pdb_file != 'None':
            obj.pdb_file = obj.pdb_file.split('\n')[0]
            res = download_single_pdb(obj.pdb_file, already_downloaded, invalid_download, pdb1, 'pdb')
            if not res:
                res = download_single_pdb(obj.pdb_file, already_downloaded, invalid_download, pdb1, 'mmCif')
                if not res:
                    invalid_download.add(obj.pdb_file)
    print("Downloaded with pdb format: " + str(len(already_downloaded)))
    print("Invalid: " + str(len(invalid_download)))
    print(invalid_download)


def get_pdb_xl_known_structures(xl_list):
    pdb1 = PDB.PDBList()
    already_downloaded = set()
    invalid_download = set()
    for sample in xl_list:
        if sample[general_utils.XL_UNIPORT_A] != sample[general_utils.XL_UNIPORT_B]:
            if sample[general_utils.XL_PDB_TYPE_A] == general_utils.XL_AVAILABLE_IN_PDB:
                name = sample[general_utils.XL_PDB_A]
                res = download_single_pdb(name, already_downloaded, invalid_download, pdb1, 'pdb')
                if not res:
                    res = download_single_pdb(name, already_downloaded, invalid_download, pdb1, 'mmCif')
                    if not res:
                        invalid_download.add(name)

    print("Downloaded with pdb format: " + str(len(already_downloaded)))
    print("Invalid: " + str(len(invalid_download)))


def get_pdb_files(xl_list, wanted_format="pdb"):
    pdb1 = PDB.PDBList()
    already_downloaded = set()
    invalid_download = set()
    for sample in xl_list:
        pdb_names = sample[-1].split(',')
        name = pdb_names[0]
        if not isfile(PDB_FILES_DIR + 'pdb' + name.lower + '.ent'):
            download_single_pdb(name, already_downloaded, invalid_download, pdb1, wanted_format)

    print("Downloaded: " + str(len(already_downloaded)))
    print("Invalid: " + str(len(invalid_download)))


def check_missing_files(xl_list, wanted_format="pdb"):
    pdb1 = PDB.PDBList()
    already_downloaded = set()
    invalid_download = set()
    invalid_first_pdb = set()
    samples_without_pdb = []
    for sample in xl_list:
        at_least_one_succ = False
        pdb_names = sample[-1].split(',')
        i = 0
        while i < len(pdb_names):
            name = pdb_names[i]
            if name not in already_downloaded:
                if name not in invalid_download:
                    ret_val = pdb1.retrieve_pdb_file(name, pdir="./pdbs", overwrite=False,
                                                     file_format=wanted_format)
                    if os.path.exists(ret_val):
                        already_downloaded.add(name)
                        at_least_one_succ = True
                        i = len(pdb_names)
                    else:
                        invalid_download.add(name)
                        if i == 0:
                            invalid_first_pdb.add(name)
                else:
                    if i == 0:
                        invalid_first_pdb.add(name)
            i += 1
        if not at_least_one_succ:
            samples_without_pdb.append(sample)

    print(len(invalid_first_pdb))
    print(len(samples_without_pdb))
    print(len(invalid_download))
    general_utils.save_obj(invalid_first_pdb, "best_missing_pdbs")
    general_utils.save_obj(samples_without_pdb, "samples_wo_pdbs")
    general_utils.save_obj(invalid_download, "missing_pdbs")


def create_cif_set(pdb_files_dir='./pdbs'):
    cif_files = general_utils.load_obj('cif_files')
    print(f"Initial set size: {len(cif_files)}")
    for filename in os.listdir(pdb_files_dir):
        if filename.endswith(".cif"):
            pdb_name = filename.split('.')[0]
            cif_files.add(pdb_name.upper())

    general_utils.save_obj(cif_files, 'cif_files')
    print(f"End set size: {len(cif_files)}")


def create_af_set():
    af_set = set()
    for filename in os.listdir('./pdbs/alpha_fold/pdb_files'):
        unicode = filename.split('.')[0]
        af_set.add(unicode)

    general_utils.save_obj(af_set, 'af_pdbs_unicodes')
    print(len(af_set))


def generate_model_name_low(pdb, unicode_a, unicode_b=None, models_dir='best_models/', pref='pdb'):
    pdb_name = pref + pdb.lower()
    name = pdb_name + '_' + unicode_a
    if unicode_b is not None:
        name += '_' + unicode_b
    name += '.B99990001.pdb'
    return models_dir + name


def find_pdb_file_from_xl_obj(xl_obj, cif_files, inter_pdbs=None):
    pdb = xl_obj.pdb_file
    file_name = None
    if xl_obj.uniport_a == xl_obj.uniport_b and \
            (inter_pdbs is None or xl_obj.pdb_file not in inter_pdbs):
        file_name = alpha_fold_files.generate_af_model_name(xl_obj.uniport_a)
    if file_name is None or not isfile(file_name):
        if pdb != '' and pdb != 'None' and pdb is not None:
            pref = 'pdb'
            suff = '.ent'
            if pdb in cif_files:
                pref = ''
                suff = '.cif'
            file_name = PDB_FILES_DIR + pref + pdb.lower() + suff
        else:
            return None
    if not isfile(file_name):
        return None
    return file_name


def update_xl_objects_with_obsolete_uniports(xl_objects):
    uni_dict = general_utils.load_obj(UNIPORTS_DICT)
    for obj in xl_objects:
        obj.pdb_file = obj.pdb_file.split('\n')[0]
        uni_a = obj.uniport_a.split('.')[0]
        if uni_a in uni_dict and uni_a[0] != 'M':
            if len(uni_dict[uni_a].split(" ")) == 1:
                obj.uniport_a = uni_dict[uni_a]
        uni_b = obj.uniport_b.split('.')[0]
        if uni_b in uni_dict and uni_b[0] != 'M':
            if len(uni_dict[uni_b].split(" ")) == 1:
                obj.uniport_b = uni_dict[uni_b]


def get_web_browser():
    options = FirefoxOptions()
    options.add_argument('--headless')
    s = Service(GeckoDriverManager().install())
    browser = webdriver.Firefox(service=s, options=options)
    return browser


def find_uniports_in_web(base_address, uniport_item_xpath, uniports, new_uniports, browser, first=True):
    for uni in uniports:
        if uni not in new_uniports:
            try:
                address = base_address + uni
                browser.get(address)
                if first:
                    table_click = browser.find_element(By.XPATH, '/html/body/form/div/span/label[2]/input')
                    table_click.click()
                    show_res_button = browser.find_element(By.XPATH, '/html/body/form/div/section/button')
                    show_res_button.click()
                    first = False
                entry = browser.find_element(By.XPATH, uniport_item_xpath)
                new_uniport = entry.text
                new_uniports[uni] = new_uniport
            except Exception as e:
                print(f"Exception with uniport {uni}")


def update_alternative_uniports_from_manual_file(new_uniports=None):
    if new_uniports is None:
        new_uniports = general_utils.load_obj(UNIPORTS_DICT)
    with open('/cs/labs/dina/seanco/xl_parser/data/alternative_uniports.txt', 'r') as f:
        for line in f:
            old, new = line.split()
            new_uniports[old] = new.split('\n')[0]
    return new_uniports


def find_alternative_uniports(uniports):
    split_uniports = [uni.split('.') for uni in uniports]
    new_uniports = dict()
    if isfile(general_utils.OBJ_DIR + UNIPORTS_DICT + '.pkl'):
        new_uniports = general_utils.load_obj(UNIPORTS_DICT)
    print(f"len of new uniports dict before: {len(new_uniports)}")
    obsolet_uniprot_xpath = '//*[@id="root"]/div/div[1]/div/section[2]/div[3]/table/tbody/tr/td[5]/ul/li[1]/a'
    dna_uniports_xpath = '//*[@id="root"]/div/div[1]/div/section[2]/div[3]/table/tbody/tr/td[2]/span/a'
    browser = get_web_browser()
    obsolete_uniports = [uni[0] for uni in split_uniports if len(uni) == 1]
    find_uniports_in_web(UNIPARC_WEB_PATH, obsolet_uniprot_xpath, obsolete_uniports, new_uniports, browser, True)
    dna_uniports = [uni[0] for uni in split_uniports if len(uni) > 1]
    find_uniports_in_web(UNIPROT_WEB_PATH, dna_uniports_xpath, dna_uniports, new_uniports, browser, False)
    update_alternative_uniports_from_manual_file(new_uniports)
    print(new_uniports)
    print(f"len of new uniports dict after: {len(new_uniports)}")

    general_utils.save_obj(new_uniports, UNIPORTS_DICT)
    return new_uniports


def missing_pdb_sequences_to_fasta():
    fasta_dict = general_utils.load_obj('fasta_files_dict_by_uniport')
    missing_fasta = 0
    large_seq = 0
    add_to_fasta = {}
    obsolete = []
    with open('/cs/labs/dina/seanco/xl_parser/uniports_to_predict.txt', 'r') as f:
        uniports = [line[:-1] for line in f]
        for uniport in uniports:
            uniport = uniport.split('-')[0].split('.')[0]
            if uniport[0] == 'M':
                uniport = uniport[0] + 'COT' + uniport[1:]
            if uniport in fasta_dict:
                seq = fasta_dict[uniport]
                if 16 < len(seq) < 2700:
                    existing_model_name = alpha_fold_files.generate_af_model_name(uniport)
                    if not isfile(existing_model_name):
                        if 'U' in seq:
                            seq = seq.replace("U", "C")
                            print("replaced U with C")
                        add_to_fasta[uniport] = seq
                else:
                    if len(seq) == 0:
                        obsolete.append(uniport)
                    print(uniport + " length: ", len(seq))
                    large_seq += 1
            else:
                if '.' in uniport:
                    obsolete.append(uniport)
                print("missing fasta: ", uniport)
                missing_fasta += 1
        general_utils.seqs_dict_to_fasta(add_to_fasta, FASTA_TO_PREDICT)
        print("Missing fasta files: ", missing_fasta)
        print(f"large sequences: {large_seq}")
        print(f"number of empty fasta sequences: {len(obsolete)}")
        with open('/cs/labs/dina/seanco/xl_parser/obsolete_uniports.txt', 'w') as o_f:
            find_alternative_uniports(obsolete)
            for uni in obsolete:
                o_f.write(uni + '\n')


def check_pdb_chain_file_exist(chain_id, pdb_file_path):
    pdb_name_split = pdb_file_path.split('/')
    pref, suff = pdb_name_split[-1].split('.')
    new_name = pref + '_' + chain_id + '.' + suff
    new_path = '/'.join(pdb_name_split[:-1]) + '/' + new_name
    if isfile(new_path):
        return True, new_path
    return False, new_path


def export_single_chain_to_pdb(chain, pdb_file_path):
    exist, new_path = check_pdb_chain_file_exist(chain.id, pdb_file_path)
    if exist:
        return
    if new_path[-3:] == 'cif':
        io = MMCIFIO()
    else:
        io = PDBIO()
    io.set_structure(chain)
    io.save(new_path)


def export_chains_to_pdb(chains, pdb_file_path, chain_a, chain_b):
    done = 0
    for chain in chains:
        if chain_a == chain.id or chain_b == chain.id:
            export_single_chain_to_pdb(chain, pdb_file_path)
            done += 1
        if done == 2:
            break


def export_chains_to_pdb_from_xl_objects(xl_objects):
    inter_pdbs = get_inter_pdbs()
    cif_files = general_utils.load_obj('cif_files')
    pdb_parser = PDBParser(PERMISSIVE=1)
    cif_parser = MMCIFParser()
    xl_objects = sorted(xl_objects, key=lambda o: o.uniport_a)
    prev_pdb = '-1'
    for obj in xl_objects:
        pdb_file_path = find_pdb_file_from_xl_obj(obj, cif_files, inter_pdbs=inter_pdbs)
        a_exist, _ = check_pdb_chain_file_exist(obj.chain_a, pdb_file_path)
        b_exist, _ = check_pdb_chain_file_exist(obj.chain_b, pdb_file_path)
        if (a_exist and b_exist) or pdb_file_path is None:
            continue
        if pdb_file_path != prev_pdb and pdb_file_path is not None:
            if pdb_file_path[-3:] == 'cif':
                structure = cif_parser.get_structure(pdb_file_path.split('/')[-1], pdb_file_path)
            else:
                structure = pdb_parser.get_structure(pdb_file_path.split('/')[-1], pdb_file_path)
            tmp_chains = list(structure.get_chains())
            prev_pdb = pdb_file_path
        export_chains_to_pdb(tmp_chains, pdb_file_path, obj.chain_a, obj.chain_b)


def find_chain_ids_cif(obj, mmcif_dict, optional_chain_a, optional_chain_b):
    for i, uni in enumerate(mmcif_dict['_struct_ref_seq.pdbx_db_accession']):
        if uni == obj.uniport_a:
            optional_chain_a.add(mmcif_dict['_struct_ref_seq.pdbx_strand_id'][i])
        if uni == obj.uniport_b:
            optional_chain_b.add(mmcif_dict['_struct_ref_seq.pdbx_strand_id'][i])


def find_chain_ids_pdb(obj, entity_chain_dict, optional_chain_a, optional_chain_b, uni_chain_dict, pdb_name):
    if pdb_name.split('.')[-1] == 'cif':
        return
    id_list = uni_chain_dict[(obj.uniport_a, obj.uniport_b)]
    for ids in id_list:
        id1, id2 = ids
        for c in entity_chain_dict[pdb_name][str(id1)]:
            optional_chain_a.add(c)
        for c in entity_chain_dict[pdb_name][str(id2)]:
            optional_chain_b.add(c)


def process_multichain_xl_object(obj, chains, fasta_dict, structure, max_dist=45, counters=None, uni_chain_dict=None,
                                 pdb_file_path=None, mmcif_dict=None, entity_chain_dict=None, error_object=None, already_processed=False):
    if counters is None:
        counters = np.zeros(5)
    if (obj.uniport_a not in fasta_dict or obj.uniport_b not in fasta_dict) and uni_chain_dict is None:
        counters[2] += 1
        return None
    tmp_chains = list(structure.get_chains())
    optional_chain_a, optional_chain_b = set(), set()
    if uni_chain_dict is None:
        print("error: process multichain object by uni_dict is none")
        return None
    elif already_processed:
        optional_chain_a = {obj.chain_a}
        optional_chain_b = {obj.chain_b}
    elif mmcif_dict is not None:
        find_chain_ids_cif(obj, mmcif_dict, optional_chain_a, optional_chain_b)
    else:
        find_chain_ids_pdb(obj, entity_chain_dict, optional_chain_a, optional_chain_b, uni_chain_dict,
                           pdb_file_path.split('/')[-1])

    if len(optional_chain_a) == 0 or len(optional_chain_b) == 0:
        counters[3] += 1
        return None
    tmp_objects = []
    for a in optional_chain_a:
        for b in optional_chain_b:
            tmp_obj = copy.deepcopy(obj)
            tmp_obj.chain_a = a
            tmp_obj.chain_b = b
            tmp_obj.process_single_xl(None, multichain=True, chains=tmp_chains, error_objects=error_object)
            if 0 < tmp_obj.distance:
                tmp_objects.append(tmp_obj)
    if len(tmp_objects) == 1:
        counters[0] += 1
    elif len(tmp_objects) == 0:
        counters[1] += 1
        return None
    else:
        tmp_objects.sort(key=lambda o: o.distance)
        counters[4] += 1
    tmp = tmp_objects[0]
    try:
        export_chains_to_pdb(tmp_chains, pdb_file_path, tmp.chain_a, tmp.chain_b)
    except Exception as e:
        print(e)
        print(f"export chain error. file: {pdb_file_path}, chains: {tmp.chain_a}, {tmp.chain_b}, uniports: {obj.uniport_a}, {obj.uniport_b}")
        return None
    return tmp


def create_single_pdb_xl_list_txt_file(xl_objects, pdb_name, cif_files, counters, max_dist=45, inter_pdbs=None,
                                       out_path=OUTPUT_XL_FILES_PATH, skip_exist=True, exist=[], dimer_dict=None):
    """
    Creates .txt file for single with all the cross links related to it
    in the format 'chain_number residue_number chain_number residue_number'
    """
    existing_lines = set()
    new_file_path = out_path + pdb_name + ".txt"
    if isfile(new_file_path) and skip_exist:
        exist[0] += 1
        f = open(new_file_path, 'a')
    else:
        f = open(new_file_path, 'w')
        # pdb_file = None
    for obj in xl_objects:
        chain_a, chain_b = obj.chain_a, obj.chain_b
        if chain_a == '' or chain_b == '':
            if obj.uniport_a == obj.uniport_b and (inter_pdbs is None or obj.pdb_file not in inter_pdbs or (obj.pdb_path != '' and obj.pdb_path[-3:] == 'pdb')):
                chain_a = chain_b = 'A'
            else:
                chain_a, chain_b = obj.uniport_a, obj.uniport_b
        if dimer_dict is not None:
            if pdb_name in dimer_dict:
                cur_dict = dimer_dict[pdb_name]
                if chain_a in cur_dict:
                    chain_a = cur_dict[chain_a]
                if chain_b in cur_dict:
                    chain_b = cur_dict[chain_b]
        line = f"{obj.res_num_a} {chain_a} {obj.res_num_b} {chain_b} 1 {LINKER_MAX_DIST[obj.linker_type]}\n"
        if line not in existing_lines:
            f.write(line)
            existing_lines.add(line)
            opposite_line = f"{obj.res_num_b} {chain_b} {obj.res_num_a} {chain_a} 1 {LINKER_MAX_DIST[obj.linker_type]}\n"
            existing_lines.add(opposite_line)
    f.close()


def single_thread_extract_xl_features(xl_file_paths, pdb_files, output_path=XL_NEIGHBORS_FILES_PATH,
                                      skip_exist=False, predict=False):
    processes = []
    for i in range(len(xl_file_paths)):
        out_file = output_path + xl_file_paths[i].split('/')[-1]
        if skip_exist and isfile(out_file):
            continue
        if isinstance(pdb_files[i], str):
            p = subprocess.Popen([XL_NEIGHBORS_EXE, xl_file_paths[i], pdb_files[i], output_path])
        else:
            p = subprocess.Popen([XL_NEIGHBORS_EXE, xl_file_paths[i]] + pdb_files[i] + [output_path])
        if predict:
            p.communicate()
        else:
            processes.append(p)
        # if i % 1000 == 0:
        #     print(f"Extracted {i} feature files")
    # exit_codes = [p.wait() for p in processes]


def complete_missing_xl_features(xl_objects, cif_files=None):
    if cif_files is None:
        cif_files = general_utils.load_obj('cif_files')
    missing_objects = []
    for obj in xl_objects:
        pdb_name = obj.pdb_path.split('/')[-1].split('.')[0]
        feat_file = XL_NEIGHBORS_FILES_PATH + pdb_name + '.txt'
        if not isfile(feat_file):
            missing_objects.append(obj)
    print(f"objects with missing feature files: {len(missing_objects)}")
    extract_xl_features_from_xl_objects(missing_objects, cif_files)


def fix_objects_pdb_path(xl_objects, cif_files, inter_pdbs):
    objects = [o for o in xl_objects if o.chain_a != '' and o.pdb_path[-3:] == 'pdb']
    for obj in objects:
        obj.chain_a = obj.chain_b = ''


def extract_xl_features_from_xl_objects(xl_objects, cif_files=None, inter_pdbs=None, xl_dir=OUTPUT_XL_FILES_PATH,
                                        out_dir=None, skip_exist=False, predict=False):
    if cif_files is None:
        cif_files = general_utils.load_obj('cif_files')
    xl_file_paths = set()
    pdb_files_dict = dict()
    for obj in xl_objects:
        pdb_name, pdb_file = get_obj_files_key(obj, cif_files, inter_pdbs, by_chain=False)
        xl_file_path = xl_dir + pdb_name + '.txt'
        xl_file_paths.add(xl_file_path)
        if obj.chain_a != '' and pdb_file[-3:] != 'pdb':
            if pdb_file not in pdb_files_dict:
                pdb_files_dict[pdb_file] = set()
            pdb_pref, suff = pdb_file.split('.')
            pdb_files_dict[pdb_file].add(pdb_pref + '_' + obj.chain_a + '.' + suff)
            pdb_files_dict[pdb_file].add(pdb_pref + '_' + obj.chain_b + '.' + suff)
        else:
            pdb_files_dict[pdb_file] = {pdb_file}
    xl_file_paths = sorted(xl_file_paths)
    pdb_files_lists = sorted(pdb_files_dict.items(), key=lambda item: item[0])
    pdb_files_lists = [sorted(tup[1]) for tup in pdb_files_lists]
    single_thread_extract_xl_features(xl_file_paths, pdb_files_lists, out_dir, skip_exist, predict)


def extract_xl_features(multi_thread=False, xl_files=None, pdb_files=None):
    cif_files = general_utils.load_obj('cif_files')
    xl_file_paths = []
    if pdb_files is None:
        pdb_files = []
        if xl_files is None:
            xl_files = os.listdir(OUTPUT_XL_FILES_PATH)
        for file in xl_files:
            xl_file_path = OUTPUT_XL_FILES_PATH + file
            uniport = file.split('.')[0]
            pdb_file = alpha_fold_files.generate_af_model_name(uniport)
            if pdb_file is None or not isfile(pdb_file):
                pref = 'pdb'
                suff = '.ent'
                if uniport in cif_files:
                    pref = ''
                    suff = '.cif'
                pdb_file = PDB_FILES_DIR + pref + uniport.lower() + suff
                if not isfile(pdb_file):
                    continue
            xl_file_paths.append(xl_file_path)
            pdb_files.append(pdb_file)
    if not multi_thread:
        single_thread_extract_xl_features(xl_file_paths, pdb_files)
        print("Finish all files")
        return
    available_cpus = len(os.sched_getaffinity(0)) - 1
    print("available cpus: ", available_cpus, flush=True)
    start_idx = 0
    step = int(len(xl_file_paths) / (available_cpus - 1))
    try:
        for pid in range(available_cpus):
            end = min(start_idx + step, len(xl_file_paths))
            multiprocessing.Process(target=single_thread_extract_xl_features,
                                    args=(xl_file_paths[start_idx: end], pdb_files[start_idx: end])).start()
            start_idx += step

    except Exception as e:
        print("Thread_Error")
        print(e)


def get_obj_files_key(obj, cif_files, inter_pdbs, by_chain=True):
    """
    Returns the file key of the xl object, which is used in different dictionaries and file name
    (feature and xl files)
    :param inter_pdbs:
    :param cif_files:
    :param obj:
    :param by_chain: if true - returns <pdb_key>_<chain_a>, <pdb_key>_<chain_b> else, returns <pdb_key>, None
    """
    pdb_file = obj.pdb_path
    if pdb_file == '':
        pdb_file = find_pdb_file_from_xl_obj(obj, cif_files, inter_pdbs)
        if pdb_file is None:
            pdb_file = find_pdb_file_from_xl_obj(obj, cif_files)
            if pdb_file is None:
                return None, None
    pdb_name = pdb_file.split('/')[-1].split('.')[0]
    if by_chain:
        if obj.chain_a != '' and obj.chain_b != '':
            pdb_name_a = f"{pdb_name}_{obj.chain_a}"
            pdb_name_b = f"{pdb_name}_{obj.chain_b}"
            return pdb_name_a, pdb_name_b
        else:
            return pdb_name, pdb_name
    return pdb_name, pdb_file


def find_feature_dict_keys_from_xl_obj(obj, feat_dict, problems, cif_files, inter_pdbs, feat_by_pdb=True, predict=False):
    if feat_by_pdb and not predict:
        uni_a, uni_b = get_obj_files_key(obj, cif_files, inter_pdbs, True)
        if uni_a not in feat_dict or uni_b not in feat_dict:
            uni_a, uni_b = obj.uniport_a.split('-')[0].split('.')[0], obj.uniport_b.split('-')[0].split('.')[0]
    else:
        uni_a, uni_b = obj.uniport_a.split('-')[0].split('.')[0], obj.uniport_b.split('-')[0].split('.')[0]
        if uni_a not in feat_dict or uni_b not in feat_dict:
            uni_a, uni_b = get_obj_files_key(obj, cif_files, inter_pdbs, True)
    return uni_a, uni_b


def read_features_to_dict_single_sample(feat_file, feat_dict, key):
    try:
        with open(feat_file, 'r') as f:
            for line in f:
                if line.split()[0] == "residue:":
                    cur_residue = line.split()[1]
                    feat_dict[key][cur_residue] = []
                else:
                    features = line.split(',')[:-1]
                    features = np.array(features, dtype=np.float64)
                    feat_dict[key][cur_residue].append(features)
    except Exception as e:
        print(f"Error reading {feat_file}")


def predict_read_features(pdb_names, feat_dict, feat_files, keys=None):
    # keys = ['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'J', 'K']
    for i, pdb in enumerate(pdb_names):
        if keys is None:
            key = pdb.split('/')[-1].split('.')[0]
            key = key.split('_')
            if key[-1] == 'tr':
                key = key[0]
            else:
                key = key[-1]
        else:
            key = keys[i]
        if key not in feat_dict:
            feat_dict[key] = dict()
        read_features_to_dict_single_sample(feat_files[i], feat_dict, key)


def read_features_to_dict(feat_dict=None, feat_files_path=XL_NEIGHBORS_FILES_PATH, out_path=XL_NEIGHBORS_FEATURE_DICT):
    """
    Read the output of C++ program xl_neighbor into dictionary of features (uniport, residue -> features)
    Saves dictionary - the input for graph dataset construction.
    :return:
    """
    print("Start reading features into dictionary")
    if feat_dict is None:
        feat_dict = {}
    print(f"Number of uniports in dict START: {len(feat_dict)}")
    for i, file in enumerate(os.listdir(feat_files_path)):
        key = file.split('.')[0]
        if key != '' and (key not in feat_dict or len(feat_dict[key]) == 0):
            feat_dict[key] = {}
            read_features_to_dict_single_sample(feat_files_path + file, feat_dict, key)
    general_utils.save_obj(feat_dict, out_path)
    print(f"Number of uniports in dict END: {len(feat_dict)}")
    print("Finish!")
    return feat_dict


def get_xl_neighbors_dict(dict_file=XL_NEIGHBORS_FEATURE_DICT):
    return general_utils.load_obj(dict_file)


def analyze_missing_files_of_objects(xl_objects):
    uniports = dict()
    inter_pdbs = get_inter_pdbs()
    pdbs = dict()
    cif_files = general_utils.load_obj('cif_files')
    for obj in xl_objects:
        pdb_file_name = find_pdb_file_from_xl_obj(obj, cif_files, inter_pdbs=inter_pdbs)
        if pdb_file_name is None:
            if obj.uniport_a == obj.uniport_b:
                if obj.uniport_a not in uniports:
                    uniports[obj.uniport_a] = 0
                uniports[obj.uniport_a] += 1
            if obj.pdb_file != '':
                if obj.pdb_file not in pdbs:
                    pdbs[obj.pdb_file] = 0
                pdbs[obj.pdb_file] += 1

    uniports = dict(sorted(uniports.items(), key=lambda item: item[1], reverse=True))
    pdbs = dict(sorted(pdbs.items(), key=lambda item: item[1], reverse=True))
    print(uniports)
    print(pdbs)
    print(f"total objects with missing uniports: {sum(uniports.values())}")
    print(f"total objects with missing pdbs: {sum(pdbs.values())}")


def read_pdb_info_jsons(dir_path='/cs/labs/dina/seanco/xl_parser/data/pdb_uniport_map/'):
    pdb_dict = dict()
    dup_pdbs = 0
    none_records = 0
    for file_name in os.listdir(dir_path):
        f = open(dir_path + file_name, )
        js = json.load(f)
        for record in js:
            pdb_key = record['data']['rcsb_id']
            if pdb_key not in pdb_dict:
                pdb_dict[pdb_key] = dict()
                for polymer in record['data']['polymer_entities']:
                    if polymer['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers'] is None:
                        none_records += 1
                        continue
                    if polymer['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers'][0]['database_name'] == 'UniProt':
                        entity_id = int(polymer['rcsb_polymer_entity_container_identifiers']['entity_id'])
                        pdb_dict[pdb_key][entity_id] = \
                            (polymer['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers'][0]['database_accession'])
            else:
                dup_pdbs += 1
    return pdb_dict


def parse_pdb_headers(dir_path='/cs/labs/dina/seanco/xl_parser/pdbs/', pdb_dict=None):
    if pdb_dict is None:
        pdb_dict = dict()
    print(f"start dict lengrh: {len(pdb_dict)}")
    for file_name in os.listdir(dir_path):
        if file_name in pdb_dict:
            continue
        if file_name.split('.')[-1] != 'cif':
            pdb_dict[file_name] = dict()
            command = f"grep COMPND {dir_path + file_name} | grep -E 'MOL_ID|CHAIN:'"
            result = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
            lines = str(result.stdout).strip().split('\\n')
            i = 0
            while i < len(lines) - 1:
                lines[i] = lines[i].strip()
                words = lines[i].split()
                try:
                    mol_id_idx = words.index('MOL_ID:')
                except ValueError as e:
                    i += 1
                    continue
                entity_id = words[mol_id_idx + 1].split(';')[0]
                if not entity_id.isdigit():
                    print(f"problem with entity_id: {entity_id} in pdb: {file_name}")
                i += 1
                chain_line = lines[i].strip().split()
                if len(chain_line) >= 2 and chain_line[2] == 'MOLECULE:':
                    chain_line = lines[i].split('\\n')[-1]
                chains = chain_line[3:]
                chains = [c.split(',')[0] for c in chains if c != '\\n\'']
                chains[-1] = chains[-1].split(';')[0]
                pdb_dict[file_name][entity_id] = chains
                i += 1
    print(f"end dict length: {len(pdb_dict)}")
    general_utils.save_obj(pdb_dict, 'pdb_entity_chain_dict')


def test_pdb_entity_dict():
    pdb_dict = get_pdb_entity_chain_dict()
    for key, val in pdb_dict.items():
        for id_, chains in val.items():
            if not id_.isdigit():
                print(f"problem with entity_id: {id_} in pdb: {key}")
            for c in chains:
                if len(c) > 2:
                    print(f"problem with chain: {c} in pdb: {key}")


def get_pdb_entity_chain_dict():
    return general_utils.load_obj('pdb_entity_chain_dict')


def get_inter_pdbs():
    return general_utils.load_obj('inter_pdbs')


def analyze_extracted_features(xl_objects):
    problems = [0, 0, 0, 0, 0]
    feat_dict = get_xl_neighbors_dict()
    cif_files = general_utils.load_obj('cif_files')
    inter_pdbs = get_inter_pdbs()
    feat_list = []
    dist_avg_list = []
    for obj in xl_objects:
        uni_a, uni_b = find_feature_dict_keys_from_xl_obj(obj, feat_dict, problems, cif_files, inter_pdbs)
        uniport_feat_dict_a = feat_dict[uni_a]
        uniport_feat_dict_b = feat_dict[uni_b]
        if obj.res_num_a not in uniport_feat_dict_a or obj.res_num_b not in uniport_feat_dict_b:
            problems[1] += 1
            # print(f"Problem type 1 with: {uni_a}, {uni_b} residues: {obj.res_num_a}, {obj.res_num_b}")
            continue
        res_a_feat = torch.from_numpy(np.stack(uniport_feat_dict_a[obj.res_num_a][:-1]))
        res_b_feat = torch.from_numpy(np.stack(uniport_feat_dict_b[obj.res_num_b][:-1]))
        feat_list.append(res_a_feat)
        feat_list.append(res_b_feat)
        dist_a, dist_b = uniport_feat_dict_a[obj.res_num_a][-1], uniport_feat_dict_b[obj.res_num_b][-1]
        dist_avg_list.append(np.mean(dist_a))
        dist_avg_list.append(np.mean(dist_b))
    feat_list = torch.cat(feat_list)
    feat_avg = torch.mean(feat_list, dim=0)
    feat_std = torch.std(feat_list, dim=0)
    dist_avg = np.mean(dist_avg_list)
    print(f"avg feature values: {feat_avg}")
    print(f"std feature values: {feat_std}")
    print(f"avg distances: {dist_avg}")
    print(problems)


def filter_cif_objects(xl_objects):
    cif_files = general_utils.load_obj('cif_files')
    inter_pdbs = get_inter_pdbs()
    cif_objects = []
    pdb_objects = []
    for obj in xl_objects:
        pdb = obj.pdb_path
        if obj.pdb_path == '':
            pdb = find_pdb_file_from_xl_obj(obj, cif_files, inter_pdbs)
        pdb = pdb.split('/')[-1]
        if pdb.split('.')[-1] == 'cif':
            cif_objects.append(obj)
        else:
            pdb_objects.append(obj)
    return pdb_objects, cif_objects


def num_of_chains_in_pdb(xl_objects, update=True):
    cif_files = general_utils.load_obj('cif_files')
    inter_pdbs = get_inter_pdbs()
    pdb_parser = PDBParser(PERMISSIVE=1)
    cif_parser = MMCIFParser()
    pdb_chains = dict()
    if update:
        pdb_chains = get_num_of_chains_in_pdb()
    for obj in xl_objects:
        if obj.uniport_a == obj.uniport_b and obj.pdb_file not in inter_pdbs:
            continue
        if obj.pdb_file not in pdb_chains:
            pdb_file_name = obj.pdb_path
            if obj.pdb_file in cif_files:
                structure = cif_parser.get_structure(pdb_file_name.split('/')[-1], pdb_file_name)
            else:
                structure = pdb_parser.get_structure(pdb_file_name.split('/')[-1], pdb_file_name)
            chains = list(structure.get_chains())
            pdb_chains[obj.pdb_file] = len(chains)
    general_utils.save_obj(pdb_chains, 'num_of_chains_by_pdb')
    return pdb_chains


def get_num_of_chains_in_pdb():
    return general_utils.load_obj('num_of_chains_by_pdb')


def get_xl_residues_list_by_chain(xl_file_name, chain_id):
    xl_residues = []
    with open(UNPROCESSED_XL_FILES_PATH + xl_file_name + ".txt", 'r') as f:
        for line in f:
            res_a, chain_a, res_b, chain_b, min_dist, max_dist = line.split()
            if chain_a == chain_id:
                xl_residues.append(int(res_a))
            if chain_b == chain_id:
                xl_residues.append(int(res_b))
    if len(xl_residues) > 1 and xl_residues[0] == 1:
        return xl_residues[1:]
    return xl_residues


def get_lysine_residues_list(pdb_path, pdb_parser, cif_parser, chain_id):
    print(pdb_path)
    if pdb_path.split('.')[-1] == 'cif':
        structure = cif_parser.get_structure(pdb_path.split('/')[-1], pdb_path)
    else:
        structure = pdb_parser.get_structure(pdb_path.split('/')[-1], pdb_path)
    chains = list(structure.get_chains())
    i = 0
    lys_list = []
    while i < len(chains):
        if chains[i].id == chain_id:
            lys_list = [int(res.id[1]) for res in chains[i] if res.resname == 'LYS']
            break
    return lys_list


def find_numbering_offset_single_chain(xl_offsets, pdb_path, xl_file_name, pdb_parser, cif_parser):
    pdb_name = pdb_path.split('/')[-1].split('.')[0]
    found_offset = INVALID_OFFSET
    if isfile(pdb_path):
        tmp_split = pdb_name.split('_')
        chain_id = 'A'
        if len(tmp_split) > 1:
            chain_id = tmp_split[-1]
        if pdb_name in xl_offsets:
            return
        xl_residues = get_xl_residues_list_by_chain(xl_file_name, chain_id)
        xl_residues.sort()
        xl_residues = np.unique(np.array(xl_residues))
        if len(xl_residues) > 1:
            lysine_residues = get_lysine_residues_list(pdb_path, pdb_parser, cif_parser, chain_id)
            lysine_residues.sort()
            lysine_residues = np.array(lysine_residues)
            gap_xl_residues = xl_residues - xl_residues[0]
            for i in range(len(lysine_residues) - len(xl_residues) + 1):
                tmp_lysine_residues = lysine_residues - lysine_residues[i]
                tmp_lysine_residues = tmp_lysine_residues[i:]
                intersection = np.intersect1d(tmp_lysine_residues, gap_xl_residues, assume_unique=True)
                if len(intersection) == len(gap_xl_residues):
                    offset = xl_residues[0] - lysine_residues[i]
                    if abs(offset) < abs(found_offset):
                        found_offset = offset  # searching for best match
                    else:
                        break  # The gap can only increase because elements sorted
    xl_offsets[pdb_name] = found_offset


def upd_obj_res_num_by_offset(obj, xl_offsets, pdb_path):
    pdb_name = pdb_path.split('/')[-1].split('.')[0]
    tmp_pdb_name = pdb_name
    if obj.uniport_a != pdb_name:
        tmp_pdb_name = pdb_name + "_" + obj.chain_a
    offset_a = xl_offsets[tmp_pdb_name]
    if offset_a != INVALID_OFFSET:
        obj.res_num_a = str(int(obj.res_num_a) - offset_a)
    if obj.chain_a == obj.chain_b:
        offset_b = offset_a
    else:
        tmp_pdb_name = pdb_name + "_" + obj.chain_b
        offset_b = xl_offsets[tmp_pdb_name]
    if offset_b != INVALID_OFFSET:
        obj.res_num_b = str(int(obj.res_num_b) - offset_b)
        return 1
    if offset_a != INVALID_OFFSET:
        return 1
    return 0


def create_xl_numbering_offset_dict(xl_objects, cif_files=None, inter_pdbs=None):
    if cif_files is None:
        cif_files = general_utils.load_obj('cif_files')
    if inter_pdbs is None:
        inter_pdbs = get_inter_pdbs()
    xl_offsets = dict()
    pdb_parser = PDBParser(PERMISSIVE=1)
    cif_parser = MMCIFParser()
    fixed_objects = []
    for obj in xl_objects:
        if obj.res_num_a == '-1' or obj.res_num_b == '-1':
            continue
        pdb_path = obj.pdb_path
        if pdb_path == '':
            pdb_path = find_pdb_file_from_xl_obj(obj, cif_files, inter_pdbs=inter_pdbs)
        pdb_name = pdb_path.split('/')[-1].split('.')[0]
        if obj.uniport_a != pdb_name:
            pref, suff = pdb_path.split('.')
            pdb_path_tmp = f"{pref}_{obj.chain_a}.{suff}"
            find_numbering_offset_single_chain(xl_offsets, pdb_path_tmp, pdb_name, pdb_parser, cif_parser)
            if obj.chain_a != obj.chain_b:
                pdb_path_tmp = f"{pref}_{obj.chain_b}.{suff}"
                find_numbering_offset_single_chain(xl_offsets, pdb_path_tmp, pdb_name, pdb_parser,cif_parser)
        else:
            find_numbering_offset_single_chain(xl_offsets, pdb_path, pdb_name, pdb_parser, cif_parser)
        res = upd_obj_res_num_by_offset(obj, xl_offsets, pdb_path)
        if res > 0:
            fixed_objects.append(obj)
    print(f"fixed objects: {len(fixed_objects)}")
    found_offset = [pdb for pdb, off in xl_offsets.items() if off != INVALID_OFFSET]
    print(f"found offsets: {len(found_offset)}")
    return fixed_objects


def filter_objects_from_list_by_pdb(xl_objects, inter_pdbs, filter_pdbs=None):
    cif_files = general_utils.load_obj('cif_files')
    filtered = []
    if filter_pdbs is None:
        filter_pdbs = {'1UJZ': {'A', 'B'}, '2BBM': {'A', 'B'}, '3WYO': {'A', 'B', 'C', 'D'}, '6JXD': {'E', 'F', 'A', 'B', 'C', 'D'},
                       '6NR8': {'L', 'P', 'O', 'J', 'K', 'H', 'E', 'G', 'M', 'N', 'D', 'B'}, '7BYI': {'A', 'B'}}
    for obj in xl_objects:
        pdb_file = obj.pdb_path
        if pdb_file is not None:
            pdb_name = pdb_file.split('/')[-1].split('.')[0]
            if pdb_name[:3] == 'pdb':
                pdb_name = pdb_name[3:]
            pdb_name = pdb_name.upper()
            if pdb_name in filter_pdbs:
                if obj.chain_a != obj.chain_b and obj.chain_a in filter_pdbs[pdb_name] and obj.chain_b in filter_pdbs[pdb_name]:
                    continue
            filtered.append(obj)
    print(f"start length: {len(xl_objects)}")
    print(f"end length: {len(filtered)}")
    return filtered


def create_dimers_dict(xl_objects, uni_chain_dict):
    entity_chain_dict = get_pdb_entity_chain_dict()
    inter_pdbs = get_inter_pdbs()
    cif_files = general_utils.load_obj('cif_files')
    prev_pdb = '-1'
    mmcif_dict = None
    cif_parser = MMCIFParser()
    dimer_dict = dict()
    for obj in xl_objects:
        if obj.pdb_file not in inter_pdbs:
            continue
        optional_chain_a, optional_chain_b = set(), set()
        pdb_file = find_pdb_file_from_xl_obj(obj, cif_files, inter_pdbs)
        if pdb_file is None:
            continue
        if pdb_file.split('.')[-1] == 'cif':
            if pdb_file != prev_pdb:
                structure = cif_parser.get_structure(pdb_file.split('/')[-1], pdb_file)
                mmcif_dict = cif_parser.get_mmcif_dict()
            find_chain_ids_cif(obj, mmcif_dict, optional_chain_a, optional_chain_b)
        else:
            find_chain_ids_pdb(obj, entity_chain_dict, optional_chain_a, optional_chain_b, uni_chain_dict, pdb_file.split('/')[-1])
        if len(optional_chain_a) > 1 or len(optional_chain_b) > 1:
            pdb_key = pdb_file.split('/')[-1].split('.')[0]
            if pdb_key not in dimer_dict:
                dimer_dict[pdb_key] = dict()
            if len(optional_chain_a) > 1:
                optional_chain_a = ''.join(sorted(optional_chain_a))
                for c in optional_chain_a:
                    dimer_dict[pdb_key][c] = optional_chain_a
            if obj.uniport_a != obj.uniport_b and len(optional_chain_b) > 1:
                optional_chain_b = ''.join(sorted(optional_chain_b))
                for c in optional_chain_b:
                    dimer_dict[pdb_key][c] = optional_chain_b
        prev_pdb = pdb_file
    general_utils.save_obj(dimer_dict, 'dimers_in_pdb_dict')
    return dimer_dict


def get_chain_by_chain_id(chains, chain_id):
    i = 0
    chain_res = 0
    while i < len(chains) and chain_res == 0:
        chain = chains[i]
        if chain.id == chain_id:
            chain_res = chain
        i += 1
    return chain_res


def get_seq_from_af_pdb(chain):
    seq = [seq1(aa.resname).lower() for aa in chain]
    return seq


def get_seq_from_chain(chain):
    seq = []
    prev_idx = 0
    for aa in chain:
        diff = int(aa.id[1]) - prev_idx
        if diff > 1:
            seq += ['_'] * (diff - 1)
        seq.append(seq1(aa.resname).lower())
        prev_idx = int(aa.id[1])
    return seq


def get_sequences_for_obj(obj, chains, seq_dict):
    if obj.chain_a not in seq_dict:
        chain_a = get_chain_by_chain_id(chains, obj.chain_a)
        if obj.pdb_path[-3:] == 'pdb':
            seq_a = get_seq_from_af_pdb(chain_a)
        else:
            seq_a = get_seq_from_chain(chain_a)
        seq_dict[obj.chain_a] = seq_a
    else:
        seq_a = seq_dict[obj.chain_a]
    if obj.chain_a == obj.chain_b:
        return seq_a, seq_a
    if obj.chain_b not in seq_dict:
        chain_b = get_chain_by_chain_id(chains, obj.chain_b)
        seq_b = get_seq_from_chain(chain_b)
    else:
        seq_b = seq_dict[obj.chain_b]
    return seq_a, seq_b

def create_scwrl_mutation_files(res_num, chain_id, closest_res, pdb_path, seq, out_dir):
    pdb_name = pdb_path.split('/')[-1].split('.')[0].split('_')[0]
    new_out_dir = f"{out_dir}{pdb_name}/{chain_id}/{res_num}/"
    print(new_out_dir)
    if os.path.isdir(new_out_dir):
        return new_out_dir
    os.makedirs(new_out_dir)
    seq[int(res_num) - 1] = seq[int(res_num) - 1].upper()
    for c in closest_res:
        seq[c[0] - 1] = seq[c[0] - 1].upper()  # to apply scwrl on these residues
    for i, c in enumerate(closest_res):
        for aa in AA_LETTERS:
            old_aa = c[1]
            if old_aa != aa:
                new_seq = copy.deepcopy(seq)
                new_seq[c[0] - 1] = aa
                seq_path = new_out_dir + str(c[0]) + '_' + old_aa + '_' + aa + '.txt'
                with open(seq_path, 'w') as f:
                    new_seq = "".join(new_seq)
                    new_seq = new_seq.replace('_', '')
                    f.write(new_seq)
                    f.close()
                new_pdb = new_out_dir + str(c[0]) + '_' + old_aa + '_' + aa + '.pdb'
                command = SCWRL_SCRIPT + f"-i {pdb_path} -o {new_pdb} -s {seq_path}"
                subprocess.run(command, shell=True)
    return new_out_dir


def create_dssp_files_from_list(to_create):
    for f in to_create:
        if not os.path.isfile(f + '.dssp'):
            dssp_dict_from_pdb_file(f, DSSP='/cs/staff/dina/software/Staccato/mkdssps')

def create_dssp_files(xl_objects):
    to_create = set()
    for obj in xl_objects:
        if obj.pdb_path[-3:] == 'pdb':
            to_create.add(obj.pdb_path)
        else:
            pdb_pref, suff = obj.pdb_path.split('.')
            to_create.add(pdb_pref + '_' + obj.chain_a + '.' + suff)
            to_create.add(pdb_pref + '_' + obj.chain_b + '.' + suff)
    create_dssp_files_from_list(to_create)



def create_fasta_from_pdb(pdb_path):
    pref, _ = pdb_path.split('.')
    fasta_path = f"{FASTA_PATH}{pref.split('/')[-1]}.fasta"
    if not os.path.isfile(fasta_path):
        cmd = f"{PDB_2_FASTA} {pdb_path} | tee {fasta_path}"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
    return fasta_path

def check_af_multimer(xl_objects):
    xl_objects = [obj for obj in xl_objects if obj.chain_a != obj.chain_b]
    couples = dict()
    for obj in xl_objects:
        if (obj.pdb_path, obj.chain_a, obj.chain_b) in couples or (obj.pdb_path, obj.chain_b, obj.chain_a) in couples:
            continue
        pref, suff = obj.pdb_path.split('.')
        pdb_a = f"{pref}_{obj.chain_a}.{suff}"
        pdb_b = f"{pref}_{obj.chain_b}.{suff}"
        fasta_a = create_fasta_from_pdb(pdb_a)
        fasta_b = create_fasta_from_pdb(pdb_b)
        merged_fasta = alpha_fold_files.merge_fastas_to_af_multimer(fasta_a, fasta_b, FASTA_PATH)
        couples[(obj.pdb_path, obj.chain_a, obj.chain_b)] = merged_fasta
    general_utils.save_obj(couples, "multimer_fastas_by_objects")

