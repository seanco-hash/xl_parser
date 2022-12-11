
import requests
import general_utils
from os.path import isfile
from os import listdir
import os
import json
import cross_link
import numpy as np
import shutil

AF_PDB_DIR = '/cs/labs/dina/seanco/xl_parser/pdbs/alpha_fold/pdb_files/'
AF_PRED_ERROR_DIR = '/cs/labs/dina/seanco/xl_parser/pdbs/alpha_fold/predict_error/'

AF_PDB_URL_PREFIX = 'https://alphafold.ebi.ac.uk/files/AF-'
AF_PDB_URL_SUFFIX = '-F1-model_v3.pdb'
AF_PRED_ERR_URL_SUFFIX = '-F1-predicted_aligned_error_v3.json'

PDB_FILE_FORMAT = '.pdb'
PRED_ERR_FILE_FORMAT = '.json'
PRED_ERR_FILE_FORMAT_PICKLE = '.pkl'


def generate_af_model_name(accession):
    acc = accession.split('-')[0]  # For accessions like p0023-1
    acc = acc.split('.')[0]
    pdb_name = AF_PDB_DIR + acc + PDB_FILE_FORMAT
    return pdb_name


def get_pae_by_pickle(accession, res_a, res_b, n_residues, file_name):
    with open(file_name, 'rb') as f:
        feat_dict = general_utils.pickle.load(f)
        pae = feat_dict['paes']['model_1']
        if pae.shape == (n_residues, n_residues) and int(res_a) <= n_residues and int(res_b) <= n_residues:
            pred_err = pae[int(res_a) - 1, int(res_b) - 1]
            return pred_err
        print(accession)
        return -1


def get_pae_by_json(accession, res_a, res_b, n_residues, file_name):
    ret_val = 0
    f = open(file_name, )
    js = json.load(f)[0]
    if 'predicted_aligned_error' in js and len(js['predicted_aligned_error']) == n_residues and \
            len(js['predicted_aligned_error'][0]) == n_residues and \
            int(res_a) <= n_residues and int(res_b) <= n_residues:
        dist_matrix = np.asarray(js['predicted_aligned_error']).reshape((n_residues, n_residues))
    elif 'distance' in js and len(js['distance']) == n_residues * n_residues and \
            int(res_a) <= n_residues and int(res_b) <= n_residues:
        dist_matrix = np.asarray(js['distance']).reshape((n_residues, n_residues))
    else:
        print(accession)
        if 'predicted_aligned_error' in js:
            print(f"problem with PAE. residues: {res_a} , {res_b} , {n_residues} , distances: {len(js['predicted_aligned_error'])}")
        else:
            print(f"problem with PAE. residues: {res_a} , {res_b} , {n_residues} , distances: unknown")

        return -1
    pred_err = ((dist_matrix[int(res_a) - 1, int(res_b) - 1]) +
                (dist_matrix[int(res_b) - 1, int(res_a) - 1]) / 2)
    return pred_err


def get_pred_align_err(accession, res_a, res_b, n_residues):
    try:
        accession = accession.split('-')[0]
        accession = accession.split('.')[0]
        file_name = AF_PRED_ERROR_DIR + accession + PRED_ERR_FILE_FORMAT
        if isfile(file_name):
            return get_pae_by_json(accession, res_a, res_b, n_residues, file_name)
        else:
            file_name = AF_PRED_ERROR_DIR + accession + PRED_ERR_FILE_FORMAT_PICKLE
            if isfile(file_name):
                return get_pae_by_pickle(accession, res_a, res_b, n_residues, file_name)
            else:
                return -1
    except OSError as e:
        return -1
    except ValueError as e:
        return -1


def download_pae(uniport):
    pred_err_url = AF_PDB_URL_PREFIX + uniport + AF_PRED_ERR_URL_SUFFIX
    pred_er_file = requests.get(pred_err_url, allow_redirects=True)
    if pred_er_file.status_code != 404:
        open(AF_PRED_ERROR_DIR + uniport + PRED_ERR_FILE_FORMAT, 'wb').write(pred_er_file.content)
        return 0
    return 1


def get__af_file(uniport='A0A096LP55'):
    try:
        future_file_name = AF_PDB_DIR + uniport + PDB_FILE_FORMAT
        if isfile(future_file_name):
            return 1
        pdb_url = AF_PDB_URL_PREFIX + uniport + AF_PDB_URL_SUFFIX
        pdb_file = requests.get(pdb_url, allow_redirects=True)
        if pdb_file.status_code != 404:
            open(future_file_name, 'wb').write(pdb_file.content)
            download_pae(uniport)
        else:
            print("Not found model for uniport: " + uniport)
            return 0

    except Exception as e:
        print(e)
        return 0

    return 1


def get_xl_objects_af_pdb():
    xl_list = cross_link.CrossLink.load_all_xl_objects_as_list()
    print("list length: ", len(xl_list))
    already_downloaded = set()
    already_exist = 0
    with open("uniports_to_predict.txt", 'r') as missing_uniports_file:
        for line in missing_uniports_file:
            already_downloaded.add(line[:-1])
    files_downloaded = 0
    dup_uniports = 0
    with open("uniports_to_predict.txt", 'a') as missing_uniports_file:
        for sample in xl_list:
            if sample[cross_link.UNIPORT_A] == sample[cross_link.UNIPORT_B]:
                uniport = sample[cross_link.UNIPORT_A]
                uni_split = uniport.split('-')
                if len(uni_split) > 1:
                    uniport = uni_split[0]
                if uniport not in already_downloaded:
                    already_downloaded.add(uniport)
                    existing_model_name = generate_af_model_name(uniport)
                    if not isfile(existing_model_name):
                        res = get__af_file(uniport)
                        files_downloaded += res
                        if res == 0:
                            # if uniport not in to_ommit:
                            missing_uniports_file.write(uniport + '\n')
                            print(files_downloaded + already_exist)
                    else:
                        already_exist += 1
            else:
                dup_uniports += 1

        print("list length: ", len(xl_list))
        print("duplicates_uniports: ", dup_uniports)
        print('Files downloaded: ' + str(files_downloaded))
        print("already exist: ", already_exist)


def get_intra_xl_af_pdb(xl_list):
    already_downloaded = set()
    total_uniports = 0
    files_downloaded = 0
    already_exist = 0
    # problematic_uni = open("alphafold/pronlematic_uniport.txt", 'r')
    # to_ommit = set()
    # for uni in problematic_uni.readlines():
    #     to_ommit.add(uni[:-1])
    with open("uniports_to_predict.txt", 'w') as missing_uniports_file:
        for sample in xl_list:
            if sample[general_utils.XL_UNIPORT_A] == sample[general_utils.XL_UNIPORT_B] and sample[
                general_utils.XL_PDB_TYPE_A] != \
                    general_utils.XL_AVAILABLE_IN_PDB:
                uniport = sample[general_utils.XL_ACCESSION_A]
                uni_split = uniport.split('-')
                if len(uni_split) > 1:
                    uniport = uni_split[0]
                if uniport not in already_downloaded:
                    already_downloaded.add(uniport)
                    total_uniports += 1
                    existing_model_name = generate_af_model_name(uniport)
                    if not isfile(existing_model_name):
                        res = get__af_file(uniport)
                        files_downloaded += res
                        if res == 0:
                            # if uniport not in to_ommit:
                            missing_uniports_file.write(uniport + '\n')
                    else:
                        already_exist += 1

        print('Files downloaded: ' + str(files_downloaded))
        print("Total uniports: ", total_uniports)
        print("already exist: ", already_exist)


def get_af_of_cross_link_objects(xl_objects):
    already_exist = set()
    missing_uniprots= set()
    for cl in xl_objects:
        # if cl.pdb_file is None or cl.pdb_file == "" or cl.pdb_file == "None":
        if cl.uniport_a not in already_exist:
            res = get__af_file(cl.uniport_a.split('-')[0])
            already_exist.add(cl.uniport_a)
            if res != 1:
                missing_uniprots.add(cl.uniport_a)
        if cl.uniport_a != cl.uniport_b and cl.uniport_b not in already_exist:
            res = get__af_file(cl.uniport_b.split('-')[0])
            already_exist.add(cl.uniport_b)
            if res != 1:
                missing_uniprots.add(cl.uniport_b)

    with open("uniports_to_predict.txt", 'w') as missing_uniports_file:
        for uniport in missing_uniprots:
            missing_uniports_file.write(uniport + '\n')

    print("missing files: ", len(missing_uniprots))


def check_missing_af_files():
    fasta_dict = general_utils.load_obj('fasta_files_dict_by_uniport')
    missing_fasta = 0
    with open('/cs/labs/dina/seanco/xl_parser/uniports_to_predict.txt', 'r') as f:
        uniports = [line[:-1] for line in f]
        for uniport in uniports:
            if uniport in fasta_dict:
                seq = fasta_dict[uniport]
                x = len(seq)
                if x <= 1400:
                    print(f"uniport: {uniport}, length: {x}")
            else:
                missing_fasta += 1
                # print(uniport)
        print(len(uniports))
        print(missing_fasta)


def remove_predicted_seqs_from_fasta(dir_path='/cs/labs/dina/seanco/xl_parser/alphafold/slurm/'):
    fasta_path = "/cs/labs/dina/seanco/xl_parser/fasta_to_predict.fasta"
    done_file_suffix = '.done.txt'
    seqs_dict = general_utils.read_fasta_files(fasta_path)
    new_dict = dict()
    # len_below = 0
    # len_above = 0
    for uniport, seq in seqs_dict.items():
        if not isfile(dir_path + uniport + done_file_suffix):
            new_dict[uniport] = seq
    # print(f"short seqs: {len_below}")
    # print(f"long seqs {len_above}")
    general_utils.seqs_dict_to_fasta(new_dict, fasta_path)


def organize_local_alphfold_output(out_dir='/cs/labs/dina/seanco/xl_parser/alphafold/slurm/', dest_dir="", files_spliit_len=6):
    to_remove = []
    copied = 0
    for file in listdir(out_dir):
        suff = file.split('.')[-1]
        words = file.split('_')
        if (len(words) == files_spliit_len and suff == 'pdb' and words[3] == '1') or \
                (len(words) == files_spliit_len - 1 and suff == 'json'):
            if words[0][0] == 'M':
                words[0] = words[0][0] + words[0][4:]
            if suff == 'pdb':
                dest_dir = AF_PDB_DIR
            else:
                dest_dir = AF_PRED_ERROR_DIR
            new_file_name = dest_dir + words[0] + '.' + suff
            shutil.copy(out_dir + file, new_file_name)
            copied += 1
        if suff != '.txt':
            to_remove.append(file)

    for file in to_remove:
        if not os.path.isdir(out_dir + file):
            os.remove(out_dir + file)


def fix_json_files_dir():
    to_remove = []
    for f in listdir(AF_PDB_DIR):
        suff = f.split('.')[-1]
        if suff == 'json':
            shutil.copy(AF_PDB_DIR + f, AF_PRED_ERROR_DIR)
            to_remove.append(f)

    for f in to_remove:
        os.remove(AF_PDB_DIR + f)


# def main():
#
#      # xl_list = load_obj('xl_with_dup_all_cols')
#      # just_test()
#      # xl_list = load_obj('xl_with_dup_all_cols')
#      # get_intra_xl_af_pdb(xl_list)
#      # get_xl_objects_af_pdb()
#
#
# if __name__ == "__main__":
#     main()
# get_pae_by_json('Q16881', 15, 30, 649, '/cs/labs/dina/seanco/xl_parser/pdbs/alpha_fold/predict_error/Q16881.json')