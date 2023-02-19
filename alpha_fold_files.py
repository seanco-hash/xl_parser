
import requests
import general_utils
from os.path import isfile
from os import listdir
import os
import json
import cross_link
import numpy as np
import shutil
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import matplotlib.pyplot as plt


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
    pred_err = (((dist_matrix[int(res_a) - 1, int(res_b) - 1]) +
                (dist_matrix[int(res_b) - 1, int(res_a) - 1])) / 2)
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


def organize_local_alphfold_output(out_dir='/cs/labs/dina/seanco/xl_parser/alphafold/slurm/', dest_dir="",
                                   files_spliit_len=6, model_rank_idx=3):
    to_remove = []
    copied = 0
    for file in listdir(out_dir):
        suff = file.split('.')[-1]
        words = file.split('_')
        if (len(words) == files_spliit_len and suff == 'pdb' and words[model_rank_idx] == '1') or \
                (len(words) == files_spliit_len - 1 and suff == 'json'):
            if words[0][0] == 'M':
                words[0] = words[0][0] + words[0][4:]
            if suff == 'pdb':
                dest_dir = AF_PDB_DIR
            else:
                dest_dir = AF_PRED_ERROR_DIR
            new_file_name = f"{dest_dir}{words[0]}_{words[1]}_{words[2]}.{suff}"
            # new_file_name = dest_dir + words[0] + '.' + suff
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


def plot_error_vs_distance(xl_objects):
    distances = [obj.distance for obj in xl_objects if obj.error != cross_link.INVALID_ERROR_VALUE]
    errors = [obj.error for obj in xl_objects if obj.error != cross_link.INVALID_ERROR_VALUE]
    print(distances)
    print(errors)
    general_utils.plot_scatter(distances, errors, '', 'distance', 'error', mark_low_th=0)


def merge_fastas_to_af_multimer(fasta_a, fasta_b, fasta_dir):
    f1 = open(fasta_a, 'r')
    _ = f1.readline()
    seq1 = ''.join([s[:-1] for s in f1])
    f2 = open(fasta_b, 'r')
    _ = f2.readline()
    seq2 = ''.join([s[:-1] for s in f2])
    f1.close()
    f2.close()
    if len(seq1) + len(seq2) < 2000:
        af_predict_fasta_dir = fasta_dir + 'af_predict/'
        seq = seq1 + ':' + seq2
        new_fasta = f"{af_predict_fasta_dir}{fasta_a.split('.')[0].split('/')[-1]}_{fasta_b.split('.')[0].split('_')[-1]}.fasta"
        with open(new_fasta, 'w') as f:
            f.write(f">{new_fasta}\n")
            f.write(seq)
        return new_fasta
    return None


def get_obj_distance_af_multimer(xl_objects):
    pdb_parser = PDBParser(PERMISSIVE=1)
    cif_parser = MMCIFParser()
    orig_dist = []
    multi_dist = []
    pae_errors = []
    for xl_obj in xl_objects:
        multi_res_b, multi_res_a = 0, 0
        multimer_pdb = f"{AF_PDB_DIR}{xl_obj.pdb_path.split('/')[-1].split('.')[0]}_{xl_obj.chain_a}_{xl_obj.chain_b}.pdb"
        if not os.path.isfile(multimer_pdb):
            multimer_pdb = f"{AF_PDB_DIR}{xl_obj.pdb_path.split('/')[-1].split('.')[0]}_{xl_obj.chain_b}_{xl_obj.chain_a}.pdb"
            if not os.path.isfile(multimer_pdb):
                continue
        structure_multimer = pdb_parser.get_structure(multimer_pdb.split('/')[-1], multimer_pdb)
        chains_multimer = list(structure_multimer.get_chains())
        if xl_obj.pdb_path[-3:] == 'cif':
            structure_orig = cif_parser.get_structure(xl_obj.pdb_path.split('/')[-1], xl_obj.pdb_path)
        else:
            structure_orig = pdb_parser.get_structure(xl_obj.pdb_path.split('/')[-1], xl_obj.pdb_path)
        chains_orig = list(structure_orig.get_chains())
        chain_a, chain_b = None, None
        i = 0
        while i < len(chains_orig) and (chain_a is None or chain_b is None):
            if chains_orig[i].id == xl_obj.chain_a:
                chain_a = chains_orig[i]
            if chains_orig[i].id == xl_obj.chain_b:
                chain_b = chains_orig[i]
            i += 1
        if chain_a is None or chain_b is None:
            print(f"problem with chains: {xl_obj.pdb_path}, {xl_obj.chain_a}, {xl_obj.chain_b}")
            continue
        if len(chain_a.child_list) == len(chains_multimer[0].child_list) and len(chain_b.child_list) == len(chains_multimer[1].child_list):
            pairs_chains = [[chain_a, chains_multimer[0], multi_res_a, int(xl_obj.res_num_a)], [chain_b, chains_multimer[1], multi_res_b, int(xl_obj.res_num_b)]]
        elif len(chain_a.child_list) == len(chains_multimer[1].child_list) and len(chain_b.child_list) == len(chains_multimer[0].child_list):
            pairs_chains = [[chain_a, chains_multimer[1], multi_res_a, int(xl_obj.res_num_a)], [chain_b, chains_multimer[0], multi_res_b, int(xl_obj.res_num_b)]]
        else:
            continue
        for p_chains in pairs_chains:
            i = 0
            while i < len(p_chains[0].child_list) and p_chains[2] == 0:
                aa = p_chains[0].child_list[i]
                if aa.id[1] == p_chains[3]:
                    p_chains[2] = p_chains[1].child_list[i]
                i += 1
        dist = cross_link.CrossLink.get_atom_distance_form_residues(pairs_chains[0][2], pairs_chains[1][2], 'CA', 'CA')
        json_file = f"{AF_PRED_ERROR_DIR}{multimer_pdb.split('/')[-1].split('.')[0]}{PRED_ERR_FILE_FORMAT}"
        err = get_pae_by_json(xl_obj.pdb_path, str(pairs_chains[0][2].id[1]), pairs_chains[1][2].id[1] + len(pairs_chains[0][1]), len(pairs_chains[0][1]) + len(pairs_chains[1][1]), json_file)
        if dist != -1:
            print(f"orig dist: {xl_obj.distance}, multimer dist: {dist}, PAE err {err}")
            orig_dist.append(xl_obj.distance)
            multi_dist.append(dist)
            pae_errors.append(err)
        else:
            print(f"problem calculating distance {xl_obj.pdb_path}")

    diff = np.abs(np.array(orig_dist) - np.array(multi_dist))
    pae_errors = np.array(pae_errors)
    res = [orig_dist, multi_dist, pae_errors]
    general_utils.save_obj(res, 'af_multimer_distances')
    print(diff)
    print(f"diff mean: {np.mean(diff)}, std: {np.std(diff)}")
    plt.scatter(diff, pae_errors)
    plt.show()


def plot_multimer_dist_err():
    res = general_utils.load_obj('af_multimer_distances')
    orig_dist, multi_dist, pae_errors = res
    orig_dist, multi_dist, pae_errors = np.array(orig_dist), np.array(multi_dist), np.array(pae_errors)
    # multi_dist = multi_dist[np.where(orig_dist < 20)]
    # pae_errors = pae_errors[np.where(orig_dist < 20)]
    # orig_dist = orig_dist[np.where(orig_dist < 20)]
    # diff = np.abs(orig_dist - multi_dist)
    # diff[diff > 30] = 30
    # multi_dist[multi_dist > 60] = 60
    fig, ax = plt.subplots()
    sc = ax.scatter(multi_dist, orig_dist, c=pae_errors, cmap=plt.cm.get_cmap('RdYlBu'))
    ax.set_ylim([0, 60])
    ax.set_xlim([0, max(multi_dist) + 5])

    # Add black line where x == y
    # min_val = min(np.min(orig_dist), np.min(multi_dist))
    # max_val = 45
    # ax.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--')
    cbar = fig.colorbar(sc, ax=ax, orientation='horizontal', pad=0.15, aspect=50)
    cbar.ax.xaxis.set_ticks_position('bottom')
    # cbar = plt.colorbar(sc, orientation='horizontal')
    cbar.ax.set_ylabel('PAE', fontsize=12)
    # Add legend to show what colors represent
    handles, labels = sc.legend_elements(prop="colors", alpha=0.6)
    # legend = ax.legend(handles, labels, loc="upper right", title="AF PAE")

    # Set axis labels and title
    ax.set_ylabel('PDB distance', fontsize=12)
    ax.set_xlabel('AFM distance', fontsize=12)
    plt.gca().set_aspect('equal', adjustable='box')
    # ax.set_title('Crosslink distance in AF multimer vs PDB, colored by PAE', fontsize=14)
    plt.show()


def fix_objects_pae_error(xl_objects):
    pdb_parser = PDBParser(PERMISSIVE=1)
    xl_objects = sorted(xl_objects, key=lambda x: x.pdb_path)
    prev_pdb = ''
    for i, obj in enumerate(xl_objects):
        if obj.pdb_path[-3:] == 'pdb':
            if obj.pdb_path != prev_pdb:
                structure = pdb_parser.get_structure(obj.pdb_path.split('/')[-1], obj.pdb_path)
                chains = list(structure.get_chains())
            # err_file = f"{AF_PRED_ERROR_DIR}{obj.pdb_path.split('/')[-1].split('.')[0]}{PRED_ERR_FILE_FORMAT}"
            # err = get_pae_by_json(obj.uniport_a, obj.res_num_a, obj.res_num_b, len(chains[0]), err_file)
            err = get_pred_align_err(obj.uniport_a, obj.res_num_a, obj.res_num_b, len(chains[0]))
            obj.error = err
            prev_pdb = obj.pdb_path
            if i % 100 == 0:
                print(i)
    plot_error_vs_distance(xl_objects)

