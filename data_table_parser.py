
import numpy as np
import urllib.request
import urllib.error
from general_utils import save_obj
from general_utils import load_obj
import general_utils
import cross_link
import pdb_files_manager


PEP_A = 0
POS_IN_PEP_A = 1
PEP_B = 12
POS_IN_PEP_B = 13
STRUCTURES = -1
ACCESSION_A = 3
UNIPORT_A = 4
ACCESSION_B = 15
UNIPORT_B = 16
RES_NUM_A = 5
RES_NUM_B = 17
DISTANCE = 24
DATASET = 28
ALL_IDX = [PEP_A, POS_IN_PEP_A, UNIPORT_A, RES_NUM_A,
           PEP_B, POS_IN_PEP_B, UNIPORT_B, RES_NUM_B, DISTANCE, STRUCTURES]
FASTA_IDX = [PEP_A, POS_IN_PEP_A, UNIPORT_A, PEP_B, POS_IN_PEP_B, UNIPORT_B]
FASTA_SUFFIX = ".fasta"
HEADLINE = "AvailableStructures\n"
UNIPORT_SERVER = "https://www.uniprot.org/uniprot/"


def upd_fasta_dict_from_file(file_path='/cs/labs/dina/seanco/xl_parser/data/ACP_Percoll_Fraction_XL_Bruce.fasta'):
    cur_dict = load_obj(general_utils.FASTA_DICT_NAME)
    new_seqs = general_utils.read_fasta_files(file_path)
    cur_dict.update(new_seqs)
    save_obj(cur_dict, general_utils.FASTA_DICT_NAME)


def get_fasta(uniport):
    cur_url = UNIPORT_SERVER + uniport + FASTA_SUFFIX
    fasta = 0
    try:
        with urllib.request.urlopen(cur_url) as f:
            lines = (f.read().decode('utf-8').split('\n'))[1::]
            fasta = ''.join(lines)
    except urllib.error.URLError as e:
        print(e.reason)
        print('\n' + uniport)

    return fasta


def update_uniport_fasta_dict_from_xl_objects():
    xl_list = cross_link.CrossLink.load_all_xl_objects()
    pdb_files_manager.update_xl_objects_with_obsolete_uniports(xl_list)
    cur_dict = load_obj(general_utils.FASTA_DICT_NAME)
    start_len = len(cur_dict)
    for obj in xl_list:
        try:
            if obj.uniport_a not in cur_dict or \
                    (obj.uniport_a in cur_dict and len(cur_dict[obj.uniport_a]) <= 0):
                fasta = get_fasta(obj.uniport_a)
                cur_dict[obj.uniport_a] = fasta
            if obj.uniport_a != obj.uniport_b:
                if obj.uniport_b not in cur_dict or \
                        (obj.uniport_b in cur_dict and len(cur_dict[obj.uniport_b]) <= 0):
                    fasta = get_fasta(obj.uniport_b)
                    cur_dict[obj.uniport_b] = fasta
        except Exception as e:
            print(e)
            continue
    print("added fastas: ", len(cur_dict) - start_len)
    save_obj(cur_dict, general_utils.FASTA_DICT_NAME)


def fix_fasta_int_values():
    cur_dict = load_obj(general_utils.FASTA_DICT_NAME)
    to_remove = []
    upd = 0
    for key, val in cur_dict.items():
        try:
            if isinstance(val, int):
                fasta = get_fasta(key)
                if not isinstance(fasta, int) and len(fasta) > 0:
                    cur_dict[key] = fasta
                    upd += 1
                else:
                    to_remove.append(key)
        except Exception as e:
            print(e)
            continue

    print("removed: ", len(to_remove))
    for key in to_remove:
        del(cur_dict[key])
    save_obj(cur_dict, general_utils.FASTA_DICT_NAME)
    print("fixed: ", upd)


def update_uniport_fasta_dict():
    cur_dict = load_obj(general_utils.FASTA_DICT_NAME)
    xl_list = load_obj('xl_with_dup_all_cols')
    for sample in xl_list:
        try:
            if sample[ACCESSION_A] not in cur_dict or (sample[ACCESSION_A] in cur_dict and len(cur_dict[sample[ACCESSION_A]]) <= 0):
                fasta = get_fasta(sample[UNIPORT_A])
                cur_dict[sample[ACCESSION_A]] = fasta
        except Exception as e:
            print(e)
            continue
    save_obj(cur_dict, 'fasta_files_dict_by_uniport')


def convert_fasta_dict_key():
    fasta_dict = load_obj('fasta_files_dict')
    new_fasta_dict = {}
    xl_list = load_obj('xl_with_dup_all_cols')
    for sample in xl_list:
        if sample[UNIPORT_A] in fasta_dict and sample[ACCESSION_A] not in new_fasta_dict:
            new_fasta_dict[sample[ACCESSION_A]] = fasta_dict[sample[UNIPORT_A]]

    print(len(fasta_dict))
    print(len(new_fasta_dict))
    save_obj(new_fasta_dict, general_utils.FASTA_DICT_NAME)


def create_fasta_dict(path='/cs/labs/dina/seanco/xl_parser/xldb_data/all_data.txt'):
    input_file = open(path)
    idx = np.array(FASTA_IDX)
    xl_fasta_dict = {}
    xl_samples = []
    i = 0
    for line in input_file:
        cols = np.array(line.split("\t"))
        if cols[STRUCTURES] != HEADLINE:
            xl_samples.append(np.asarray(cols[idx]))
            if cols[ACCESSION_A] not in xl_fasta_dict:
                fasta_a = get_fasta(cols[ACCESSION_A])
                xl_fasta_dict[cols[ACCESSION_A]] = fasta_a
            if cols[ACCESSION_B] not in xl_fasta_dict:
                fasta_b = get_fasta(cols[ACCESSION_B])
                xl_fasta_dict[cols[ACCESSION_B]] = fasta_b
        else:
            i += 1
            print(line)

    print("number of fastas: " + str(len(xl_fasta_dict)) + '\n')
    print("number of samples: " + str(len(xl_samples)) + '\n')
    print(f"number of problems: {i} ")
    save_obj(xl_fasta_dict, "fasta_files_dict")
    # save_obj(xl_samples, "all_samples")


def read_xl(path):
    input_file = open(path)
    idx = np.array(ALL_IDX)
    xl_with_struct = []
    for line in input_file:
        cols = np.array(line.split("\t"))
        if cols.size > STRUCTURES and cols[STRUCTURES] != "\n":
            xl_with_struct.append(np.asarray(cols[idx]))
    return xl_with_struct


def read_datasets_names(path='/cs/labs/dina/seanco/xl_parser/xldb_data/all_data.txt'):
    with open(path, 'r') as input_file:
        datasets_dict = {}
        for line in input_file:
            cols = np.array(line.split("\t"))
            if cols[STRUCTURES] != "AvailableStructures\n":
                datasets = cols[DATASET].split(',')
                for dataset in datasets:
                    if dataset not in datasets_dict:
                        datasets_dict[dataset] = 1
                    else:
                        datasets_dict[dataset] += 1

        with open('/cs/labs/dina/seanco/xl_parser/data/dataset_names.txt', 'w') as out_file:
            for key, val in sorted(datasets_dict.items(), key=lambda item: item[1]):
                out_file.write(key + '\t' + str(val) + '\n')
            out_file.close()
        input_file.close()


def read_all_clear_dup(path='/cs/labs/dina/seanco/xl_parser/data/xl_db_data_update.txt'):
    input_file = open(path)
    xl_dict = {}
    # idx = np.array(ALL_IDX)
    xl_with_struct = []
    num_of_structures = 0
    unredundant_xl = 0
    dup = 0
    for line in input_file:
        cols = np.array(line.split("\t"))
        if cols[STRUCTURES] != "AvailableStructures\n":
            # if cols[PEP_A] == cols[PEP_B] and cols[POS_IN_PEP_A] == cols[POS_IN_PEP_B]:
            #     continue
            structures = cols[STRUCTURES][:-1]
            key_a = cols[UNIPORT_A] + '+' + cols[PEP_A] + '+' + cols[RES_NUM_A]
            key_b = cols[UNIPORT_B] + '+' + cols[PEP_B] + '+' + cols[RES_NUM_B]
            if key_a in xl_dict:
                if key_b in xl_dict[key_a]:
                    if structures in xl_dict[key_a]:
                        dup += 1
                        xl_with_struct.append(np.asarray(cols))
                        continue
                    xl_dict[key_a].add(structures)
                    xl_dict[key_b].add(structures)
                    xl_with_struct.append(np.asarray(cols))
                    # xl_with_struct.append(np.asarray(cols[idx]))
                    xl_with_struct[-1][-1] = xl_with_struct[-1][-1][:-1]
                else:
                    unredundant_xl += 1
                    xl_with_struct.append(np.asarray(cols))
                    # xl_with_struct.append(np.asarray(cols[idx]))
                    xl_with_struct[-1][-1] = xl_with_struct[-1][-1][:-1]
                    num_of_structures += len(structures.split(","))
                    if key_b not in xl_dict:
                        xl_dict[key_b] = {key_a, structures}
                    else:
                        xl_dict[key_b].add(key_a)
                        xl_dict[key_b].add(structures)
                    xl_dict[key_a].add(key_b)
                    xl_dict[key_a].add(structures)
            else:
                unredundant_xl += 1
                xl_with_struct.append(np.asarray(cols))
                # xl_with_struct.append(np.asarray(cols[idx]))
                xl_with_struct[-1][-1] = xl_with_struct[-1][-1][:-1]
                num_of_structures += len(structures.split(","))
                xl_dict[key_a] = {key_b, structures}
                xl_dict[key_b] = {key_a, structures}

    print(f"unredundant: {unredundant_xl} duplicates: {dup}")
    print("num of samples: " + str(len(xl_with_struct)))
    save_obj(xl_with_struct, "xl_with_dup_all_cols")
    return xl_with_struct


def linker_types_parser():
    types_dict = {}
    with open('/cs/labs/dina/seanco/xl_parser/xldb_data/linker_types.txt', 'r') as f:
        first_line = True
        for line in f:
            if first_line:
                first_line = False
                continue
            data = line.split('\t')
            types_dict[data[0]] = data[1]

        save_obj(types_dict, 'linker_types')


# def main():
#     # data_path = '/cs/labs/dina/seanco/xl_parser/data/xl_db_data_update.txt'
#     # input_path = sys.argv[1]
#     # output_path = sys.argv[2]
#     # read_all_clear_dup()
#     # read_all_clear_dup('/cs/labs/dina/seanco/xl_parser/xldb_data/xl_data_of_af_organisms.txt')
#     # create_fasta_dict(data_path)
#     # convert_fasta_dict_key()
#     # update_uniport_fasta_dict()
#     # linker_types_parser()
#     # read_datasets_names(data_path)
#     # update_uniport_fasta_dict_from_xl_objects()
#     # fix_fasta_int_values()
#     upd_fasta_dict_from_file()
#
#
# if __name__ == "__main__":
#     main()



