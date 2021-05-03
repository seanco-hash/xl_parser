import time
import sys
import numpy as np
import urllib.request
import urllib.error
from general_utils import save_obj

PEP_A = 0
POS_IN_PEP_A = 1
PEP_B = 12
POS_IN_PEP_B = 13
STRUCTURES = -1
UNIPORT_A = 4
UNIPORT_B = 16
ALL_IDX = [PEP_A, POS_IN_PEP_A, UNIPORT_A, PEP_B, POS_IN_PEP_B, UNIPORT_B, STRUCTURES]
FASTA_IDX = [PEP_A, POS_IN_PEP_A, UNIPORT_A, PEP_B, POS_IN_PEP_B, UNIPORT_B]
FASTA_SUFFIX = ".fasta"
HEADLINE = "AvailableStructures\n"
UNIPORT_SERVER = "https://www.uniprot.org/uniprot/"


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


def create_fasta_dict(path):
    start = time.time()
    input_file = open(path)
    idx = np.array(FASTA_IDX)
    xl_fasta_dict = {}
    xl_samples = []
    i = 0
    for line in input_file:
        cols = np.array(line.split("\t"))
        if cols[STRUCTURES] != HEADLINE:
            xl_samples.append(np.asarray(cols[idx]))
            if cols[UNIPORT_A] not in xl_fasta_dict:
                fasta_a = get_fasta(cols[UNIPORT_A])
                xl_fasta_dict[cols[UNIPORT_A]] = fasta_a
            if cols[UNIPORT_B] not in xl_fasta_dict:
                fasta_b = get_fasta(cols[UNIPORT_B])
                xl_fasta_dict[cols[UNIPORT_B]] = fasta_b
        else:
            print(line)

    end = time.time()
    print("number of fastas: " + str(len(xl_fasta_dict)) + '\n')
    print("number of samples: " + str(len(xl_samples)) + '\n')
    print("time: " + str(end - start) + '\n')
    save_obj(xl_fasta_dict, "fasta_files_all_samples")
    save_obj(xl_samples, "all_samples")


def read_xl(path):
    input_file = open(path)
    idx = np.array(ALL_IDX)
    xl_with_struct = []
    for line in input_file:
        cols = np.array(line.split("\t"))
        if cols.size > STRUCTURES and cols[STRUCTURES] != "\n":
            xl_with_struct.append(np.asarray(cols[idx]))
    return xl_with_struct


def read_all_clear_dup(path):
    input_file = open(path)
    xl_dict = {}
    idx = np.array(ALL_IDX)
    xl_with_struct = []
    num_of_structures = 0
    unredundant_xl = 0
    for line in input_file:
        cols = np.array(line.split("\t"))
        if cols[STRUCTURES] != "\n" and cols[STRUCTURES] != "AvailableStructures\n":
            # if cols[PEP_A] == cols[PEP_B] and cols[POS_IN_PEP_A] == cols[POS_IN_PEP_B]:
            #     continue
            structures = cols[STRUCTURES][:-1]
            if cols[PEP_A] in xl_dict:
                if cols[PEP_B] in xl_dict[cols[PEP_A]]:
                    if structures in xl_dict[cols[PEP_A]]:
                        continue
                    xl_dict[cols[PEP_A]].add(structures)
                    xl_dict[cols[PEP_B]].add(structures)
                    xl_with_struct.append(np.asarray(cols[idx]))
                    xl_with_struct[-1][-1] = xl_with_struct[-1][-1][:-1]
                else:
                    unredundant_xl += 1
                    xl_with_struct.append(np.asarray(cols[idx]))
                    xl_with_struct[-1][-1] = xl_with_struct[-1][-1][:-1]
                    num_of_structures += len(structures.split(","))
                    if cols[PEP_B] not in xl_dict:
                        xl_dict[cols[PEP_B]] = {cols[PEP_A], structures}
                    else:
                        xl_dict[cols[PEP_B]].add(cols[PEP_A])
                        xl_dict[cols[PEP_B]].add(structures)
                    xl_dict[cols[PEP_A]].add(cols[PEP_B])
                    xl_dict[cols[PEP_A]].add(structures)
            else:
                unredundant_xl += 1
                xl_with_struct.append(np.asarray(cols[idx]))
                xl_with_struct[-1][-1] = xl_with_struct[-1][-1][:-1]
                num_of_structures += len(structures.split(","))
                xl_dict[cols[PEP_A]] = {cols[PEP_B], structures}
                xl_dict[cols[PEP_B]] = {cols[PEP_A], structures}

    avg_num_of_structures = num_of_structures / unredundant_xl
    save_obj(xl_with_struct, "xl_no_dup")
    return xl_with_struct


def main():
    input_path = sys.argv[1]
    # output_path = sys.argv[2]
    read_all_clear_dup(input_path)
    # create_fasta_dict(input_path)


if __name__ == "__main__":
    main()



