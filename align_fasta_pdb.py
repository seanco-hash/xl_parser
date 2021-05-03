import sys

from Bio.SubsMat.MatrixInfo import blosum62

import general_utils
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio import Align
from os import listdir
from os.path import isfile, join

# TODO: Update to real indexes, update the xl reader to include residue positions
UNI_A = 2
UNI_B = 6
POS_RES_A = 3
POS_RES_B = 7
CONSERVED_AA = {'K', 'R', 'Q', 'E'}


# TODO: run to download missing pdb files
def check_missing_files(xl_list):
    downloaded = set()
    missing_xl = []
    for f in listdir("./pdbs"):
        if isfile(join("./pdbs", f)):
            downloaded.add(f)

    for sample in xl_list:
        pdb_names = sample[-1].split(',')
        for name in pdb_names:
            file_name = "pdb" + name.lower() + ".ent"
            if file_name not in downloaded:
                missing_xl.append(name)

    get_pdb_files(missing_xl, "mmCif")


def get_pdb_files(xl_list, wanted_format="pdb"):
    pdb1 = PDB.PDBList()
    already_downloaded = set()
    for sample in xl_list:
        pdb_names = sample[-1].split(',')
        for name in pdb_names:
            if name not in already_downloaded:
                pdb1.retrieve_pdb_file(name, pdir="./pdbs", overwrite=False, file_format=wanted_format)
                already_downloaded.add(name)


def align_binary_search(alignment, pos, start, end):
    if end <= start:
        return -1
    half = (start + end) / 2
    if pos < alignment[half][0]:
        return align_binary_search(alignment, pos, start, half)
    if pos > alignment[half][1]:
        return align_binary_search(alignment, pos, half+1, end)
    return half


def check_align(alignments, pos, target_seq):
    for alignment in alignments:
        fasta_align = alignment.aligned[0]
        seq_align = alignment.aligned[1]
        idx = align_binary_search(fasta_align, pos, 0, len(fasta_align))
        map_align = seq_align[idx]
        shift = pos - fasta_align[idx][0]
        new_residue_pos = map_align[0] + shift
        if target_seq[new_residue_pos] in CONSERVED_AA:
            return alignment

    return False


def do_align(fasta_seq, polypep_list, pos):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = blosum62
    best_score = 0
    best_chain_num = 0
    best_align = 0
    for i, pep in enumerate(polypep_list):
        pep_seq = pep.get_sequence()
        alignments = aligner.align(fasta_seq, pep_seq)
        if alignments.score > best_score:
            res_alignment = check_align(alignments, pos, pep_seq)
            if res_alignment is not False:
                best_align = res_alignment
                best_chain_num = i
                best_score = alignments.score
    return best_align, best_chain_num


def align(fasta_dict, xl_list):
    pdb_parser = PDBParser(PERMISSIVE=1)
    ppb = PDB.CaPPBuilder
    for sample in xl_list:
        pdb_names = sample[-1].split(',')
        for name in pdb_names:
            structure = pdb_parser.get_structure(name, "pdbs/pdb" + name + ".ent")
            polypeptide_list = ppb.build_peptides(structure, 1)
            pep_a_align = do_align(fasta_dict[sample[UNI_A]], polypeptide_list, sample[POS_RES_A])
            pep_b_align = do_align(fasta_dict[sample[UNI_B]], polypeptide_list, sample[POS_RES_B])


def main():
    fasta_dict = general_utils.load_obj(sys.argv[1])
    xl_list = general_utils.load_obj(sys.argv[2])
    get_pdb_files(xl_list)


if __name__ == "__main__":
    main()