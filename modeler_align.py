import sys

from modeller import *

import general_utils

ALI_PATH = './ali_files/'
ALI_SUFF = '.ali'
PAP_SUFF = '.pap'
ALI_OUT_PATH = './alignments/ali/'
PAP_OUT_PATH = './alignments/pap/'
UNI_A = 2
UNI_B = 6
POS_RES_A = 3
POS_RES_B = 7

log.verbose()
env = Environ()
env.io.atom_files_directory = ['./pdbs']
env.io.hetatm = True


def align_single(pdb_name, uni_code, error_file):
    try:
        mdl = model(env, file=pdb_name)
        aln = Alignment(env)
        aln.append_model(mdl, atom_files=pdb_name, align_codes=pdb_name)
        aln.append(file=ALI_PATH + uni_code + ALI_SUFF, align_codes=(uni_code, '_fix_pos'))
        aln.salign(rr_file='$(LIB)/as1.sim.mat', # Substitution matrix used
                   fix_offsets=(0, -10, -20, -30, -40),
                   output='', max_gap_length=10,
                   gap_function=True,  # to use structure-dependent gap penalty
                   alignment_type='PAIRWISE',
                   # feature_weights=(1., 0., 0., 0., 0., 0.), overhang=99,
                   feature_weights=(1., 0., 0., 0., 0., 0.), overhang=99,
                   gap_penalties_1d=(-450, 0),
                   gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
                   similarity_flag=True,
                   local_alignment=True)

        out_name = ALI_OUT_PATH + uni_code + '_' + pdb_name + ALI_SUFF
        aln.write(file=out_name, alignment_format='PIR')
        out_name = PAP_OUT_PATH + uni_code + '_' + pdb_name + PAP_SUFF
        aln.write(file=out_name, alignment_format='PAP')
    except OSError:
        error_file.write(pdb_name + '\n')


def modeler_align(xl_list, bad_samples):
    error_file = open('missing_pdb_files', 'w')
    already_aligned = set()
    for sample in xl_list:
        if sample[UNI_A] in bad_samples or sample[UNI_B] in bad_samples:
            continue
        pdb_names = sample[-1].split(',')
        for name in pdb_names:
            key_name = sample[UNI_A] + name
            if key_name not in already_aligned:
                align_single('pdb'+name, sample[UNI_A], error_file)
                already_aligned.add(key_name)
            key_name = sample[UNI_B] + name
            if key_name not in already_aligned:
                align_single('pdb'+name, sample[UNI_B], error_file)
                already_aligned.add(key_name)

    error_file.close()


def main():
    xl_list = general_utils.load_obj(sys.argv[1])
    bad_samples_file = open(sys.argv[2])
    bad_samples = set()
    for line in bad_samples_file:
        bad_samples.add(line[:-1])
    modeler_align(xl_list, bad_samples)


if __name__ == "__main__":
    main()
