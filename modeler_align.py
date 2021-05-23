import os
import shutil
import sys

from modeller import *
from modeller import soap_loop
from modeller.automodel import automodel, assess
from modeller.scripts import complete_pdb

import general_utils

AUTO_PDB_SUFF = '.B999900'
ALI_PATH = './ali_files_try/'
ALI_SUFF = '.ali'
PAP_SUFF = '.pap'
ALI_OUT_PATH = './alignments/ali/'
PAP_OUT_PATH = './alignments/pap/'
UNI_A = 2
UNI_B = 6
POS_RES_A = 3
POS_RES_B = 7
NUM_OF_MODELS_TO_CREATE = 1

log.verbose()
env = Environ()
env.io.atom_files_directory = ['./pdbs']
env.io.hetatm = True


def remove_files(out_name, n):
    file_to_remove = out_name + '.V999900%02d' % n
    os.remove(file_to_remove)
    file_to_remove = out_name + '.D000000%02d' % n
    os.remove(file_to_remove)


def create_models(pdb_name, unicode, aln):
    # sp = soap_loop.Scorer()
    os.chdir('./tmp_models')
    out_name = unicode + '_' + pdb_name
    a_model = automodel(env, alnfile=aln, knowns=pdb_name, sequence=unicode,
                        assess_methods=assess.DOPE, root_name=out_name)
    a_model.starting_model = 1
    a_model.ending_model = NUM_OF_MODELS_TO_CREATE
    a_model.make()

    best_score = 0
    best_model = "x"
    for n in range(1, NUM_OF_MODELS_TO_CREATE + 1):
        model_name = out_name + AUTO_PDB_SUFF + '%02d.pdb' % n
        mdl = complete_pdb(env, model_name)
        s = selection(mdl)  # all atom selection
        score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',  # file='TCR.B9999000' + str(n) + '.profile',
                         normalize_profile=True, smoothing_window=15)
        remove_files(out_name, n)
        if score < best_score:
            best_score = score
            best_model = model_name

    shutil.copy(best_model, '../best_models/')
    os.chdir('..')
    return best_score, best_model
    # print(best_model, best_score)


def align_single(pdb_name, error_file, uni_code_a, uni_code_b=None, is_together=False):
    try:
        mdl = model(env, file=pdb_name)
        aln = Alignment(env)
        aln.append_model(mdl, atom_files=pdb_name, align_codes=pdb_name)
        if is_together:
            unicode = uni_code_a + '_' + uni_code_b
        else:
            unicode = uni_code_a
        aln.append(file=ALI_PATH + unicode + ALI_SUFF, align_codes=(unicode, '_fix_pos'))
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

        ali_out_name = ALI_OUT_PATH + unicode + '_' + pdb_name + ALI_SUFF
        aln.write(file=ali_out_name, alignment_format='PIR')
        pap_out_name = PAP_OUT_PATH + unicode + '_' + pdb_name + PAP_SUFF
        aln.write(file=pap_out_name, alignment_format='PAP')
        score, _model = create_models(pdb_name, unicode, aln)
        return score, _model

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
            if sample[UNI_A] != sample[UNI_B]:
                tmp_name = sample[UNI_A] + '.' + sample[UNI_B]
                key_name = tmp_name + name
                key_name2 = sample[UNI_B] + '.' + sample[UNI_A] + name
                if key_name not in already_aligned and key_name2 not in already_aligned:
                    score1, model1 = align_single('pdb' + name, error_file, sample[UNI_A],
                                                  sample[UNI_B], True)
                    score2, model2 = align_single('pdb' + name, error_file, sample[UNI_B],
                                                  sample[UNI_A], True)
                    if score1 < score2:
                        to_del = model2
                    else:
                        to_del = model1
                    os.remove('./best_models/' + to_del)
                    already_aligned.add(key_name)
                    already_aligned.add(key_name2)
            else:
                key_name = sample[UNI_A] + name
                if key_name not in already_aligned:
                    align_single('pdb'+name, error_file, sample[UNI_A])
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
