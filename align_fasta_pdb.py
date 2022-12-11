#
# import shutil
# from Bio import PDB
# from Bio.PDB.PDBParser import PDBParser
# from os.path import isfile
# from general_utils import *
# from alpha_fold_files import generate_af_model_name
# from alpha_fold_files import get_pred_align_err
#
# # PEP_A = 0
# # POS_IN_PEP_A = 1
# # UNI_A = 2
# # POS_RES_A = 3
# # PEP_B = 4
# # POS_IN_PEP_B = 5
# # UNI_B = 6
# # POS_RES_B = 7
# # DISTANCE = 8
# CONSERVED_AA = {'K', 'R', 'Q', 'E'}
# DIFFERENT_CHAINS = -1
# INVALID_SEQ = -1
# INVALID_XL_POS = -2
# MISSING_FASTA = -3
# SUCCESS = 1
# WORKDIR = "/cs/labs/dina/seanco/xl_parser/"
# # WORKDIR = 'C:/Users/seanco/Desktop/University/Master/Research/XL_DB_parser/'
# INTRA_MODELS_DIR = "/cs/labs/dina/seanco/xl_parser/best_models/"
# INTER_MODELS_DIR = "/cs/labs/dina/seanco/xl_parser/best_models/inter/"
#
#
# def do_write_double(uni_code_a, seq_a, pos_a, uni_code_b, seq_b, pos_b):
#     file_name = 'ali_files/' + uni_code_a + '_' + uni_code_b + '.ali'
#     if not isfile(file_name):
#         f_out = open(file_name, 'w')
#         f_out.write(">P1;" + uni_code_a + '_' + uni_code_b + '\n')
#         f_out.write("sequence:" + uni_code_a + '_' + uni_code_b + '::::::::' + '\n')
#         f_out.write(seq_a + '/' + '\n')
#         f_out.write(seq_b + '*' + '\n')
#         f_out.write('\n')
#         f_out.write('>P1;_fix_pos\n')
#         f_out.write('sequence:x: :: :: : :-1.00:-1.00\n')
#         arr = ['0' for i in range(len(seq_a))]
#         arr[int(pos_a)-1] = '4'
#         s = "".join(arr)
#         f_out.write(s + '/' + '\n')
#         arr = ['0' for i in range(len(seq_b))]
#         arr[int(pos_b)-1] = '4'
#         s = "".join(arr)
#         f_out.write(s + '*' + '\n')
#         f_out.close()
#
#
# def write_double_ali(uni_code_a, seq_a, pos_a, uni_code_b, seq_b, pos_b, errors_file):
#     if isinstance(seq_a, int):
#         errors_file.write(uni_code_a + '\n')
#         return MISSING_FASTA
#     if isinstance(seq_b, int):
#         errors_file.write(uni_code_b + '\n')
#         return MISSING_FASTA
#     if len(seq_a) < 1:
#         errors_file.write(uni_code_a + '\n')
#         return INVALID_SEQ
#     if len(seq_b) < 1:
#         errors_file.write(uni_code_b + '\n')
#         return INVALID_SEQ
#     if int(pos_a) >= len(seq_a) or int(pos_b) >= len(seq_b):
#         errors_file.write(uni_code_a + '_' + uni_code_b + '_' + str(pos_a) + '_' + str(pos_b) + '\n')
#         return INVALID_XL_POS
#     do_write_double(uni_code_a, seq_a, pos_a, uni_code_b, seq_b, pos_b)
#     do_write_double(uni_code_b, seq_b, pos_b, uni_code_a, seq_a, pos_a)
#     return SUCCESS
#
#
# def write_single_ali(uni_code, seq, pos_a, pos_b, errors_file):
#     if isinstance(seq, int):
#         errors_file.write(uni_code + '\n')
#         return MISSING_FASTA
#     if len(seq) < 1:
#         errors_file.write(uni_code + '\n')
#         return INVALID_SEQ
#     if int(pos_a) > len(seq) or int(pos_b) > len(seq):
#         errors_file.write(uni_code + '_' + str(pos_a) + '_' + str(pos_b) + '\n')
#         return INVALID_XL_POS
#     file_name = 'ali_files/' + uni_code + '.ali'
#     if not isfile(file_name):
#         f_out = open(file_name, 'w')
#         f_out.write(">P1;" + uni_code + '\n')
#         f_out.write("sequence:" + uni_code + '::::::::' + '\n')
#         f_out.write(seq + '*' + '\n')
#         f_out.write('\n')
#         f_out.write('>P1;_fix_pos\n')
#         f_out.write('sequence:x: :: :: : :-1.00:-1.00\n')
#         arr = ['0' for i in range(len(seq))]
#         arr[int(pos_a)-1] = '4'
#         if pos_b != DIFFERENT_CHAINS:
#             arr[int(pos_b)-1] = '4'
#         s = "".join(arr)
#         f_out.write(s + '*' + '\n')
#         f_out.close()
#     return SUCCESS
#
#
# def write_ali_files(fasta_dict, xl_list):
#     f = open('bad_samples.txt', 'a')
#     already_created = set()
#     errors = 0
#     missing_fasta = 0
#     for sample in xl_list:
#         if sample[XL_UNIPORT_A] == sample[XL_UNIPORT_B]:
#             key = sample[XL_UNIPORT_A]
#         else:
#             key = sample[XL_UNIPORT_A] + '.' + sample[XL_UNIPORT_B]
#         if key not in already_created:
#             if sample[XL_UNIPORT_A] != sample[XL_UNIPORT_B]:
#                 res = write_double_ali(sample[XL_UNIPORT_A], fasta_dict[sample[XL_UNIPORT_A]],
#                                        sample[XL_RES_NUM_A], sample[XL_UNIPORT_B],
#                                        fasta_dict[sample[XL_UNIPORT_B]], sample[XL_RES_NUM_B], f)
#                 if res != INVALID_XL_POS:
#                     already_created.add(key)
#                     key = sample[XL_UNIPORT_B] + '.' + sample[XL_UNIPORT_A]
#                     already_created.add(key)
#             else:
#                 res = write_single_ali(sample[XL_UNIPORT_A], fasta_dict[sample[XL_UNIPORT_A]],
#                                        sample[XL_RES_NUM_A], sample[XL_RES_NUM_B], f)
#                 if res != INVALID_XL_POS:
#                     already_created.add(key)
#             if res == INVALID_XL_POS or res == INVALID_SEQ:
#                 errors += 1
#             elif res == MISSING_FASTA:
#                 missing_fasta += 1
#
#     print("errors with xl:" + str(errors))
#     print("missing fasta: " + str(missing_fasta))
#     f.close()
#
#
# def generate_model_name_low(pdb, unicode_a, unicode_b=None, models_dir='best_models/', pref='pdb'):
#     pdb_name = pref + pdb.lower()
#     name = pdb_name + '_' + unicode_a
#     if unicode_b is not None:
#         name += '_' + unicode_b
#     name += '.B99990001.pdb'
#     return models_dir + name
#
#
# def generate_model_name(sample, af_set, pdb, cif_files, file_type=REGULAR_PDBS):
#     file_name = None
#     if file_type == REGULAR_PDBS or file_type == KNOWN_STRUCTURE_PDB:
#         pref = 'pdb'
#         if pdb in cif_files:
#             pref = ''
#         if sample[XL_UNIPORT_A] != sample[XL_UNIPORT_B]:
#             file_name = generate_model_name_low(pdb, sample[XL_UNIPORT_A], sample[XL_UNIPORT_B],
#                                                 INTER_MODELS_DIR, pref)
#             if not isfile(file_name):
#                 file_name = generate_model_name_low(pdb, sample[XL_UNIPORT_B], sample[XL_UNIPORT_A],
#                                                     INTER_MODELS_DIR, pref)
#         else:
#             file_name = generate_model_name_low(pdb, sample[XL_UNIPORT_A], None, INTRA_MODELS_DIR)
#     elif file_type == AF_PDB:
#         file_name = generate_af_model_name(pdb)
#
#     if not isfile(file_name):
#         return None
#     return file_name
#
#
# def calc_dist_single_xl(sample, pdb_chains):
#     target_atoms = 'CA'
#     atom_ca_a = 0
#     atom_ca_b = 0
#     dist = 0
#     chain_a, chain_b = 0, 0
#     all_seq = ''
#     for chain in pdb_chains:
#         all_seq += chain.get_sequence()
#
#     idx_pep_a = all_seq.find(sample[XL_PEP_A])
#     if idx_pep_a != -1:
#         lys_idx_a = idx_pep_a + int(sample[XL_POS_IN_PEP_A])
#         idx_pep_b = all_seq.find(sample[XL_PEP_B])
#
#         if idx_pep_b != -1:
#             lys_idx_b = idx_pep_b + int(sample[XL_POS_IN_PEP_B])
#             cur_idx = 0
#             i = 0
#             while i < len(pdb_chains) and (cur_idx <= lys_idx_a or cur_idx <= lys_idx_b):
#                 chain = pdb_chains[i]
#                 if cur_idx <= lys_idx_a < cur_idx + len(chain):
#                     res_a = chain[lys_idx_a - cur_idx]
#                     if target_atoms in res_a:
#                         atom_ca_a = res_a[target_atoms]
#                     chain_a = i
#                 if cur_idx <= lys_idx_b < cur_idx + len(chain):
#                     res_b = chain[lys_idx_b - cur_idx]
#                     if target_atoms in res_b:
#                         atom_ca_b = res_b[target_atoms]
#                     chain_b = i
#                 cur_idx += len(chain)
#                 i += 1
#
#             if isinstance(atom_ca_a, int) or isinstance(atom_ca_b, int):
#                 return [0, 0, 0, 0, 0]
#             dist = atom_ca_a - atom_ca_b
#             if dist > 60:
#                 dist = 0
#             return [chain_a, res_a.id[1], chain_b, res_b.id[1], dist]
#
#     return [0, 0, 0, 0, 0]
#
#
# def print_tests_results(tot_samples, tot_dif, dif_categories, errors):
#     lows = {0: 0, 1: 5, 2: 10, 3: 20}
#     highs = {0: 5, 1: 10, 2: 20, 3: 60}
#     res_file = open(WORKDIR + "distances_dif.txt", 'w')
#     avg = tot_samples / tot_dif
#     res_file.write("average diff in distances: " + str(avg) + "\n")
#     for i in range(len(dif_categories)):
#         res_file.write("samples with dif of " + str(lows[i]) + " to " + str(highs[i]) + " : " + str(dif_categories[i]) + "\n")
#     res_file.write("Exception raised: " + str(errors))
#     res_file.close()
#
#
# def check_valid_sample_dist(dist):
#     try:
#         x = float(dist)
#         if x == 0 or x > 60:
#             return False
#         return x
#     except ValueError:
#         return False
#
#
# def categories_dif_distance(dif, categories):
#     if dif <= 5:
#         categories[0] += 1
#     elif dif <= 10:
#         categories[1] += 1
#     elif dif <= 20:
#         categories[2] += 1
#     else:
#         categories[3] += 1
#
#
# def validate_model_xl_distance_af(model_dist, sample, pp_list):
#     n_res = len(pp_list[0])
#     pred_error = get_pred_align_err(sample[XL_ACCESSION_A], sample[XL_RES_NUM_A], sample[XL_RES_NUM_B], n_res)
#     return pred_error
#
#
# def validate_model_xl_distance_regular(model_dist, sample):
#     sample_dist = check_valid_sample_dist(sample[XL_DISTANCE])
#     if sample_dist is not False:
#         dif = abs(model_dist - float(sample[XL_DISTANCE]))
#         return dif
#     return -1
#
#
# def validate_model_xl_distance(sample, pdb_type, pp_list, linker_types):
#     model_xl_data = calc_dist_single_xl(sample, pp_list)
#     model_dist = model_xl_data[-1]
#     if model_dist > 0:
#         if pdb_type == AF_PDB:
#             pred_or_diff = validate_model_xl_distance_af(model_dist, sample, pp_list)
#         else:
#             pred_or_diff = validate_model_xl_distance_regular(model_dist, sample)
#         model_xl_data.append(pdb_type)
#         model_xl_data.append(pred_or_diff)
#         if sample[XL_DATASETS] in linker_types:
#             linker = linker_type_dict[linker_types[sample[XL_DATASETS]]]
#         else:
#             linker = 0
#         model_xl_data.append(linker)
#     return model_xl_data, model_dist
#
#
# def filter_files(file_name, xl_data):
#     if xl_data[5] != AF_PDB:
#         is_accurate = xl_data[6] < 10
#     else:
#         is_accurate = xl_data[6] <= 15
#     if is_accurate and xl_data[-1] == DSSO_LINKER:
#         if xl_data[4] <= 35:
#             shutil.copy(file_name, WORKDIR + 'best_models/dsso/')
#     elif is_accurate and xl_data[-1] == BDP_NHP_LINKER:
#         if xl_data[4] <= 40:
#             shutil.copy(file_name, WORKDIR + 'best_models/bdp_nhp/')
#
#
# def test_models(xl_list):
#     af_set = load_obj('af_pdbs_unicodes')
#     cif_files = load_obj('cif_files')
#     linker_types = load_obj('linker_types')
#     good_models = {}
#     pdb_parser = PDBParser(PERMISSIVE=1)
#     ppb = PDB.CaPPBuilder()
#     prev_sample = [""] * 34
#     prev_pdb = ""
#     parser_exceptions = 0
#     no_file_error = 0
#     pdb_file_name = None
#     counter = 0
#     for sample in xl_list:
#         if sample[XL_UNIPORT_A] != sample[XL_UNIPORT_B]:
#             continue
#         pdb_names, pdb_type = pdb_from_sample(sample, af_set)
#         if len(pdb_names) > 1:
#             pdb_names = [pdb_names[0]]
#         for cur_pdb_name in pdb_names:
#             if sample[XL_UNIPORT_A] != prev_sample[XL_UNIPORT_A] or sample[XL_UNIPORT_B] != prev_sample[XL_UNIPORT_B] \
#                     or cur_pdb_name != prev_pdb:
#                 prev_pdb = cur_pdb_name
#                 prev_sample = sample
#                 pdb_file_name = generate_model_name(sample, af_set, cur_pdb_name, cif_files, pdb_type)
#                 if pdb_file_name is None:
#                     no_file_error += 1
#                     continue
#                 try:
#                     structure = pdb_parser.get_structure(cur_pdb_name, pdb_file_name)
#                     polypeptide_list = ppb.build_peptides(structure, 1)
#                 except:
#                     parser_exceptions += 1
#                     continue
#             if pdb_file_name is not None:
#                 model_xl_data, model_dist = validate_model_xl_distance(sample, pdb_type, polypeptide_list,
#                                                                        linker_types)
#                 if model_dist > 0:
#                     key = pdb_file_name + '+' + sample[XL_PEP_A] + '_' + sample[XL_PEP_B] + '_' +  \
#                           str(model_xl_data[1]) + '_' + str(model_xl_data[3])
#                     good_models[key] = model_xl_data
#                     # filter_files(pdb_file_name, model_xl_data)
#                     counter += 1
#                     if counter > 1000:
#                         print(len(good_models))
#                         save_obj(good_models, "all_models_dict")
#                         counter = 0
#
#     save_obj(good_models, "all_models_dict")
#     # print("missed files: ", no_file_error)
#     print(len(good_models))
#
#
# def test_data(xl_list):
#     counter = 0
#     count_no_dup = 0
#     already_created = set()
#     intra = 0
#     for sample in xl_list:
#         pdbs = sample[-1].split(',')
#         counter += len(pdbs)
#         if sample[XL_UNIPORT_A] == sample[XL_UNIPORT_B]:
#             intra += 1
#             key = sample[XL_UNIPORT_A]
#         else:
#             key = sample[XL_UNIPORT_A] + '.' + sample[XL_UNIPORT_B]
#         if key not in already_created:
#             if sample[XL_UNIPORT_A] != sample[XL_UNIPORT_B]:
#                 count_no_dup += len(pdbs)
#                 already_created.add(key)
#                 key = sample[XL_UNIPORT_B] + '.' + sample[XL_UNIPORT_A]
#                 already_created.add(key)
#             else:
#                 count_no_dup += len(pdbs)
#                 already_created.add(key)
#     print("total pdb available: " + str(counter) + "\n")
#     print("no duplications available pdbs: " + str(count_no_dup))
#
#
# def main():
#     # fasta_dict = general_utils.load_obj(sys.argv[1])
#     # xl_list = general_utils.load_obj(sys.argv[2])
#     # xl_list = load_obj("xl_with_dup_all_cols")
#     xl_list = load_obj("xl_with_dup_all_cols_af_organs")
#     # get_pdb_files(xl_list)
#     # download_missing_pdbs()
#     # fasta_dict = load_obj('fasta_files_dict')
#     # write_ali_files(fasta_dict, xl_list)
#     # test_data(xl_list)
#     test_models(xl_list)
#     # align(fasta_dict, xl_list)
#
#
# if __name__ == "__main__":
#     main()
