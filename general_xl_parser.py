import numpy as np
import general_utils
import yaml
import cross_link
import data_table_parser


def load_config(cfg_file):
    with open(cfg_file, "r") as f:
        cfg = yaml.safe_load(f)
        return cfg


class XlParser:
    def __init__(self, cfg):
        self.deliminator = cfg['DELIMINATOR']
        self.parse_uniport_needed = cfg['IS_UNIPORT_NEEDS_PARSE']

    @staticmethod
    def parse_single_sample_uniport(cfg, full_uniport):
        if cfg['OUT_FILE_NAME'] == "max_linker_xl":
            return XlParser.parse_maxlinker_uniport(full_uniport)
        pattern = cfg['UNIPORT_PARSE_PATTERN']
        i = 0
        while i < len(pattern):
            slices = full_uniport.split(str(pattern[i]))
            i += 1
            full_uniport = slices[int(pattern[i])]
            i += 1
        return full_uniport

    def parse_uniport(self, full_uniport_a, full_uniport_b, cfg):
        if self.parse_uniport_needed:
            uniport_a = XlParser.parse_single_sample_uniport(cfg, full_uniport_a)
            if full_uniport_a == full_uniport_b:
                return uniport_a, uniport_a
            else:
                uniport_b = XlParser.parse_single_sample_uniport(cfg, full_uniport_b)
                return uniport_a, uniport_b
        return full_uniport_a, full_uniport_b

    @staticmethod
    def parse_maxlinker_uniport(line):
        uniports = line.split(';')
        sliced_uniports = [item.split('-')[0] for item in uniports]
        i = 1
        while i < len(sliced_uniports):
            if sliced_uniports[0] != sliced_uniports[i]:
                return ""
            i += 1
        return sliced_uniports[0]

    @staticmethod
    def parse_maxlinker_res_num(line):
        candidates = line.split(';')
        for candidate in candidates:
            unipot_pos = candidate.split('_')
            only_uniport = unipot_pos[0].split('-')
            if len(only_uniport) == 1:
                return unipot_pos[1]
        return '-1'

    @staticmethod
    def parse_res_num(cols, cfg):
        if cfg['RES_A_POS'] >= 0:
            if not cfg['IS_RES_POS_NEEDS_PARSE']:
                res_num_a = cols[cfg['RES_A_POS']]
                res_num_b = cols[cfg['RES_B_POS']]
            elif cfg['OUT_FILE_NAME'] == "max_linker_xl":
                res_num_a = XlParser.parse_maxlinker_res_num(cols[cfg['RES_A_POS']])
                res_num_b = XlParser.parse_maxlinker_res_num(cols[cfg['RES_B_POS']])
            else:
                return '-1', '-1'
        else:
            res_num_a = int(cols[cfg['PEP_A_POS']]) + int(cols[cfg['POS_IN_PEP_A']]) - 1
            res_num_b = int(cols[cfg['PEP_B_POS']]) + int(cols[cfg['POS_IN_PEP_B']]) - 1
        return res_num_a, res_num_b

    @staticmethod
    def parse_maxlinker_pos_in_pep(line):
        i = 0
        while i < len(line):
            if str.islower(line[i]):
                return i
            i += 1
        return -1

    @staticmethod
    def parse_pos_in_pep(line, cfg):
        if cfg['POS_IN_PEP_A'] == -1:
            return -1, -1
        if cfg['OUT_FILE_NAME'] == "max_linker_xl":
            a = XlParser.parse_maxlinker_pos_in_pep(line[cfg['PEP_A']])
            b = XlParser.parse_maxlinker_pos_in_pep(line[cfg['PEP_B']])
            return a, b
        return line[cfg['POS_IN_PEP_A']], line[cfg['POS_IN_PEP_B']]
    #
    # @staticmethod
    # def parse_peps_klykov(line, cfg):
    #     pep_a = line[cfg['PEP_A']].replace('[', '')
    #     pep_b = line[cfg['PEP_B']].replace('[', '')
    #     pep_a = pep_a.replace(']', '')
    #     pep_b = pep_b.replace(']', '')
    #     return pep_a, pep_b

    # @staticmethod
    # def parse_peptides(line, cfg):
    #     if cfg['IS_PEP_NEEDS_PARSE']:
    #         if cfg['OUT_FILE_NAME'] == "klykov_et_al_20":
    #             return XlParser.parse_peps_klykov(line, cfg)
    #     return line[cfg['PEP_A']], line[cfg['PEP_B']]

    def parse(self, cfg, xl_file):
        xl_dict = {}
        cross_links = []
        with open(xl_file, 'r') as f:
            first_line = f.readline()
            for line in f:
                cols = np.array(line.split(self.deliminator))
                uniport_a, uniport_b = self.parse_uniport(cols[cfg['UNIPORT_A']], cols[cfg['UNIPORT_B']], cfg)
                if uniport_a == "" or uniport_b == "":
                    continue
                res_num_a, res_num_b = XlParser.parse_res_num(cols, cfg)
                if res_num_a == -1 or res_num_b == -1:
                    continue
                if cfg['PEP_A'] != -1 and cfg['POS_IN_PEP_A'] != -1:
                    key_a = uniport_a + '+' + cols[cfg['PEP_A']] + '+' + cols[cfg['POS_IN_PEP_A']]
                    key_b = uniport_b + '+' + cols[cfg['PEP_B']] + '+' + cols[cfg['POS_IN_PEP_B']]
                else:
                    key_a = uniport_a + '+' + res_num_a
                    key_b = uniport_b + '+' + res_num_b
                if key_a in xl_dict and key_b in xl_dict[key_a]:
                    continue
                pos_in_pep_a, pos_in_pep_b = XlParser.parse_pos_in_pep(cols, cfg)
                if cfg['PEP_A'] != -1:
                    xl_obj = cross_link.CrossLink(cols[cfg['PEP_A']], pos_in_pep_a, res_num_a,
                                       uniport_a, cols[cfg['PEP_B']],
                                       pos_in_pep_b, res_num_b, uniport_b,
                                       cfg['LINKER_TYPE'], cfg['PDB_NAME'], cfg['OUT_FILE_NAME'])
                else:
                    xl_obj = cross_link.CrossLink("", cfg['POS_IN_PEP_A'], res_num_a,
                                       uniport_a, "",
                                       cfg['POS_IN_PEP_B'], res_num_b, uniport_b,
                                       cfg['LINKER_TYPE'], cfg['PDB_NAME'], cfg['OUT_FILE_NAME'])
                cross_links.append(xl_obj)
                if key_a in xl_dict:
                    if key_b not in xl_dict:
                        xl_dict[key_b] = {key_a}
                    else:
                        xl_dict[key_b].add(key_a)
                    xl_dict[key_a].add(key_b)
                else:
                    xl_dict[key_a] = {key_b}
                    xl_dict[key_b] = {key_a}
        print("found cross links: ", len(cross_links))
        return cross_links

    @staticmethod
    def save_cross_links(cross_links, cfg=None, name="xl_objects_xldb_data", dir_="xl_objects/"):
        if cfg is not None:
            file_name = dir_ + "xl_objects_" + cfg['OUT_FILE_NAME']
        else:
            file_name = dir_ + name
        general_utils.save_obj(cross_links, file_name)

    @staticmethod
    def convert_old_xl_list_to_cross_link_objects():
        xl_list = general_utils.load_obj('xl_with_dup_all_cols')
        linker_types = general_utils.load_obj('linker_types')
        object_list = []
        already_created = set()
        for sample in xl_list:
            pdb_file = ''
            if sample[general_utils.XL_UNIPORT_A] == sample[general_utils.XL_UNIPORT_B]:
                if sample[general_utils.XL_PDB_TYPE_A] == general_utils.XL_AVAILABLE_IN_PDB:
                    pdb_file = sample[general_utils.XL_PDB_A]
            else:
                if sample[general_utils.XL_STRUCTURES] != '':
                    pdb_file = sample[general_utils.XL_STRUCTURES].split(',')[-1]
            if sample[general_utils.XL_DATASETS] in linker_types:
                linker = linker_types[sample[general_utils.XL_DATASETS]]
            else:
                linker = ''
            cl = cross_link.CrossLink(sample[general_utils.XL_PEP_A], sample[general_utils.XL_POS_IN_PEP_A],
                           sample[general_utils.XL_RES_NUM_A], sample[general_utils.XL_UNIPORT_A],
                           sample[general_utils.XL_PEP_B], sample[general_utils.XL_POS_IN_PEP_B],
                           sample[general_utils.XL_RES_NUM_B], sample[general_utils.XL_UNIPORT_B],
                           linker, pdb_file, sample[general_utils.XL_DATASETS])
            key1, key2 = cl.get_keys()
            if not (key1 in already_created or key2 in already_created):
                object_list.append(cl)
                already_created.add(key1)

        XlParser.save_cross_links(object_list)


# def main():
#     cfg_path = '/cs/labs/dina/seanco/xl_parser/configs/parser_max_linker.yaml'
#     cfg = load_config(cfg_path)
#     p = XlParser(cfg)
#     xl_obj_list = p.parse(cfg, cfg['DATA_FILE'])
#     XlParser.save_cross_links(xl_obj_list, cfg)
#     # XlParser.convert_old_xl_list_to_cross_link_objects()
#
#
# if __name__ == "__main__":
#     main()
