import os
import sys
import subprocess
import os.path as osp
from argparse import ArgumentParser


WORKDIR = "/cs/labs/dina/seanco/xl_parser/"

# trans_file = "sorted_pos"
structdir = "dockground_1"

out_dir = "/cs/labs/dina/seanco/scripts"

ask_time = "1-0"
mem = "64000M"
mem_cpu = "64000M"
cpus = "6"
cpus_only = "2"
gpus = "1"
gpu_mem = "16g"

INTRO = "#!/bin/tcsh \n#SBATCH --mem=" + mem + "\n#SBATCH -c" + cpus + " \n#SBATCH --time=" + ask_time + \
        "\n#SBATCH --killable\n" + \
      "cd /cs/labs/dina/seanco/xl_parser\n"

PYTHON = "python3 "

ALPHAFOLD_SCRIPT = " \"/cs/labs/dina/bshor/local_alphafold/seanco_run_alphafold_sequence.sh\""
COLABFOLD_RUN = "/cs/labs/dina/seanco/xl_parser/alphafold/colabfold/run_alphafold_from_uniport_list.py"
CROSS_LINK_RUN = "/cs/labs/dina/seanco/xl_parser/cross_link.py"
COLABFOLD_BATCH_RUN = "colabfold_batch "
COLABFOLD_SCRIPT = "/cs/labs/dina/seanco/xl_parser/alphafold/run_colabfold_script.sh"
CROSS_LINK_SCRIPT = "/cs/labs/dina/seanco/xl_parser/scripts/run_cross_link.sh"
COLABFOLD_BATCH_SCRIPT = "/cs/labs/dina/seanco/xl_parser/alphafold/run_colabfold_batch_script.sh"
VERSION_CHECK_SCRIPT = "/cs/labs/dina/seanco/xl_parser/scripts/version_check.sh"
INLINE_INTRO = "--mem=" + mem + " -c" + cpus + " --time=" + ask_time + " --gres=gpu:" + gpus + "," \
               "vmem:" + gpu_mem + " "
INLINE_INTRO_CPU = "--mem=" + mem_cpu + " -c" + cpus_only + " --exclude=sm-[01-04,07-08],sm-[17-20],sm-16,sm-15,ohm-55,ohm-61,ohm-62,ohm-64,ohm-57,ohm-58,ohm-59,ohm-60,cb-06,cb-05,cb-11,ohm-63,cb-16,cb-08,cb-13,cb-14,cb-17" + " --time=" + ask_time + " "


def line_gen(lst, factor):
    for i in range(0,len(lst),factor):
        yield lst[i:i+factor]


def veifyDir(path):
    if not osp.isdir(path):
        os.mkdir(path)


def run_colabfold():
    with open(COLABFOLD_SCRIPT, 'w') as script:
        script.write("#!/bin/tcsh\n")
        script.write("cd /cs/labs/dina/seanco/xl_parser/alphafold/colabfold\n")
        # script.write("module load cuda")
        script.write("module load cuda/11.0\n")
        script.write("module load cudnn/8.0.2\n")
        # script.write("module load cudnn\n")
        # script.write("module load tensorflow/2.5.0\n")
        script.write(PYTHON + COLABFOLD_RUN + "\n")
        script.close()
    subprocess.run("sbatch " + INLINE_INTRO + COLABFOLD_SCRIPT, shell=True)


def run_cross_link():
    with open(CROSS_LINK_SCRIPT, 'w') as script:
        script.write("#!/bin/tcsh\n")
        script.write("cd /cs/labs/dina/seanco/xl_parser/\n")
        script.write("source /cs/labs/dina/seanco/xl_parser/xl_db_parser_venv/bin/activate.csh\n")
        script.write(PYTHON + CROSS_LINK_RUN + "\n")
        script.close()
    subprocess.run("sbatch " + INLINE_INTRO_CPU + CROSS_LINK_SCRIPT, shell=True)


def check_versions():
    subprocess.run("sbatch " + INLINE_INTRO + VERSION_CHECK_SCRIPT, shell=True)


def run_colabfold_batch(input_file_or_dir, output_dir, num_recycle=3, gpu_memory='16g', _ask_time='1:0:0'):
    sbatch_args = "--mem=" + mem + " -c" + cpus + " --time=" + _ask_time + " --gres=gpu:" + gpus + \
                  ",vmem:" + gpu_memory + " "
    print(sbatch_args)
    with open(COLABFOLD_BATCH_SCRIPT, 'w') as script:
        script.write("#!/bin/tcsh\n")
        script.write("cd /cs/labs/dina/seanco/xl_parser/alphafold/colabfold/multichain/\n")
        # script.write("module load cuda\n")
        script.write("source /cs/labs/dina/seanco/xl_parser/xl_db_parser_venv/bin/activate.csh\n")
        script.write("module load cuda/11.0\n")
        script.write("module load cudnn/8.0.2\n")
        # script.write("module load cudnn\n")
        # script.write("module load tensorflow/2.5.0\n")
        script.write(COLABFOLD_BATCH_RUN + input_file_or_dir + " " + output_dir + " --data=./weights/" + " --num-recycle " + str(num_recycle) + "\n")
        script.close()
    subprocess.run("sbatch " + sbatch_args + COLABFOLD_BATCH_SCRIPT, shell=True)


def run_alphafold():
    subprocess.run("sbatch " + INLINE_INTRO + ALPHAFOLD_SCRIPT, shell=True)


def parse_input(parser):
    parser.add_argument(
        "input",
        help="Can be one of the following: "
             "Directory with fasta/a3m files, a csv/tsv file, a fasta file or an a3m file")

    parser.add_argument("results",
                        help="Directory to write the results to")

    parser.add_argument(
        "--num_recycle",
        help="Number of prediction cycles."
             "Increasing recycles can improve the quality but slows down the prediction.",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--gpu_memory",
        help="Amount of GPU to use",
        type=str,
        default='16',
    )

    parser.add_argument(
        "--time",
        help="Time to run on cluster. for 1 hour: 1:0:0. for 1 day: 1-0",
        type=str,
        default='1:0:0',
    )

    parser.add_argument("--num-models", type=int, default=5, choices=[1, 2, 3, 4, 5])


def main():
    # run_colabfold()
    parser = ArgumentParser()
    parse_input(parser)
    args = parser.parse_args()
    if len(sys.argv) < 3:
        print("usage: python3 cluster_script.py <input file / dir> <output dir> <optional - num cycles> <optional - gpu memory> <optional - time>")
        exit(1)
    run_colabfold_batch(args.input, args.results, num_recycle=args.num_recycle, gpu_memory=args.gpu_memory,
                        _ask_time=args.time)


if __name__ == "__main__":
    run_cross_link()
