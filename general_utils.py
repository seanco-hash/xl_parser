import pickle
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gc
import json
import matplotlib.pyplot as plt
import seaborn as sb
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels
import sklearn.metrics as metrics
from matplotlib.ticker import StrMethodFormatter

PROJ_DIR = '/cs/labs/dina/seanco/xl_parser/'
OBJ_DIR = '/cs/labs/dina/seanco/xl_parser/obj/'
# PROJ_DIR = 'C:/Users/seanco/Desktop/University/Master/Research/XL_DB_parser/'

SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 24
VERY_SMALL_SIZE = 12

XL_PEP_A = 0
XL_POS_IN_PEP_A = 1
XL_PEP_ID_A = 2
XL_UNIPORT_A = 3
XL_ACCESSION_A = 4
XL_RES_NUM_A = 5
XL_GENE_NAME_A = 6
XL_PEP_START_POS_A = 7
XL_PDB_A = 8
XL_PDB_TYPE_A = 9
XL_SITE_A = 10
XL_ATOM_NUM_A = 11
XL_PEP_B = 12
XL_POS_IN_PEP_B = 13
XL_PEP_ID_B = 14
XL_UNIPORT_B = 15
XL_ACCESSION_B = 16
XL_RES_NUM_B = 17
XL_GENE_NAME_B = 18
XL_PEP_START_POS_B = 19
XL_PDB_B = 20
XL_PDB_TYPE_B = 21
XL_SITE_B = 22
XL_ATOM_NUM_B = 23
XL_DISTANCE = 24
XL_NETWORK_DISTANCE = 25
XL_NUM_IDS = 26
XL_CONFIDENCE = 27
XL_DATASETS = 28
XL_DATASETS_COUNT = 29
XL_STRUCTURES = -1

XL_AVAILABLE_IN_PDB = 'Structure'

REGULAR_PDBS = 0
AF_PDB = 1
KNOWN_STRUCTURE_PDB = 2

DSSO_LINKER = 1
BDP_NHP_LINKER = 2

FASTA_DICT_NAME = 'fasta_files_dict_by_uniport'

linker_type_dict = {'DSSO' : DSSO_LINKER, 'BDP-NHP' : BDP_NHP_LINKER}


def pdb_from_sample(sample, af_accessions):
    if sample[XL_PDB_TYPE_A] == XL_AVAILABLE_IN_PDB:
        return [sample[XL_PDB_A]], KNOWN_STRUCTURE_PDB
    elif af_accessions is not None:
        cur_accession = sample[XL_ACCESSION_A].split('-')[0]
        if cur_accession in af_accessions:
            return [cur_accession], AF_PDB
    if sample[XL_STRUCTURES] == '':
        return [], REGULAR_PDBS
    return sample[-1].split(','), REGULAR_PDBS


def save_np_obj(obj, name, dir_=OBJ_DIR):
    with open(dir_ + name + '.npy', 'wb') as f:
        np.save(f, obj)


def load_np_obj(name, dir_=OBJ_DIR):
    with open(dir_ + name + '.npy', 'rb') as f:
        return np.load(f)


def save_obj(obj, name, dir_=OBJ_DIR):
    with open(dir_ + name + '.pkl', 'wb') as f:
        gc.disable()
        pickle.dump(obj, f, protocol=4)
        gc.enable()


def load_obj(name, dir_=OBJ_DIR):
    with open(dir_ + name + '.pkl', 'rb') as f:
        gc.disable()
        obj = pickle.load(f)
        gc.enable()
    return obj


def save_obj_json(obj, name, dir_=OBJ_DIR):
    with open(dir_ + name + '.json', 'w') as fp:
        json.dump(obj, fp)


def load_obj_json(name, dir_=OBJ_DIR):
    with open(dir_ + name + '.json', 'r') as fp:
        obj = json.load(fp)
    return obj


def seqs_dict_to_fasta(seqs, fasta_name):
    records = list()
    for uni, seq in sorted(seqs.items(), key=lambda x: len(x[1])):
        print(f"{uni}: {len(seq)}")
        record = SeqRecord(Seq(seq), id=uni, name=uni, description="")
        records.append(record)
    with open(fasta_name, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")


def initialize_plt_params():
    plt.rcParams.update({'axes.facecolor':'white'})
    # plt.gcf().subplots_adjust(bottom=0.15)
    # plt.gcf().subplots_adjust(left=0.15)
    # plt.gcf().subplots_adjust(top=0.8)
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=VERY_SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
    plt.grid(color='gray', linewidth=0.1)
    # plt.figure(dpi=300)


def find_complexes_rank(rmsd, score, title):
    rmsd_prob_score = np.array([x for x, _ in sorted(zip(rmsd, score), key=lambda pair: pair[1])]).flatten()
    # idx = np.argsort(rmsd_prob_score)[:20]
    idx = np.argmax(rmsd_prob_score < 10)
    print(f"{title} First good solution: {idx}")


def plot_scatter(x_data, y_data, title, x_title, y_title, dot_size=0.7, color_arr=None, mark_low_th=10, save_path=None):
    initialize_plt_params()
    cmap = plt.get_cmap('Blues', 1)
    cmap.set_under('red')
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    if mark_low_th >= 0:
        plt.scatter(x_data, y_data, s=dot_size, c=color_arr, cmap=cmap, vmin=mark_low_th)
    else:
        plt.scatter(x_data, y_data, s=dot_size, c=color_arr, cmap=cmap)
    # plt.title(title, fontsize=MEDIUM_SIZE)
    # plt.xlim(0, 120)
    # plt.xlabel(x_title)
    # plt.ylabel(y_title)
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def interpolate_scores(scores, prob_weight, reg_weight, filter_th=0, soap_score=None, soap_weight=0.2):
    score_reg = scores[0]
    score_prob = scores[1]
    score_reg = (score_reg - np.mean(score_reg)) / np.std(score_reg)
    score_prob = (score_prob - np.mean(score_prob)) / np.std(score_prob)
    filter_idx = np.argwhere(score_reg < filter_th)
    score_reg = score_reg[filter_idx]
    score_prob = score_prob[filter_idx]
    if soap_score is not None:
        soap_score = np.array(soap_score)[filter_idx]
    score = (reg_weight * score_reg) - (prob_weight * score_prob)
    if soap_score is not None and soap_weight > 0:
        score += soap_score * soap_weight
    return score, filter_idx, soap_score


def get_color_arr(rmsd, interface_rmsd, filter_idx=None):
    color_arr = np.ones(rmsd.shape)
    color_arr *= 20
    color_arr[rmsd <= 10] = 0
    if len(interface_rmsd) > 0:
        if filter_idx is not None:
            interface_rmsd = interface_rmsd[filter_idx]
        color_arr[interface_rmsd <= 4] = 0
    return color_arr


def plot_funnel(r=4, s=3, v=2, t="TRIC",
                file_path='/cs/labs/dina/seanco/Tric/tric_align/cmake-build-release/stats_prob_tric.txt',
                negate_score=True, color='c', mark_low_th=10, save=False, prob_score_weight=None, normalize=True,
                filter_th=0, soap=True, soap_weight=0.2, reg_score_weight=0.8, i_rmsd=False, plot_soap=False, dir_path=None):
    if prob_score_weight is None:
        file_path = [file_path]
    filter_idx = None
    violations = []
    scores = []
    rmsds = []
    soap_score = None
    if soap:
        soap_score = []
    interface_rmsd = []
    for file in file_path:
        f = open(file, 'r')
        violation = []
        score = []
        rmsd = []
        for line in f:
            sp = line.split(" ")
            violation.append(int(sp[v]))
            rmsd.append(float(sp[r]))
            if i_rmsd:
                interface_rmsd.append(float(sp[-2]))
            if soap:
                soap_score.append(float(sp[-1]))
            if negate_score:
                score.append(-float(sp[s]))
            else:
                score.append(float(sp[s]))
        violations.append(np.array(violation))
        scores.append(np.array(score))
        rmsds.append(np.array(rmsd))
    s_v, v_r, s_r = None, None, None
    if prob_score_weight is None:
        if save:
            s_r = dir_path + t + "_svdist.png"
        violation, rmsd, score = violations[0], rmsds[0], scores[0]
        if len(rmsd) == 0 or min(rmsd) > 15:
            print("no good models for " + t)
            return
        if normalize:
            score = (score - np.mean(score)) / np.std(score)
            if filter_th is not None:
                filter_idx = np.argwhere(score < filter_th)
                score, violation, rmsd = score[filter_idx], violation[filter_idx], rmsd[filter_idx]
                if soap:
                    soap_score = np.array(soap_score)[filter_idx]
    else:
        if save:
            if soap_weight > 0 and prob_score_weight is not None and prob_score_weight > 0:
                  s_r = dir_path + t + "_all.png"
            else:
                  s_r = dir_path + t + "_ours.png"
        violation, rmsd = violations[1], rmsds[1]
        if len(rmsd) == 0 or min(rmsd) > 15:
            print("no good models for " + t)
            return
        score, filter_idx, soap_score = interpolate_scores(scores, prob_score_weight, reg_score_weight,
                                                     filter_th, soap_score, soap_weight)
        rmsd, violation = rmsd[filter_idx], violation[filter_idx]
    if len(score) > 0:
        color_arr = get_color_arr(rmsd, np.array(interface_rmsd), filter_idx)
        plot_scatter(rmsd, score, t, "rmsd (Å)", "score", 6, color_arr, mark_low_th, s_r)
        find_complexes_rank(color_arr, score, t)
        if soap and plot_soap:
            s_r = dir_path + t + "_soap.png"
            plot_scatter(rmsd, soap_score, t + "_soap", "rmsd (Å)", "soap score", 6, color_arr, mark_low_th, s_r)
            find_complexes_rank(color_arr, soap_score, "Soap ")


def plot_single_prediction_distribution(probas, targets, n=10):
    for i in range(n):
        y = probas[i]
        plot_line(range(len(y)), y, f"distribution, y={targets[i]}", "class",
                  "probability", False)


def plot_line(x, y, title, x_labe, y_label, plot_values=True, bar=False):
    # initialize_plt_params()
    if bar:
        sb.barplot(data=x, x='labels', y='lens')
        # plt.bar(x, y)
    else:
        plt.xticks(x)
        plt.plot(x, y, marker='o')
    plt.title(title, fontsize=BIGGER_SIZE)
    plt.xlabel(x_labe)
    plt.ylabel(y_label)
    plt.legend(["H: Helix", "B: Beta", "L: Loop"])
    if plot_values:
        for index in range(len(x)):
            plt.text(x[index], y[index], y[index], size=12)
    plt.show()


def read_fasta_files(fasta_path):
    uni_seq_dict = dict()
    for record in SeqIO.parse(fasta_path, "fasta"):
        uni_seq_dict[record.id.split('.')[0]] = str(record.seq)
    return uni_seq_dict


def plot_histogram(data, title_="", x_label='Value', y_label='Frequency', xticks_values=None, xticks_names=None, bins=None,
                   data_labels=['x'], normalize=False, colors=['cadetblue']):
    initialize_plt_params()
    # bins = list(range(50)) + [52, 55, 60, 65, 70, 80, 90, 100, 125, 150, 200]
    if bins is None:
        bins = list(range(46)) #+ [105, 110, 125, 135, 150, 175, 200]
    ns = []
    for i in range(len(data)):
        # n, bins, patches = plt.hist(x=data[i], bins=bins, color=colors[i], alpha=0.7, rwidth=0.85,
        #                             label=data_labels[i], density=normalize)
        axes2 = sb.kdeplot(data=data[i], shade=True, color=colors[i])
        axes2.set(yticklabels=[])
        # ymin, ymax = axes2.get_lim()
        # axes2.vlines(x=[45], ymin=ymin, ymax=ymax, colors=['tab:red'], ls='--', lw=2, alpha=0.5)
        # ns.append(n.max())
    # maxfreq = max(ns)
    # plt.bar_label(patches)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    # plt.title(title_, fontsize=BIGGER_SIZE)
    # plt.suptitle(title_, fontsize=BIGGER_SIZE)
    # plt.text(23, 45, f"$\mu={round(np.std(data), 1)}, u={round(np.mean(data), 1)}$")
    # Set a clean upper y-axis limit.
    # plt.ylim(ymax=maxfreq + (0.1 * maxfreq))
    # plt.xticks([0, 15, 30, 50, 100, 150, 200])
    if xticks_values is not None:
        plt.xticks(xticks_values, xticks_names)
    if len(data) > 1:
        plt.legend(loc='upper right', labels=data_labels)
    plt.axvline(x=45, color='black', ls=':', lw=2)
    # plt.xlim(0, bins[-1])
    plt.show()


def plot_confusion_matrix(y_true, y_pred, classes, normalize=False, title=None, cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    # classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True Distance (Å)',
           xlabel='Predicted Distance (Å)')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    plt.xlim(-0.5, len(np.unique(y_true))-0.5)
    plt.ylim(len(np.unique(y_true))-0.5, -0.5)
    print("Finish plot cm")
    return ax


def plot_roc(labels, predictions, title, n_classes, class_names):
    initialize_plt_params()
    one_hot_labels = np.zeros((labels.size, labels.max() + 1))
    one_hot_labels[np.arange(labels.size), labels] = 1
    try:
        auc = metrics.roc_auc_score(one_hot_labels, predictions, multi_class='ovo', average=None)
        for label in range(n_classes):
            preds = predictions[:, label]
            fpr, tpr, _ = metrics.roc_curve(labels, preds, pos_label=label)
            a = auc[label]
            plt.plot(fpr, tpr, label=f"{class_names[label]}, auc={a:.3f}")
        print(auc)
        plt.legend(loc=4)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title, fontsize=BIGGER_SIZE)
        # plt.savefig(os.path.join(self.fig_save_path, self.title))
        # plt.show()
        return plt, auc
    except Exception as e:
        print(e)
        return None, None


def chatgpt_ablation_plot():
    # define the data for the plot
    components = ['A', 'B', 'C', 'D']
    baseline_performance = 0.8
    ablation_results = [0.7, 0.75, 0.6, 0.78]

    # create a bar plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(components, ablation_results, color='gray')

    # add a horizontal line for the baseline performance
    ax.axhline(y=baseline_performance, color='black', linestyle='--')

    # set the axis labels and title
    ax.set_xlabel('Components', fontsize=14)
    ax.set_ylabel('Performance', fontsize=14)
    ax.set_title('Ablation Study Results', fontsize=16)

    # add text labels for each bar
    for i, result in enumerate(ablation_results):
        ax.text(i, result + 0.01, f'{result:.2f}', ha='center', fontsize=12)

    # add a legend for the baseline and modified bars
    ax.legend(['Baseline', 'Modified'], fontsize=12)

    # adjust the layout of the plot
    fig.tight_layout()

    # show the plot
    plt.show()

def af_plot_ablation():
    # data for the plot
    names = ['Baseline', 'SA', 'MSA', 'Recycle', 'Recycle+MSA']
    accuracy = [57.7, 58.2, 59.6, 61.1, 62.6]
    runtime = [1.0, 1.7, 23.6, 7.9, 166.1]

    # create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

    # create the accuracy plot on the first subplot
    ax1.bar(names, accuracy)
    ax1.set_ylim([0, 70])
    ax1.set_ylabel('CA', fontsize=12)
    ax1.tick_params(axis='y', labelsize=10)
    ax1.set_title('Accuracy (%)', fontsize=14)

    # create the runtime plot on the second subplot
    ax2.bar(names, runtime)
    ax2.set_yscale('log')
    ax2.set_ylim([0.1, 1000])
    ax2.set_ylabel('Time (s)', fontsize=12)
    ax2.tick_params(axis='y', labelsize=10)
    ax2.set_title('Runtime (s)', fontsize=14)

    # adjust the layout of the subplots and save the figure
    fig.tight_layout()
    plt.show()


def plot_ablation():
    # bar_width = 0.35
    initialize_plt_params()
    x_single = [1, 2, 3, 4, 5, 6, 7, 8]
    y_single = [0.6195, 0.6115, 0.6119, 0.6290, 0.6797, 0.7267, 0.8078, 0.6862]

    x_ticks_single = ['anchors', 'radius', 'charge', 'id', 'asa', 'mol2', 'res_type', 'ss']
    plt.bar(x_single, y_single, color='cornflowerblue')
    # plt.bar([p + bar_width for p in x_single], auc_inter_single, color='navy')
    plt.xticks(x_single, x_ticks_single,rotation='vertical')
    plt.tight_layout()
    plt.ylabel("Average AUC")
    plt.axhline(y=0.857519, color='black',  linestyle='--')
    plt.ylim([0, 1])
    plt.show()

    auc_inter_single = [0.6213, 0.6264, 0.6269, 0.6334, 0.5690, 0.631, 0.7177, 0.6088]
    plt.axhline(y=0.6958, color='black', linestyle='--')
    plt.bar(x_single, auc_inter_single, color='cornflowerblue')
    plt.xticks(x_single, x_ticks_single, rotation='vertical')
    plt.ylim([0, 1])
    plt.tight_layout()
    plt.show()

    x_all_but = [7, 6, 4, 5, 3, 2, 1, 8]
    y_all_but = [0.8012, 0.8454, 0.8312, 0.8627, 0.84531, 0.8397, 0.8208, 0.8799]
    x_all_but_ticks = ['res-type', 'mol2', 'id', 'asa', 'charge', 'radius', 'anchors', 'ss']
    plt.bar(x_all_but, y_all_but, color='cornflowerblue')
    plt.axhline(y=0.857519, color='black', linestyle='--')
    plt.xticks(x_all_but, x_all_but_ticks, rotation='vertical')
    plt.tight_layout()
    plt.ylim([0, 1])
    plt.ylabel("Average AUC")
    plt.show()

    auc_inter_all_but = [0.6334, 0.6691, 0.6414, 0.7021, 0.6770, 0.6677, 0.6764, 0.5686]
    plt.bar(x_all_but, auc_inter_all_but, color='cornflowerblue')
    plt.axhline(y=0.6958, color='black', linestyle='--')
    plt.xticks(x_all_but, x_all_but_ticks, rotation='vertical')
    plt.tight_layout()
    plt.ylim([0, 1])
    plt.show()
    # plt.bar([p + bar_width for p in x_all_but], auc_inter_all_but, color='darkred')

    y_all = [0.8575195567241383]
    inter_auc_all = [0.6958]
    # inter_auc_all = [0.6483661204104554]


def box_plot(arr, labels, title='Residue SS types to distance'):
    fig7, ax7 = plt.subplots()
    ax7.set_title(title)
    ax7.boxplot(arr, notch=True, labels=labels)
    # plt.savefig("/cs/labs/dina/seanco/xl_parser/plots/ss_dist_box_plot.png")
    plt.show()
