from sklearn import preprocessing
from collections import defaultdict
import re
import time
import scipy.stats as stats
import copy
import networkx as nx
import matplotlib.pyplot as plt

def get_gene_expression_data(file_name='GSE10072_series_matrix.txt', invalid_probe=[], apply_z_transform=True):
    '''
    :param fileName: { string }
    :param invalid_probe: { list }
    :return: { dict - {probe: [ gene expression, ...]}
    '''
    dict_temp = {}
    with open(file_name) as f:
        for line in f:
            if re.match(r'(.*)_at', line):
                line_split = line.split()
                probe = line_split[0].strip('"')
                if probe not in invalid_probe:
                    if apply_z_transform:
                        dict_temp[probe] = _z_transform(map(float, line_split[1:]))
                    else:
                        dict_temp[probe] = map(float, line_split[1:])
    return dict_temp


def get_invalid_probe(probe_gene_matching_file_name='GPL96.txt'):
    # Remove probes represented in zero or more than one genes
    '''
    :param probe_gene_matching_file_name: { string }
    :return: { list, dict }
    '''
    invalid_probe = []
    valid_probe = {}
    with open(probe_gene_matching_file_name) as f:
        for line in f:
            if re.match(r'(.*)_at', line):
                line_split = line.split()
                if len(line_split) == 1 or '///' in line_split[1]:
                    invalid_probe.append(line_split[0])
                else:
                    valid_probe[line_split[0]] = line_split[1]
    return invalid_probe, valid_probe


def _z_transform(gene_expression):
    '''
    :param gene_expression: { list }
    :return: z_transform gene expression : { list }
    '''
    return preprocessing.scale(gene_expression)


def get_gene_set(file_name='pathway-api.symbols.gmt.txt'):
    '''
    :param file_name: { string }
    :return: Gene set list { [['GALM', 'G6PC', ...], [ ... ], [ ... ], ...] }
    '''
    gene_set_list = []
    with open(file_name) as f:
        for line in f:
            gene_set_list.append(line.split()[2:])

    return gene_set_list


def map_probe_into_gene(probe_expression, probe_gene_matching):
    gene_probe_expression = defaultdict(dict)

    for probe in probe_expression:
        if probe in probe_gene_matching:
            gene_probe_expression[str(probe_gene_matching[probe])][probe] = probe_expression[probe]

    return gene_probe_expression


def gene_network_mapping(gene_set, ppi_file_name='PPI_04052016.sif.txt'):
    G = nx.Graph()
    G.add_nodes_from(gene_set)
    pos = nx.spring_layout(G) # positions for all nodes

    # print G.nodes()

    gene_pair_ppi = []
    with open(ppi_file_name) as f:
        for line in f:
            gene_pair = line.split()
            gene_pair_ppi.append((gene_pair[0], gene_pair[1]))

    # each gene set
    gene_set_pair = []
    
    labels = {}
    for i in range(len(gene_set) - 1):
        gene_set_pair.append((gene_set[i], gene_set[i+1]))
        gene_set_pair.append((gene_set[i+1], gene_set[i]))
        labels[i] = gene_set[i]
    
    labels[len(gene_set) - 1] = gene_set[len(gene_set)-1]

    for i in gene_set_pair:
        if i in gene_pair_ppi:
            G.add_edge(i[0], i[1])

    nx.draw(G, edge_color='#909090', node_size=50)
    plt.show()


def one_way_anova(tumor_list, normal_list):
    '''
    :param tumor_list: { list }
    :param normal_list: { list }
    :return: f-value: { float }, p-value: { float }
    '''
    f_value, p_value = stats.f_oneway(tumor_list, normal_list)
    return f_value, p_value


# Modify later
def get_sample_disease_pairs(file_name='GSE10072_label.txt'):
    '''
    :param file_name: { string }
    :return: list of pair of sample and disease: { list - [(sample, disease), (), ...] }
    '''
    sample = []
    disease = []
    sample_with_disease = []
    with open(file_name) as f:
        for line in f:
            if re.match(r'Sample', line):
                sample = line.split()[1:]
            elif re.match(r'Disease', line):
                disease = line.split()[1:]

    for i in range(len(sample)):
        sample_with_disease.append((sample[i], disease[i]))

    return sample_with_disease


# Modify later
def get_tumor_normal_index(sample_disease_pairs_list):
    '''
    :param sample_disease_pairs_list: { list - [(sample, disease), (), ...]}
    :return: tumor's index list: { list }, normal's index list: { list }
    '''
    tumor_index = []
    normal_index = []
    for i in range(len(sample_disease_pairs_list)):
        if sample_disease_pairs_list[i][1] == 'Tumor':
            tumor_index.append(i)
        else:
            normal_index.append(i)

    return tumor_index, normal_index


def get_tumor_normal_list(probe_expression):
    tumor_index = [0, 2, 4, 5, 6, 8, 11, 12, 14, 16, 17, 20, 22, 23, 25,
                   27, 29, 31, 32, 34, 36, 38, 39, 41, 43, 45, 47, 49, 50,
                   53, 55, 57, 59, 61, 62, 63, 65, 67, 69, 71, 72, 73, 75,
                   76, 79, 80, 82, 84, 89, 91, 93, 95, 96, 97, 99, 101, 103, 104]

    normal_index = [1, 3, 7, 9, 10, 13, 15, 18, 19, 21, 24, 26, 28, 30, 33,
                    35, 37, 40, 42, 44, 46, 48, 51, 52, 54, 56, 58, 60, 64,
                    66, 68, 70, 74, 77, 78, 81, 83, 85, 86, 87, 88, 90, 92,
                    94, 98, 100, 102, 105, 106]

    tumor_list = []
    normal_list = []
    for i in tumor_index:
        tumor_list.append(probe_expression[i])

    for i in normal_index:
        normal_list.append(probe_expression[i])

    return tumor_list, normal_list


def get_gene_network_mapping(file_name='ppi.txt'):
    interactive_gene = defaultdict(list)

    with open(file_name) as f:
        for line in f:
            genes_list = line.split()
            interactive_gene[genes_list[0]].append(genes_list[1])

    return interactive_gene


def GSNFS():
    #---------------------- All needed files ----------------------#
    #Gene probe expression data
    raw_gene_probe_expression_data = 'GSE10072_series_matrix.txt'

    #Sample disease - significant gene
    sample_disease = 'GSE10072_label.txt'

    #Gene-probe matching - to check invalid probe
    gene_probe_matching = 'GPL96.txt'

    #Gene set
    gene_set = 'pathway-api.symbols.gmt.txt'

    #PPI - is used for Gene-network mapping
    ppi = 'PPI_04052016.sif.txt'


    #---------------------- Computational Time ----------------------#
    invalid_probe, valid_probe = get_invalid_probe(gene_probe_matching)

    probe_expression_dict = get_gene_expression_data(raw_gene_probe_expression_data, invalid_probe, True)
    #---------------------- Remove probes having p-value > 0.05 ----------------------#
    copy_probe_expression = probe_expression_dict.copy()
    for p in copy_probe_expression:
        tumor, normal = get_tumor_normal_list(copy_probe_expression[p])
        f_v, p_v = one_way_anova(tumor, normal)

        if p_v > 0.05:
            del probe_expression_dict[p]
    
    gene_set_list = []
    with open(gene_set) as f:
        for line in f:
            gene_set_list.append(line.split()[2:])
    
    gene_network_mapping(gene_set_list[1])

    #gene_probes = map_probe_into_gene(probe_expression_dict, valid_probe)
    #print gene_probes

    #gene_set_list = get_gene_set(gene_set)
    # Have to loop through for each gene set
    #for gsl in  gene_set_list:
        #print gsl
    #print gene_set_list[0]
    print "-----------------------------------------------------"
    sample_disease_pairs = get_sample_disease_pairs(sample_disease)
    tumor_index, normal_index = get_tumor_normal_index(sample_disease_pairs)
    print tumor_index, normal_index
    f, p = one_way_anova(tumor_index, normal_index)
    print f, p
GSNFS()