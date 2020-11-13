"""
GRIT 1.0
Genome Recoding Informatics Toolbox

S. Casper
scasper@college.harvard.edu
PI: George Church, mentor: Eriona Hysolli
"""

# see key params at the top of GRIT_util to set
import GRIT_utils
import argparse

# genetic code, bases, stops, and chromosomes
gen_code = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N',
            'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R', 'CTA': 'L', 'CTC': 'L',
            'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q',
            'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G',
            'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F',
            'TTA': 'L', 'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*',
            'TGG': 'W', '': 'X'}
stop_codons = ['TAA', 'TAG', 'TGA']
bases = ['A', 'C', 'G', 'T']
chromosomes = list(range(1, 23)) + ['X', 'Y']


########################################################################################################################
# DEMO
########################################################################################################################

def demo(data_dir='data/'):

    out = ''

    # set parameters and put them into options dict
    cbe_window = (1, 13)  # what editing window to consider an NG CBE to have
    abe_window = (1, 13)  # what editing window to consider an NG ABE to have
    primer_size_range = [200, 220]  # this will make PCR products of 200-220 bases

    # Here's now to instantiate a chromosome object. This will also automatically instantiate gene objects for every
    # gene in the chromosome. You probably won't want to instantiate gene objects yourself. Let the chromosome class
    # do it for you.
    c2 = GRIT_utils.Chromosome(2, data_dir, cbe_window, abe_window, primer_size_range)

    # Here's how to get the gene name
    c_name = c2.name
    out += 'Chromosome name (str): ' + c_name + '\n\n'

    # Here's how to get a list of all + and - strand TA->G<- sites to edit in the chromosome:
    sites = c2.all_sites
    out += 'Chromosome ' + c_name + ' all positive strand sites (list of ints): ' + str(sites[0]) + '\n\n'
    out += 'Chromosome ' + c_name + ' all negative strand sites (list of ints): ' + str(sites[1]) + '\n\n'

    # Here's how to get a list of all + and - strand TA->G<- sites to edit in the chromosome with NG C base editors:
    be_sites = c2.be_sites
    out += 'Chromosome ' + c_name + ' positive strand NG be sites (dict{int index: list [ABE guides, CBE guides]}): ' \
           + str(be_sites[0]) + '\n\n'
    out += 'Chromosome ' + c_name + ' negative strand NG be sites (dict{int index: list [ABE guides, CBE guides]}): ' \
           + str(be_sites[1]) + '\n\n'

    # If you want to find sites in the chromosome that, when edited, could cause trouble because they are also part of
    # another coding sequence but aren't a stop, you can try this. The result is a dictionary with site keys and list
    # values that give a pair of genes. The first is the one with the stop codon to edit. The second is the one with
    # the conflict.
    multifunctional = c2.multifunctional_sites
    out += 'Chromosome ' + c_name + ' multifunctional sites (dict{index:[str recoding gene, str conflict gene]}): ' + \
           str(multifunctional) + '\n\n'

    # The chromosome obj has a dict attribute called genes with key names and gene values. Here's how to get a gene obj.
    bcl = c2.genes['BCL2L11']

    # The gene object contains useful attributes and functions for recoding.

    # Here's how to get the gene name. BCL2L11 is a great example.
    g_name = bcl.name
    out += 'Gene name (str): ' + g_name + '\n\n'

    # Here's how to find which strand the gene is in. True means positive, and False means reverse
    g_strand = bcl.strand
    out += g_name +  ' strand (bool): ' + str(g_strand) + '\n\n'

    # Here's how to find gene essentiality data on the gene from the OGEE database
    # http://ogee.medgenius.info/browse/Homo%20sapiens)
    g_essentiality = bcl.essentiality
    out += g_name + ' essentiality (OGEE database) (list: [essential_results, nonessential_results]): ' + \
           str(g_essentiality) + '\n\n'

    # Here's how to get the exon start and stop indices for a gene's isoforms. Let's do the first isoform here.
    # Isoforms are sorted by length.
    g_iso1 = bcl.isos[0]
    out += g_name + ' isoform 1 (list of lists of int index pairs): ' + str(g_iso1) + '\n\n'

    # Here's how to get a list of all TA->G<- sites to edit in the gene:
    g_sites = bcl.all_sites
    out += g_name + ' all sites (list of ints): ' + str(g_sites) + '\n\n'

    # Here's how to get a list of all TA->G<- sites to edit in the gene with NG C base editors:
    g_be_sites = bcl.be_sites
    out += g_name + ' NG be sites (dict{int index: str target}): ' + str(g_be_sites) + '\n\n'

    # Here's how to get the natural and recoded coding dna sequence for the first isoform of a gene
    g_cds = bcl.cds
    g_recoded_cds = bcl.recoded_cds
    out += g_name + ' wt cds of isoform 1 (str): ' + str(g_cds) + '\n\n'
    out += g_name + ' recoded cds of isoform 1 (str): ' + str(g_recoded_cds) + '\n\n'

    # Here's how to get info for primers to verify all of the sites in this gene to edit.
    g_primers = bcl.get_primer_info()
    g_primers_results = {}
    for p_name, result in g_primers.items():
        if 'PAIR_1' in result:
            g_primers_results['left'] = result['PAIR_1']['LEFT']
            g_primers_results['right'] = result['PAIR_1']['RIGHT']
            g_primers_results['product_len'] = result['PAIR_1']['PRODUCT_LEN']
            g_primers_results['reference'] = result['PAIR_1']['REFERENCE']
    out += g_name + ' primer info (dict of various key and value types): ' + str(g_primers_results)

    return out

########################################################################################################################
# RESULT REPLICATION
########################################################################################################################


def count_total_sites(data_dir='data/'):  # replicate results for finding all sites in the genome

    # set parameters and put them into options dict
    cbe_window = (1, 13)  # what editing window to consider an NG CBE to have
    abe_window = (1, 13)  # what editing window to consider an NG ABE to have
    primer_size_range = [200, 220]  # this will make PCR products of 200-220 bases

    # init counters
    site_count = 0

    for c_name in chromosomes:
        print('processing chromosome', c_name, '...')
        c = GRIT_utils.Chromosome(c_name, data_dir, cbe_window, abe_window, primer_size_range)
        # update counters
        site_count += len(c.all_sites[0]) + len(c.all_sites[1])

    # display info
    print('total TAG sites:', site_count)

    return site_count


def count_editing_sites(data_dir='data/'):

    # returns the total number of essential sites and the number that can be edited with C/A base editors

    # set parameters and put them into options dict
    cbe_window = (1, 13)  # what editing window to consider an NG CBE to have
    abe_window = (1, 13)  # what editing window to consider an NG ABE to have
    primer_size_range = [200, 220]  # this will make PCR products of 200-220 bases

    # for counting sites
    total = 0
    to_edit = 0
    to_skip = 0

    for c_name in chromosomes:
        print('processing chromosome', c_name, '...')
        c = GRIT_utils.Chromosome(c_name, data_dir, cbe_window, abe_window, primer_size_range)

        for name, g in c.genes.items():
            # don't worry about genes that are consensus or near-consensus nonessential
            if g.essentiality != 'unavailable':
                ess = g.essentiality[0]
                noness = g.essentiality[1]
                if ess == 0 or (ess == 1 and noness >= 6):
                    continue
            if len(g.all_sites) > 0:
                total += len(g.all_sites)
                for v in g.be_sites.values():
                    if v == [-1, -1]:
                        to_skip += 1
                    else:
                        to_edit += 1

    return total, to_edit, to_skip


def find_genes_to_recode(data_dir='data/'):

    # returns a list of all essential names that are essential

    # set parameters and put them into options dict
    cbe_window = (1, 13)  # what editing window to consider an NG CBE to have
    abe_window = (1, 13)  # what editing window to consider an NG ABE to have
    primer_size_range = [200, 220]  # this will make PCR products of 200-220 bases

    genes_to_recode = []  # list to store names as strings

    for c_name in chromosomes:
        print('processing chromosome', c_name, '...')
        c = GRIT_utils.Chromosome(c_name, data_dir, cbe_window, abe_window, primer_size_range)

        for name, g in c.genes.items():
            if g.essentiality != 'unavailable':
                ess = g.essentiality[0]
                noness = g.essentiality[1]
                if ess == 0 or (ess == 1 and noness >= 6):  # if consensus noness or one ess with at leasst 6 noness
                    continue
            if len(g.be_sites) > 0:
                genes_to_recode.append(g.name)

    return genes_to_recode


def get_all_site_data(data_dir='data/'):

    # returns a list off 24 lists which each give the indices of all TAG sites in the chromosomes

    all_sites = []

    # set parameters and put them into options dict
    cbe_window = (1, 13)  # what editing window to consider an NG CBE to have
    abe_window = (1, 13)  # what editing window to consider an NG ABE to have
    primer_size_range = [200, 220]  # this will make PCR products of 200-220 bases

    for c_name in chromosomes:
        print('processing chromosome', c_name, '...')
        c = GRIT_utils.Chromosome(c_name, data_dir, cbe_window, abe_window, primer_size_range)

        chromosome_sites = c.all_sites[0] + c.all_sites[1]
        chromosome_sites.sort()
        chromosome_sites.append(len(c.seq))
        all_sites.append(chromosome_sites)

    return all_sites

########################################################################################################################
# MAIN
########################################################################################################################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genome Recoding Informatics Toolbox (GRIT)')
    parser.add_argument('-f', '--function', action='store', dest='function',
                        help="function to run in [\'w\', \'count_total_sites\', \'count_editing_sites\', "
                             "\'find_genes_to_recode\'], or custom")
    parser.add_argument('-w', '--write_file', action='store', dest='write_file',
                        help=".txt file to write to")
    args = parser.parse_args()

    if args.write_file[-4:] != '.txt':
        args.write_file += '.txt'

    if args.function == 'demo':
        output = 'GRIT DEMO\n\n' + demo() + '\n'

    elif args.function == 'count_total_sites':
        output = 'GRIT COUNT_TOTAL_SITES \n\n' + str(count_total_sites()) + '\n'

    elif args.function == 'count_editing_sites':
        counts_out = count_editing_sites()
        output = 'GRIT COUNT_EDITING_SITES \n\n' + 'Total essential: ' + str(counts_out[0]) + \
                 ', Total editable: ' + str(counts_out[1]) + ', Total not editable: ' + str(counts_out[1]) + '\n'

    elif args.function == 'find_genes_to_recode':
        output = 'GRIT FIND_GENES_TO_RECODE \n\n' + str(find_genes_to_recode()) + '\n'

    else:
        raise ValueError('Specified function must be in [\'demo\', \'count_edit_sites\', '
                         '\'find_genes_to_recode\'] or custom added function.')

    with open(args.write_file, 'w') as f:
        f.write(output)
