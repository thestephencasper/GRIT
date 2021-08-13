from array import array  # this will be needed for finding the reverse complement of strings and recoding chromosomes
import primer3 as p3  # this will be used in the get_cas9_site_primers function in the Gene class
import pandas as pd  # this will be used to process gene essentiality data retrieved from a .csv

########################################################################################################################

# This isn't generally recommended, but set this to True if you want to find primers for every single gene when
# instantiating them, however, this will slow down runtime, and if you only need primers for some genes, you can call
# the get_cas9_site_primers() function in the Gene class for certain genes AFTER instantiating a chromosome.
find_all_primers = False

# genetic code, bases, stops
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

########################################################################################################################


def get_essentiality_data(ess_data_dir):

    # Get genomic gene essentiality data downloaded from the OGEE database
    # http://ogee.medgenius.info/browse/Homo%20sapiens

    ogee_data = pd.read_csv(ess_data_dir + 'Homo sapiens_consolidated.csv')  # csv from OGEE with essentiality data
    essentiality_data = {}  # dict to store gene keys and essentiality data values
    symbols = ogee_data['symbols'].tolist()  # gene names
    status = ogee_data['essentiality status'].tolist()  # essentiality data
    for i in range(len(symbols)):  # for each gene from the csv
        e_count = status[i].count('E')  # get number of 'E's
        ne_count = status[i].count('NE')  # get number of 'NE's
        essentiality_data[symbols[i]] = [e_count-ne_count, ne_count]  # fill in data into dict
    return essentiality_data


def get_name(line):
    """
    :param line: a str representing a header line from a genbank file with information on a gene CDS
    :return: str giving the name of that gene
    """
    name_start = line.find('gene=') + 5
    name_end = line.find(']', name_start)
    gene_name = line[name_start:name_end]
    return gene_name


def get_strand(line):
    """
    :param line: a str representing a header line from a genbank file with information on a gene CDS
    :return: bool giving the strand of that gene
    """
    strand_start = line.find('location=') + 9
    if line[strand_start:strand_start+4] == 'comp':  # (complement) if negative
        return False
    else:  # if positive
        return True


def get_idxs(line):
    """
    :param line: a str representing a line from a genbank file with information on a gene CDS
    :return: a list of lists where each inner list is a pair of inds giving the start and stop of a coding segment
    """
    inds = []  # the list of lists to be returned
    curr_start = ''  # string storing the numbers for the start being added
    curr_stop = ''  # string storing the numbers for the stop being added
    on_start = True  # a switch that tells if we are adding to a start or stop
    i = line.find('location=') + 9  # start i as index near where inds begin

    while True:  # this small loop positions i at the first number
        if line[i].isnumeric():
            break
        elif line[i] in '<>':
            return -1  # when there is an error in genbank file that doesn't give hard bounds on segments
        i += 1

    while True:  # this loop builds the list of start / stop pairs
        if line[i].isnumeric() and on_start:  # add digit to start
            curr_start += line[i]
        elif line[i].isnumeric() and not on_start:  # add digit to stop
            curr_stop += line[i]
        elif not line[i].isnumeric() and on_start:  # skip over the '..'
            if line[i] == '.':  # if this start index is followed by a stop
                i += 1
                on_start = False
            else:  # if the start index stands alone; reset and append
                inds.append([int(curr_start)-1, int(curr_start)])  # -1 for pythonic indexing
                curr_start = ''
                curr_stop = ''
        elif line[i] in '<>':  # when there is an error in genbank file that doesn't give hard bounds on segments
            return -1
        else:  # reset and append the start and stop to the list
            on_start = True
            inds.append([int(curr_start)-1, int(curr_stop)])  # -1 for pythonic indexing
            curr_start = ''
            curr_stop = ''
        if line[i] == ')' or line[i] == ']':  # signifies the end of the frame
            break
        i += 1  # increment the position

    return inds


def rc(sequence):
    """
    :param sequence: a string made of chars ACGTBDHVRYSW
    :return: str giving reverse complement of seq
    """
    reverse = array('B')  # array mod imported above
    for b in reversed(sequence):  # build reverse complement base by base
        if b == 'A':
            reverse.append(84)
        elif b == 'T':
            reverse.append(65)
        elif b == 'G':
            reverse.append(67)
        elif b == 'C':
            reverse.append(71)
        elif b == 'B':
            reverse.append(86)
        elif b == 'D':
            reverse.append(72)
        elif b == 'H':
            reverse.append(68)
        elif b == 'V':
            reverse.append(66)
        elif b == 'R':
            reverse.append(89)
        elif b == 'Y':
            reverse.append(82)
        elif b == 'S':
            reverse.append(87)
        elif b == 'W':
            reverse.append(83)
        else:
            reverse.append(ord(b))
    reverse = str(reverse.tobytes())[2:-1]  # convert to a string and trim off the extra characters
    return reverse


def order_around_mid(sites, mid):
    """
    :param sites: a list
    :param middle of window
    :return: a list with the same elements as the input but ordered with the middlemost elements first and last last
    """

    if not sites:
        return []

    return_list = []  # list to build and return

    diffs = [abs(mid - pos) for pos in sites]  # get the distances from each to the middles

    m = max(diffs) + 1

    for i in range(len(sites)):
        closest_idx = diffs.index(min(diffs))  # find smallest diff from the middle
        return_list.append(sites[closest_idx])  # add to return_list
        diffs[closest_idx] += m  # add the max to the index so that it's not found again to be the closest

    return return_list


class Chromosome:

    def __init__(self, chromosome_id, data_dir, cbe_window, abe_window, primer_size_range):
        """
        :param chromosome_id: str or int identifying the chromosome
        :param data_dir: str giving directory
        :param cbe_window: tuple giving edit range for c base editors
        :param abe_window: tuple giving edit range for a base editors
        :param primer_size_range: tuple giving acceptable range of primer sizes
        """
        self.name = str(chromosome_id).upper()  # chromosome name
        self.data_dir = data_dir
        self.cbe_window = cbe_window
        self.abe_window = abe_window
        self.primer_sizes = primer_size_range
        self.c_mid = sum(self.cbe_window) / 2
        self.a_mid = sum(self.abe_window) / 2
        self.seq = self.get_seq()  # original chromosome
        self.genes = self.get_genes_as_lists()  # dict with gene names as keys and strand and isos in a list as values
        self.all_sites, self.gene_all_sites = self.get_all_sites()  # get all edit sites
        self.be_sites, self.gene_be_sites = self.get_editor_sites()  # get base editable sites
        self.recoded_seq = self.get_recoded_seq()  # recoded chromosome
        self.essentiality_data = get_essentiality_data(self.data_dir)
        self.genes = self.get_genes_as_objects()  # get gene objects, and store by overwriting the old self.genes
        self.multifunctional_sites = self.get_multifunctional_sites()  # get dictionary of multifunctional sites

        # standard initialization
        p3.bindings.setP3Globals({'PRIMER_OPT_SIZE': 20,
                                  'PRIMER_PICK_INTERNAL_OLIGO': 1,
                                  'PRIMER_INTERNAL_MAX_SELF_END': 8,
                                  'PRIMER_MIN_SIZE': 18,
                                  'PRIMER_MAX_SIZE': 25,
                                  'PRIMER_OPT_TM': 60.0,
                                  'PRIMER_MIN_TM': 57.0,
                                  'PRIMER_MAX_TM': 63.0,
                                  'PRIMER_MIN_GC': 20.0,
                                  'PRIMER_MAX_GC': 80.0,
                                  'PRIMER_MAX_POLY_X': 100,
                                  'PRIMER_INTERNAL_MAX_POLY_X': 100,
                                  'PRIMER_SALT_MONOVALENT': 50.0,
                                  'PRIMER_DNA_CONC': 50.0,
                                  'PRIMER_MAX_NS_ACCEPTED': 0,
                                  'PRIMER_MAX_SELF_ANY': 12,
                                  'PRIMER_MAX_SELF_END': 8,
                                  'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                                  'PRIMER_PAIR_MAX_COMPL_END': 8,
                                  'PRIMER_PRODUCT_SIZE_RANGE': [self.primer_sizes]})

    def get_seq(self):
        """
        :return: a string giving the chromosome sequence
        """
        path = self.data_dir + 'chr' + self.name + '.fasta'
        with open(path, "r") as f:  # open and read chromosome seq
            seq = f.read()

        seq = seq[seq.index('\n'):]  # cut off first line
        seq = seq.replace('\n', '')  # remove newline chars
        seq = seq.upper()  # make upper case
        return seq

    def get_genes_as_lists(self):
        """
        :return: a dict with keys as gene names, values as lists with [0][0] = a str telling the chr name, [0][1]
        = bool to indicate if its on the positive strand, and [1] as another list. [1][x] gives the chromosomal indices
        for an isoform of the gene. [1][x][0] and [1][x][1] are the start  and stop indices of a segment in that
        isoform. The dictionary returned will have all "unique" gene isoforms--unique in the sense that they all will
        have some exon  that isn't contained and in the same frame as another, longer isoform. However, nonunique
        isoforms will pass through if they have a segment covered from both sides by 2 alternative longer isoforms in
        the same frame. But this is not common, and not really a problem.
        """

        genes = {}  # the dictionary that will be returned.

        path = (self.data_dir + 'chr' + self.name + '_cds.fasta')  # this is how the gene data is named locally
        with open(path, "r") as f:  # open and read gene info line by line
            lines = f.readlines()  # get lines as list
        num_lines = len(lines)  # get number of lines

        initialize = True  # set to be true just for the starting case, and then permanently made false.
        curr_gene = ''  # a str to represent the gene currently being parsed
        curr_strand = None  # a bool to indicate the current gene's strand
        curr_inds = []  # a list of lists with each inner list being indices for a unique isoform of the gene.

        # the code below puts genes into its first state with long iso lists
        for i in range(num_lines):  # go through the lines of the file

            if lines[i][0] != '>':
                continue  # skip over lines that aren't headers

            if 'pseudo=true' in lines[i]:
                continue  # skip over pseudogenes

            found_gene = get_name(lines[i])  # get name as string through function call

            if curr_gene == found_gene:  # if not new, add the inds to curr_inds
                inds = get_idxs(lines[i])  #
                if inds != -1:  # if data on this gene is valid, add it
                    curr_inds.append(inds)  # append -1 to curr_inds to show the error
                continue  # skip to next line
            elif not initialize:  # if a new gene, add the last to the dict
                genes[curr_gene] = [curr_strand, curr_inds]
            else:
                initialize = False  # this is permanent

            # update curr_gene, curr_strand, and curr_inds with function calls
            curr_gene = found_gene
            curr_strand = get_strand(lines[i])
            inds = get_idxs(lines[i])
            if inds != -1:  # if indices not valid
                curr_inds = [inds]
            else:
                curr_inds = [[]]
        genes[curr_gene] = [curr_strand, curr_inds]  # get the last gene in the file after the loop is done

        # Now, we remove isoforms that are nonunique from each gene.
        # However, nonunique isoforms will pass through if they have a segment covered from both sides by 2 alternate
        # isoforms in the same frame.
        # Here, frames are relative, so we don't need to correct for strand
        for g in genes:

            # this section sorts the isoforms long to short
            iso_lens = []  # parallel to genes[g][1] and stores lens of the isos
            sorted_isos = []  # to build the list of sorted ones
            for iso in genes[g][1]:  # build the list of lengths
                iso_len = 0
                for seg in iso:  # calculate isoform length
                    iso_len += seg[1] - seg[0]
                iso_lens.append(iso_len)
            for i in range(len(iso_lens)):  # builds the iso list from the lens
                longest = iso_lens.index(max(iso_lens))
                sorted_isos.append(genes[g][1][longest])
                iso_lens[longest] = 0

            # this section creates a a parallel list of each segment's reading frames
            sorted_isos_frames = []  # parallel to sorted_isos; stores seg frames the number refers to the base in the
            # seg the first codon starts at
            for iso in sorted_isos:
                frames = [0]  # the first frame is always 0
                for i in range(1, len(iso)):  # iterate through segments; start at pair 1
                    frames.append((iso[i-1][1] - iso[i-1][0] - frames[i-1]) % 3)
                sorted_isos_frames.append(frames)  # add iso frames to the list

            # this section makes a final list without most nonunique isos (but some special cases could have one slip
            # through such as if an isoform seg is nonunique but only because it is covered on both sides by overlapping
            # segments of the same frame in 2 other isos)
            if len(sorted_isos) > 0:  # if this gene had valid isos and len > 0
                sorted_isos_final = [sorted_isos[0]]  # initialize list--will go into dict eventually
            else:  # if not, it's just an empty list
                sorted_isos_final = []
            for i in range(1, len(sorted_isos)):  # for each sorted iso from [1]
                nonunique = 0  # count of nonunique segments that tells to add or not
                for j in range(len(sorted_isos[i])):  # for each seg
                    for k in range(i):  # for each sorted isoform UP TO this one
                        keep_going = True  # for breaking out of the outer loop
                        for l in range(len(sorted_isos[k])):  # segments in prev isos
                            # if this segment is nonunique in scope and frame
                            if (sorted_isos[i][j][0] >= sorted_isos[k][l][0] and
                                    sorted_isos[i][j][1] <= sorted_isos[k][l][1] and
                                    sorted_isos_frames[i][j] ==
                                    sorted_isos_frames[k][l]):
                                nonunique += 1  # add count
                                keep_going = False  # to break outer loop
                                break  # to break inner loop
                        if not keep_going:  # to break out of the outer loop
                            break
                if nonunique < len(sorted_isos[i]):  # if a seg was unique
                    sorted_isos_final.append(sorted_isos[i])  # add to unique isos
            genes[g][1] = sorted_isos_final  # make this list of lists the [1] element in the dict value for this gene

        return genes

    def get_all_sites(self):
        """
        :return: first, a list of genomic edit sites with two inner lists. The first is for the + strand, and the second
        for the - strand. Second, a dictionary with genes as keys and a list of sites and values.
        Note: this function doesn't just take the positions of the last bases in isoforms. It only does so if it
        confirms that the last codon is a TAG. Data from genbank on some genes alleges that the coding sequence ends
        with something other than a stop codon.
        """

        # sites[0] is a list of target G locations in the genome in the + strand
        # sites[1] is a list of target C locations in the genome in the - strand
        sites = [[], []]

        # the keys are gene names and the values are lists of sites
        gene_sites = {}

        for gene_name, gene_list in self.genes.items():
            gene_sites[gene_name] = []  # initialize value list
            isos = gene_list[1]  # get isoform list

            if gene_list[0]:  # for + strand gene case
                for iso in isos:
                    if iso and self.seq[iso[-1][-1]-3: iso[-1][-1]] == 'TAG':  # if this is a TAG stop codon
                        gene_sites[gene_name].append(iso[-1][-1]-1)  # store G location
            else:  # for a - strand case
                for iso in isos:
                    if iso and self.seq[iso[0][0]: iso[0][0]+3] == 'CTA':  # if this is a TAG stop codon
                        gene_sites[gene_name].append(iso[0][0])  # store C location

            gene_sites[gene_name] = sorted(list(set(gene_sites[gene_name])))  # sort and remove duplicates, add to dict

            if gene_list[0]:  # add to genomic sites list
                sites[0] += gene_sites[gene_name]
            else:  # add to genomic sites list
                sites[1] += gene_sites[gene_name]

        sites[0] = sorted(list(set(sites[0])))  # sort and uniquify
        sites[1] = sorted(list(set(sites[1])))  # sort and uniquify

        return sites, gene_sites

    def get_editor_sites(self):
        """
        :return: lists of lists giving + and - genomic nuclease and base editor sites and dicts giving nuclease and base
        editor sites by gene
        This assumes an NG pam
        """

        be_sites = [{}, {}]  # keys to be indices and values to be cas9 targets + pams [0] is + sites and [1] is -

        for site in self.all_sites[0]:
            a_guides, c_guide = self.get_be_edit_guides(site, True)
            be_sites[0][site] = [a_guides, c_guide]

        for site in self.all_sites[1]:
            a_guides, c_guide = self.get_be_edit_guides(site, False)
            be_sites[1][site] = [a_guides, c_guide]

        # initialize dictionaries to store gene be and nuclease sites
        gene_be_sites = {}  # keys are gene names, values are dicts with index keys and guide pair list values

        # build the gene be sites and gene nuclease sites dicts from the gene_all_sites dict and the be_sites dict
        for gene in self.gene_all_sites:
            gene_be_sites[gene] = {}
            if self.genes[gene][0]:  # + strand case
                for site in self.gene_all_sites[gene]:
                    gene_be_sites[gene][site] = be_sites[0][site]
            else:
                for site in self.gene_all_sites[gene]:  # - strand case
                    gene_be_sites[gene][site] = be_sites[1][site]

        return be_sites, gene_be_sites

    def get_recoded_seq(self):
        """
        :return: str of the recoded seq of the GRCh38 version of the chromosome
        """
        # get the chromosome as a char array
        chromosome_array = array('B')
        chromosome_array.frombytes(self.seq.encode())

        for i in self.all_sites[0]:
            chromosome_array[i] = 65  # change to an A to make TAG->TAA
        for i in self.all_sites[1]:
            chromosome_array[i] = 84  # change complement to a T to make CTA->TTA

        # get the string back and trim off the b' and '
        seq = str(chromosome_array.tobytes())[2:-1]

        return seq

    def get_genes_as_objects(self):
        """
        :return: a dict with gene name keys and gene object values
        This is the final call of the init function. Prior to this call, gene info is stored in dicts with heterogeneous
        strings as values. This actually uses that info to make classes.
        """
        genes = {}  # to be filled and returned

        # fill a new dictionary with genes as objects instead of lists, passing all params into the Gene init function
        for name, gene_list in self.genes.items():
            genes[name] = Gene(self, name, gene_list[0], gene_list[1], 0, self.gene_all_sites[name],
                               self.gene_be_sites[name])
        return genes

    def get_multifunctional_sites(self):
        """
        Looks for edit sites that have another coding function somewhere else
        :return: a dict of int site keys and list values with [0] as the stop codon gene and [1] as the conflicting one
        """
        multifunctional = {}  # this will be filled with site keys and list values where the lists are gene names

        for name1, gene1 in self.genes.items():  # for each gene
            for site in gene1.all_sites:  # for each edit site in that gene
                for name2, gene2 in self.genes.items():  # for each gene
                    for iso in gene2.isos:  # for each iso in that second gene
                        for exon in iso:  # for each exon in that iso
                            if exon[0] <= site < exon[1]:  # if  site inside the exon
                                same_strand = gene1.strand == gene2.strand  # to help determine if both are stops
                                both_stops = False  # to be turned true if they are
                                if same_strand and gene1.strand and site == iso[-1][-1]-1:  # if both + stops
                                    both_stops = True
                                elif same_strand and not gene1.strand and site == iso[0][0]:  # if both - stops
                                    both_stops = True
                                if not both_stops:  # if this site is truly multifunctional
                                    multifunctional[site] = [name1, name2]  # add to dict
        return multifunctional

    def get_be_edit_guides(self, site, strand):
        """
        site is position of G to edit to an A
        strand is the strand as a bool
        Return a list where [0] is a list of ABE guides for daisy chaining, and [1] is a CBE guide
        # if this is not possible, then this return [-1, -1]
        """

        if strand:
            right_pam = []
            for pos in range(self.cbe_window[0], self.cbe_window[1]):  # for each position in the cbe_window
                if self.seq[site + pos - 22] == 'C':  # if the pam is correct (C complements G)
                    right_pam.append(pos)
            right_pam = order_around_mid(right_pam, self.c_mid)
            if right_pam:  # if a valid position exists
                position = right_pam[0]  # use most central site
                c_guide = rc(self.seq[site + position - 22: site + position])
                return [], c_guide

        else:
            right_pam = []
            for pos in range(self.cbe_window[0], self.cbe_window[1]):  # for each position in the cbe_window
                if self.seq[site - pos + 22] == 'G':
                    right_pam.append(pos)
            right_pam = order_around_mid(right_pam, self.c_mid)
            if right_pam:
                position = right_pam[0]  # use most central site
                c_guide = self.seq[site - position + 1: site - position + 23]
                return [], c_guide  # no ABE guides needed

        # executing the code below means that no PAM was found, and we have to to make a daisy chain
        a_sites = []
        if strand:
            for pos in range(self.cbe_window[0], self.cbe_window[1]):
                if self.seq[site + pos - 22] == 'T':  # if there's an A
                    a_sites.append(site+pos-22)
            mid_idx = site + self.c_mid - 22
        else:
            for pos in range(self.cbe_window[0], self.cbe_window[1]):  # for each position in the cbe_window
                if self.seq[site - pos + 22] == 'A':  # if there's an A
                    a_sites.append(site-pos+22)
            mid_idx = site - self.c_mid + 22

        if not a_sites:
            return -1, -1  # failed to find daisy chain-able sites

        else:
            # order the sites from middle outward because of higher fidelity in the middle of base editor windows
            a_sites = order_around_mid(a_sites, mid_idx)
            a_guides, a_idx = self.get_be_daisy_chain_guides(a_sites, strand)

            if a_guides == -1:  # if it is not possible ot daisy-chain edit
                return -1, -1
            else:
                if strand:  # add the last guide in the daisy chain with the altered PAM
                    c_guide = rc(self.seq[a_idx+1: a_idx+22]) + 'G'
                else:
                    c_guide = self.seq[a_idx-21: a_idx] + 'G'

                return a_guides, c_guide

    def get_be_daisy_chain_guides(self, sites, strand):
        """
        site is positions of As to edit
        strand is the strand as a bool
        """

        a_guides = []

        prevs = [-1 for i in range(len(sites))]  # make a parallel list to sites that gives the previous chain links

        site_idx = 0  # index
        for site in sites:
            # first, we want to see if this A can be made to a G using an ABE
            right_pam = []
            if strand:
                for pos in range(self.abe_window[0], self.abe_window[1]):  # for each position in the cbe_window
                    if self.seq[site + pos - 22] == 'C':  # if the pam is correct
                        right_pam.append(pos)
            else:
                for pos in range(self.abe_window[0], self.abe_window[1]):  # for each position in the cbe_window
                    if self.seq[site - pos + 22] == 'G':  # if the pam is correct
                        right_pam.append(pos)

            right_pam = order_around_mid(right_pam, self.a_mid)  # order middle outward for best fidelity

            if right_pam:  # if one or more found
                curr_idx = site_idx
                curr_pos = right_pam[0]

                while True:  # while tracing back the chain to the first A
                    # get the guide for each parent
                    if strand:
                        guide = rc(self.seq[sites[curr_idx]+curr_pos-22: sites[curr_idx]+curr_pos])
                    else:
                        guide = self.seq[sites[curr_idx]-curr_pos+1: sites[curr_idx]-curr_pos+23]
                    a_guides.append(guide)
                    if prevs[curr_idx] == -1:  # break when you get to the final site
                        break
                    last_idx = curr_idx
                    curr_idx = prevs[last_idx]
                    curr_pos = 22 - abs(sites[last_idx] - sites[curr_idx])

                for i in range(1, len(a_guides)):  # change A's in guides for G's (except for the first in the chain)
                    a_guides[i] = a_guides[i][:-1] + 'G'

                return a_guides, sites[curr_idx]

            else:  # if none were found, look for A's and add sites to queue
                a_sites = []
                if strand:
                    for pos in range(self.abe_window[0], self.abe_window[1]):  # for each position in the abe_window
                        if self.seq[site + pos - 22] == 'T':  # if and A
                            a_sites.append(site + pos - 22)
                    mid_idx = site + self.a_mid - 22
                else:
                    for pos in range(self.abe_window[0], self.abe_window[1]):  # for each position in the abe_window
                        if self.seq[site - pos + 22] == 'A':  # if an A
                            a_sites.append(site - pos + 22)
                    mid_idx = site - self.a_mid + 22
                sites += order_around_mid(a_sites, mid_idx)
                prevs += [site_idx for i in range(len(a_sites))]

            site_idx += 1  # increment index

        return -1, -1  # if no potential daisy chain scheme was found


class Gene:

    def __init__(self, chromosome, gene_name, strand, isos, active_iso, all_sites, be_sites):
        """
        :param chromosome: the paranet chromosome object in who's initialization this gene object is created
        :param name: string giving gene name
        :param strand: True if +, False if -
        :param isos: list of listpairs giving isoforms and exon starts and stops
        :param all_sites: list of all editing sites in this gene
        :param be_sites: dict of site keys and guide values
        """
        self.chromosome = chromosome  # the chromosome this object is part of
        self.name = gene_name  # name
        self.strand = strand  # True for + strand, False for -
        if self.name in self.chromosome.essentiality_data:  # if this gene in essentiality_data
            self.essentiality = self.chromosome.essentiality_data[self.name]  # get data as [num_ess, num_noness]
        else:  # if this gene in essentiality_data, usually because we're trying to find it with the wrong alias
            self.essentiality = 'unavailable'
        self.isos = isos  # the list of "unique" isoforms
        self.active_iso = active_iso  # the index of the iso to use for this gene (defaults at 0 changes if data erorr)
        self.all_sites = all_sites  # get all sites
        self.be_sites = be_sites  # get be edit sites
        self.cds = self.get_cds(recode=False)  # cds of isoform
        self.recoded_cds = self.get_cds(recode=True)  # recoded cds of isoform
        self.gene_region = self.get_gene_region(recode=False)  # get whole gene region
        self.recoded_gene_region = self.get_gene_region(recode=True)  # get whole gene region, recoded

    def get_cds(self, iso=0, recode=False):
        """
        :param iso: which isoform (when sorted longest to shortest) to return the cds of
        :param recode: bool telling whether or not to return the recoded cds
        :return: a string giving the coding sequence
        """
        if len(self.isos[0]) == 0:  # if no indices available, return empty str
            return ''

        segments = self.isos[iso]  # get segment list
        cds = ''  # initialize the string to build

        if not recode:  # build nonrecoded sequence
            for seg in segments:
                cds += self.chromosome.seq[seg[0]:seg[1]]
        else:  # build recoded sequence
            for seg in segments:
                cds += self.chromosome.recoded_seq[seg[0]:seg[1]]

        if len(cds) % 3 != 0 and iso < len(self.isos)-1:  # if error and cds wrong length, try to get another iso
            return self.get_cds(iso=iso+1, recode=recode)

        elif len(cds) % 3 != 0:  # if error and cds wrong length and no alternatives, return empty string
            self.active_iso = None
            return ''

        if not self.strand:  # rc is neg strand
            cds = rc(cds)

        self.active_iso = iso

        return cds

    def get_gene_region(self, recode=False, flank=0):
        """
        :param recode: bool which tells whether or not to return the recoded gene
        :param flank: int which tells how many extra bases to return flanking the region on either side
        :return: a str giving the whole gene region from first start codon to last stop codon w/ flanking bases
        """

        if len(self.isos[0]) == 0 or self.active_iso is None:  # if no isos for gene because of gbk data error
            return ''

        leftmost = self.isos[0][0][0]  # 5' most index for an exon of this gene
        rightmost = self.isos[0][-1][-1]  # 3' most index for an exon of this gene
        for iso in self.isos:  # find leftmost and rightmost indices
            if iso[0][0] < leftmost:
                leftmost = iso[0][0]
            if iso[-1][-1] > rightmost:
                rightmost = iso[-1][-1]

        # add flank
        leftmost -= flank
        rightmost += flank

        if recode:
            if self.strand:
                return self.chromosome.recoded_seq[leftmost: rightmost]
            else:
                return rc(self.chromosome.recoded_seq[leftmost: rightmost])
        else:
            if self.strand:
                return self.chromosome.seq[leftmost: rightmost]
            else:
                return rc(self.chromosome.seq[leftmost: rightmost])

    def get_primer_info(self, window_size=350):
        """
        :param window_size: int giving how big a window to look for primers in
        :return: dict whose keys are edit sites and whose values are dictionaries that contain the
        output of primer3's designPrimers and for pairs of primers, the left, right, and product length in list form.
        Note that this code doesn't guarantee that primers won't match to parts of the recoded gene with silent
        mutations from nuclease/HR editing.
        """

        site_primers = {}  # dict to return with site keys and primer-dict values
        mid = round(window_size/2)  # precalculate half window_size

        for site in self.all_sites:
            if self.strand:  # + strand
                search_window = self.chromosome.seq[site-mid: site+mid]
            else:  # - strand
                search_window = rc(self.chromosome.seq[site-mid: site+mid])

            # get primer3 results on window
            primer_results = p3.bindings.designPrimers({'SEQUENCE_ID': self.name+'_site_'+str(site),
                                                        'SEQUENCE_TEMPLATE': search_window,
                                                        'SEQUENCE_INCLUDED_REGION': [mid-150, 300]})  # start, length

            # initialize site_dict and get the number of primer pairs
            site_dict = {'WINDOW_SIZE': window_size, 'P3_RESULTS': primer_results}
            pair_num = primer_results['PRIMER_PAIR_NUM_RETURNED']

            for j in range(pair_num):  # for each primer pair
                j_str = str(j)  # string giving the number of the pair
                left = primer_results['PRIMER_LEFT_'+j_str+'_SEQUENCE']  # left pair
                right = primer_results['PRIMER_RIGHT_'+j_str+'_SEQUENCE']  # right pair
                product_len = primer_results['PRIMER_PAIR_'+j_str+'_PRODUCT_SIZE']  # product length
                site_dict['PAIR_'+j_str] = {'LEFT': left, 'RIGHT': right, 'PRODUCT_LEN': product_len,
                                            'REFERENCE': search_window}  # entry for pair

            site_primers[site] = site_dict  # add entry for this site

        return site_primers