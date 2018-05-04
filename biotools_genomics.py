def dna_to_rna(dna_str):
    ''' Takes a DNA string, returns corresponding RNA string. '''
    rna_str = (''.join(dna_str.upper()).replace("T", "U"))
    return rna_str


def rev_compl(str, type):
    ''' Takes a DNA or RNA string and returns its reverse complement.
    DNA or RNA type must be specified in argument as "dna" or "rna".
    '''
    len_str = len(str)
    s_c = ""
    for i in range(len_str - 1, -1, -1):
        curr_base = str[i]
        new_base = ""
        if curr_base.upper() == "G":
            new_base = "C"
        elif curr_base.upper() == "C":
            new_base = "G"
        elif curr_base.upper() == "A":
            if type.lower() == "dna":
                new_base = "T"
            elif type.lower() == "rna":
                new_base = "U"
        else:
            new_base = "A"
        s_c += new_base
    return s_c


def gc_content(dna_str):
    ''' Calculates the GC content of a given DNA string. '''
    count_gc = dna_str.upper().count("G") + dna_str.count("C")
    perc_gc = round((100 * count_gc / len(dna_str)), 2)
    return perc_gc


def single_rna_to_aa(rna_str):
    ''' Takes a 3-character rna string and returns its amino acid. '''
    from biotools_refdicts import mk_rna_aa_dict
    if not isinstance(rna_str, str) or len(rna_str) != 3:
        return "Arg error: takes one 3-character RNA string as an argument"
    else:
        rna_codon_dict = mk_rna_aa_dict()
        protein = rna_codon_dict.get(rna_str.upper())
        return protein


def rna_to_prot(rna_str):
    ''' Takes a string of RNA and returns a protein string.
    Protein string terminates, and is returned when "Stop" codon occurs.
    '''
    from biotools_refdicts import mk_rna_aa_dict
    from biotools_getdata import walk_str
    rna_codon_dict = mk_rna_aa_dict()
    rna_codon_list = walk_str(rna_str, 3, 3)
    protein_str = ""
    for i in range(0, len(rna_codon_list)):
        protein_str0 = rna_codon_dict.get(rna_codon_list[i].upper())
        if protein_str0 == "Stop":
            return protein_str
        else:
            protein_str += protein_str0


def all_rna_to_prot(rna_str):
    ''' Takes a string of RNA and returns a protein string.
    Returns total protein string, including all "Stop" codons.
    '''
    from biotools_refdicts import mk_rna_aa_dict
    from biotools_getdata import walk_str
    rna_codon_dict = mk_rna_aa_dict()
    rna_codon_list = walk_str(rna_str, 3, 3)
    protein_str = ""
    for i in range(0, len(rna_codon_list)):
        protein_str += rna_codon_dict.get(rna_codon_list[i])
    return protein_str


def splice(s, introns):
    ''' Takes a DNA string, s, and a list of its introns. Returns protein string.
    DNA with introns and exons --> DNA with exons --> RNA --> Protein
    '''
    for intron in introns:
        check_intron = s.count(intron)
        if check_intron == 0:
            continue
        elif check_intron == 1:
            s = s.replace(intron, "")
        # if intron is found more than once, remove all instances
        else:
            for instance in range(check_intron):
                s = s.replace(intron, "")
    rna_s = dna_to_rna(s)
    prot_str = rna_to_prot(rna_s)
    return prot_str


def trans_ratio(s1, s2):
    ''' Takes two DNA strings of equal length and returns Ti/Tv ratio.
    Where Ti is the number of transitions (C <--> T or A <--> G)
    and Tv is the number of transversions (A <--> C or G <--> T).
    '''
    pur = ["A", "G"]
    pyr = ["C", "T"]
    transitions = 0
    transversions = 0
    if not s1.isalpha() or not s2.isalpha() or len(s1) != len(s2):
        return "arguments must be two DNA strings of equal length"
    else:
        len_str = len(s1)
        for i in range(0, len_str):
            base1 = s1[i].upper()
            base2 = s2[i].upper()
            if base1 == base2:
                continue
            else:
                if base1 in pur and base2 in pur \
                  or base1 in pyr and base2 in pyr:
                    transitions += 1
                else:
                    transversions += 1
        if transversions > 0:
            ratio = transitions / transversions
            print("Number of transitions:", transitions)
            print("Number of transversions:", transversions)
            print("Ratio of transitions to transversions is:", ratio)
            return ratio
        else:
            print("Number of transistions:", transitions)
            print("Number of transversions: 0")
            print("Ratio of transistions to transversions not calculatable.")
            return transitions, transversions
