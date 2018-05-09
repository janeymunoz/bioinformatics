def dna_to_rna(dna_str):
    ''' Takes a DNA string, returns corresponding RNA string. '''
    from biotools_getdata import check_seq
    if check_seq(dna_str, "DNA"):
        rna_str = (''.join(dna_str.upper()).replace("T", "U"))
        return rna_str
    else:
        raise ValueError('dna_str has invalid characters, not in set("ATGC"))')


def rev_compl(str, type):
    ''' Takes a DNA or RNA string and returns its reverse complement.
    DNA or RNA type must be specified in argument as "dna" or "rna".
    '''
    from biotools_getdata import check_seq
    len_str = len(str)
    # Check if string characters are of DNA or RNA set.
    if type.upper() != "DNA" and type.upper() != "RNA":
        raise ValueError('type argument must be either "DNA" or "RNA".')
    if check_seq(str, type):
        s_c = ""
        for i in range(len_str - 1, -1, -1):
            curr_base = str[i]
            new_base = ""
            if curr_base.upper() == "G":
                new_base = "C"
            elif curr_base.upper() == "C":
                new_base = "G"
            elif curr_base.upper() == "A":
                if type.upper() == "DNA":
                    new_base = "T"
                elif type.upper() == "RNA":
                    new_base = "U"
            else:
                new_base = "A"
            s_c += new_base
        return s_c
    else:
        rna_set = 'RNA set("AUGC")'
        dna_set = 'DNA set("ATGC")'
        if type.upper() == "RNA":
            raise ValueError('str has invalid characters, not in ' + rna_set)
        else:
            raise ValueError('str has invalid characters, not in ' + dna_set)


def gc_content(str):
    ''' Calculates the GC content of a given DNA  or RNA string. '''
    from biotools_getdata import check_seq
    if check_seq(str, "rna") or check_seq(str, "dna"):
        count_gc = str.upper().count("G") + str.count("C")
        perc_gc = round((100 * count_gc / len(str)), 2)
        return perc_gc
    else:
        raise ValueError('str has invalid characters (not dna or rna bases).')


def single_rna_to_aa(rna_str):
    ''' Takes a 3-character rna string and returns its amino acid. '''
    from biotools_refdicts import mk_rna_aa_dict
    from biotools_getdata import check_seq
    if len(rna_str) > 3:
        raise ValueError('function takes a 3-character string as argument.')
    if check_seq(rna_str, "rna"):
        rna_codon_dict = mk_rna_aa_dict()
        protein = rna_codon_dict.get(rna_str.upper())
        return protein
    else:
        raise ValueError('arg. has invalid characters, not in set("AUGC").')


def rna_to_prot(rna_str):
    ''' Takes a string of RNA and returns a protein string.
    Protein string terminates, and is returned when "Stop" codon occurs.
    '''
    from biotools_refdicts import mk_rna_aa_dict
    from biotools_getdata import walk_str, check_seq
    if check_seq(rna_str, "rna"):
        rna_codon_dict = mk_rna_aa_dict()
        rna_codon_list = walk_str(rna_str, 3, 3)
        protein_str = ""
        for i in range(0, len(rna_codon_list)):
            protein_str0 = rna_codon_dict.get(rna_codon_list[i].upper())
            if protein_str0 == "Stop":
                return protein_str
            else:
                protein_str += protein_str0
    else:
        raise ValueError('arg. has invalid characters, not in set("AUGC").')


def all_rna_to_prot(rna_str):
    ''' Takes a string of RNA and returns a protein string.
    Returns total protein string, including all "Stop" codons.
    '''
    from biotools_refdicts import mk_rna_aa_dict
    from biotools_getdata import walk_str, check_seq
    if check_seq(rna_str, "rna"):
        rna_codon_dict = mk_rna_aa_dict()
        rna_codon_list = walk_str(rna_str, 3, 3)
        protein_str = ""
        for i in range(0, len(rna_codon_list)):
            protein_str += rna_codon_dict.get(rna_codon_list[i])
        return protein_str
    else:
        raise ValueError('arg. has invalid characters, not in set("AUGC").')


def splice(s, introns):
    ''' Takes a DNA string, s, and a list of its introns. Returns protein string.
    DNA with introns and exons --> DNA with exons --> RNA --> Protein
    '''
    from biotools_getdata import check_seq
    if check_seq(s, "dna"):
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
    else:
        raise ValueError('arg. has invalid characters, not in set("AUGC").')


def trans_ratio(s1, s2):
    ''' Takes two DNA strings of equal length and returns Ti/Tv ratio.
    Where Ti is the number of transitions (C <--> T or A <--> G)
    and Tv is the number of transversions (A <--> C or G <--> T).
    '''
    from biotools_getdata import check_seq
    pur = ["A", "G"]
    pyr = ["C", "T"]
    transitions = 0
    transversions = 0
    if not check_seq(s1, "dna") or not check_seq(s2, "dna") \
            or len(s1) != len(s2):
        raise ValueError("arguments must be two DNA strings of equal length.")
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
