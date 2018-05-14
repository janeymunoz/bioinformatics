def dna_to_rna(dna_seq):
    '''Takes a DNA string, returns corresponding RNA string. '''
    from biotools_getdata import check_seq
    if check_seq(dna_seq, "DNA"):
        rna_seq = (''.join(dna_seq.upper()).replace("T", "U"))
        return rna_seq
    else:
        raise ValueError('dna_str has invalid characters, not in set("ATGC"))')


def rev_compl(seq, seq_type):
    '''Takes a DNA or RNA string and returns its reverse complement.
    DNA or RNA type must be specified in argument as "dna" or "rna".
    '''
    from biotools_getdata import check_seq
    len_seq = len(seq)
    # Check if string characters are of DNA or RNA set.
    if seq_type.upper() != "DNA" and seq_type.upper() != "RNA":
        raise ValueError('type argument must be either "DNA" or "RNA".')
    if check_seq(seq, seq_type):
        s_c = ""
        for i in range(len_seq - 1, -1, -1):
            curr_base = seq[i]
            new_base = ""
            if curr_base.upper() == "G":
                new_base = "C"
            elif curr_base.upper() == "C":
                new_base = "G"
            elif curr_base.upper() == "A":
                if seq_type.upper() == "DNA":
                    new_base = "T"
                elif seq_type.upper() == "RNA":
                    new_base = "U"
            else:
                new_base = "A"
            s_c += new_base
        return s_c
    else:
        rna_set = 'RNA set("AUGC")'
        dna_set = 'DNA set("ATGC")'
        if seq_type.upper() == "RNA":
            raise ValueError('str has invalid characters, not in ' + rna_set)
        else:
            raise ValueError('str has invalid characters, not in ' + dna_set)


def gc_content(seq):
    '''Calculates the GC content of a given DNA  or RNA string. '''
    from biotools_getdata import check_seq
    if check_seq(seq, "rna") or check_seq(seq, "dna"):
        count_gc = seq.upper().count("G") + seq.count("C")
        perc_gc = round((100 * count_gc / len(seq)), 2)
        return perc_gc
    else:
        raise ValueError('str has invalid characters (not dna or rna bases).')


def single_rna_to_aa(rna_seq):
    '''Takes a 3-character rna string and returns its amino acid. '''
    from biotools_refdicts import mk_rna_aa_dict
    from biotools_getdata import check_seq
    if len(rna_seq) > 3:
        raise ValueError('function takes a 3-character string as argument.')
    if check_seq(rna_seq, "rna"):
        rna_codon_dict = mk_rna_aa_dict()
        protein = rna_codon_dict.get(rna_seq.upper())
        return protein
    else:
        raise ValueError('arg. has invalid characters, not in set("AUGC").')


def rna_to_prot(rna_seq):
    '''Takes a string of RNA and returns a protein string.
    Protein string terminates, and is returned when "Stop" codon occurs.
    '''
    from biotools_refdicts import mk_rna_aa_dict
    from biotools_getdata import walk_str, check_seq
    if check_seq(rna_seq, "rna"):
        rna_codon_dict = mk_rna_aa_dict()
        rna_codon_list = walk_str(rna_seq, 3, 3)
        protein_str = ""
        for i in range(0, len(rna_codon_list)):
            protein_str0 = rna_codon_dict.get(rna_codon_list[i].upper())
            if protein_str0 == "Stop":
                return protein_str
            else:
                protein_str += protein_str0
    else:
        raise ValueError('arg. has invalid characters, not in set("AUGC").')


def all_rna_to_prot(rna_seq):
    '''Takes a string of RNA and returns a protein string.
    Returns total protein string, including all "Stop" codons.
    '''
    from biotools_refdicts import mk_rna_aa_dict
    from biotools_getdata import walk_str, check_seq
    if check_seq(rna_seq, "rna"):
        rna_codon_dict = mk_rna_aa_dict()
        rna_codon_list = walk_str(rna_seq, 3, 3)
        protein_str = ""
        for i in range(0, len(rna_codon_list)):
            protein_str += rna_codon_dict.get(rna_codon_list[i])
        return protein_str
    else:
        raise ValueError('arg. has invalid characters, not in set("AUGC").')


def splice(s, introns):
    '''Takes a DNA string, s, and a list of its introns. Returns protein string.
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
    '''Takes two DNA strings of equal length and returns Ti/Tv ratio.
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


def count_motif(seq, check_motif):
    '''Counts the number of times a motif appears in a sequence.'''
    from biotools_getdata import walk_str
    # Count number of times check_motif found and log indices.
    len_motif = len(check_motif)
    motifs = walk_str(seq, len_motif, 1)
    num_motifs = len(motifs)
    count = 0
    indx = []
    for i in range(0, num_motifs):
        motif = motifs[i]
        if motif == check_motif:
            count += 1
            indx.append(i)
        else:
            continue
    return count, indx


def count_motif_fr(seq, check_motif):
    '''Checks a sequence for a motif in forwards and reverse direction .'''
    # Forward check.
    for_data = count_motif(seq, check_motif)
    for_count = for_data[0]
    for_indx = for_data[1]
    # Reverse check.
    rev_seq = seq[::-1]
    rev_data = count_motif(rev_seq, check_motif)
    rev_count = rev_data[0]
    rev_indx = rev_data[1]
    # Return counts.
    return for_count, for_indx, rev_count, rev_indx


def get_rfs(seq, seq_type):
    '''Get reading frames of a sequence.
    Returns a dicitonary of the amino acid codons of 6 reading frames,
    3 for sequence as is (forward), and 3 for its reverse complement.
    Keys in dictionary are "rf_0", "rf_1", etc.
    '''
    from biotools_getdata import walk_str
    all_rf = {}
    keys = []
    rf_offset_list = [None, -2, -1]
    len_codon = 3
    if seq_type.upper() == "DNA" or seq_type.upper() == "RNA":
        # For sequence in "forward" direction.
        for i in range(0, 3):
            str_i = str(i)
            rf_offset = rf_offset_list[i]
            key = "rf_" + str_i
            keys.append(key)
            rf_str = seq[i:rf_offset]
            # Get list of codons at this reading frame and add to dictionary.
            all_rf[key] = walk_str(rf_str, len_codon, len_codon)
        # For sequence's reverse complement.
        rev_c = rev_compl(seq, seq_type)
        for i in range(0, 3):
            str_i = str(i + 3)
            rf_offset = rf_offset_list[i]
            key = "rf_" + str_i
            keys.append(key)
            rf_str = rev_c[i:rf_offset]
            # Get list of codons at this reading frame and add to dictionary.
            all_rf[key] = walk_str(rf_str, len_codon, len_codon)
    else:
        raise ValueError('"type" argument must be either "DNA" or "RNA"')
    return all_rf, keys


def get_orf(seq, seq_type):
    '''Takes a DNA or RNA sequence and returns all candidate protein strings
    in a list that can be translated from open reading frames (ORFs) of the
    sequence.
    '''
    # Get DNA to RNA converter function, dna_to_rna. Needed if type == "DNA".
    # Get dictionary with RNA seqs as keys and their amino acids as values.
    from biotools_refdicts import mk_rna_aa_dict
    rna_aa_dict = mk_rna_aa_dict()
    # Use get_rfs() function to get dictionary of all potential reading frame
    # codon combinations. Keys are "rf_#" and values are lists of codons.
    rf_data = get_rfs(seq, seq_type)
    all_rf = rf_data[0]
    keys = rf_data[1]
    all_prot_str = []
    # If type == "DNA", must convert sequence to RNA first.
    if seq_type.upper() == "DNA":
        for i in range(0,  6):
            key = keys[i]
            rf = all_rf[key]
            num_dna_trip = len(rf)
            for j in range(0, num_dna_trip):
                dna_trip = rf[j]
                codon = dna_to_rna(dna_trip)
                # Replacen DNA nuecleotide triplet with RNA codon.
                rf[j] = codon
    # Check all six reading frames, now with RNA codons in groups, for ORFs.
    for i in range(0, 6):
        key = keys[i]
        rf = all_rf[key]
        num_codons = len(rf)
        for j in range(0, num_codons):
            codon = rf[j]
            aa = rna_aa_dict[codon]
            prot_str = ""
            if aa == "M":
                for k in range(j, num_codons):
                    codon = rf[k]
                    aa = rna_aa_dict[codon]
                    if aa != "Stop":
                        prot_str += aa
                    else:
                        # Append strings to list, but not duplicates.
                        if prot_str not in all_prot_str:
                            all_prot_str.append(prot_str)
                        else:
                            continue
                        break
            else:
                continue
    return all_prot_str


def get_dna_orf(dna_seq):
    '''Takes a DNA string and returns all candidate protein strings.
    More efficient than the "get_orf" funciton, but less flexible in inputs.
    '''
    from biotools_getdata import walk_str
    from biotools_refdicts import mk_rna_aa_dict
    # Get RNA sequence from input DNA sequence.
    rna_seq = dna_to_rna(dna_seq)
    # Get reverse complement of RNA sequence.
    rev_c = rev_compl(rna_seq, "RNA")
    # Get RNA codon to amino acid dictionary.
    rna_aa_dict = mk_rna_aa_dict()
    # Store details for parsing RNA sequences.
    rf_offset_list = [None, -2, -1, None, -2, -1]
    len_codon = 3
    # Container for  all candidate protein strings.
    all_prot_str = []
    # Get all reading frames (w/ codons grouped).
    for i in range(0, 6):
        rf_offset = rf_offset_list[i]
        # For sequence in "forward" direciton.
        if i < 3:
            rf_str = rna_seq[i:rf_offset]
        # For reverse complement reading frames.
        else:
            rf_str = rev_c[(i - 3):rf_offset]
        # Get list of codons at this reading frame and add to dictionary.
        rf_codons = walk_str(rf_str, len_codon, len_codon)
        num_codons = len(rf_codons)
        # Get all candidate protein strings.
        for j in range(0, num_codons):
            codon = rf_codons[j]
            aa = rna_aa_dict[codon]
            prot_str = ""
            if aa == "M":
                for k in range(j, num_codons):
                    codon = rf_codons[k]
                    aa = rna_aa_dict[codon]
                    if aa != "Stop":
                        prot_str += aa
                    else:
                        # Append strings to list, but not duplicates.
                        if prot_str not in all_prot_str:
                            all_prot_str.append(prot_str)
                        else:
                            continue
                        break
            else:
                continue
    return all_prot_str
