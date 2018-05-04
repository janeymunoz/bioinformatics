def mk_rna_aa_dict():
    ''' Returns a dictinoary of RNA codons to amino acids.
    key = 3 RNA nucleotides
    value = corresponding amino acid
    '''
    from biotools_getdata import open_read
    data = open_read('rna_codon_table.txt')
    rna_codon_dict = {}
    for i in range(0, len(data)):
        line = data[i].split()
        for j in range(0, len(line), 2):
            rna_codon_dict[line[j]] = line[j + 1]
    return rna_codon_dict


def mk_peptide_mass_dict():
    ''' Returns a dictionary of peptides and their mass value in Daltons. '''
    from biotools_getdata import open_read
    data = open_read('monoisotopic_mass_table.txt')
    monomass_dict = {}
    for i in range(0, len(data)):
        cur_pair = data[i][:-1]
        pair_key = cur_pair[:1]
        pair_val = float(cur_pair[4:])
        monomass_dict[pair_key] = pair_val
    return monomass_dict
