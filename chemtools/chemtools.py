def parse_mol(mol_str):
    ''' Parses a string representing a molecule into symbols and subscripts.
    Takes one string as an argument. Must be formatted appropriately:
        - no spaces.
        - elements are represented per IUPAC standards (camel-case).
        - subscripts for elements directly follow the element.
    Returns two lists as tuple: one with a list of element symbols,
    other with floats representing the element subscripts.
    '''
    len_str = len(mol_str)
    subscripts = []
    atom_syms = []
    for i in range(0, len_str):
        curr_char = mol_str[i]
        # Compare adjacent chars. in molecule, excluding last char.
        if i < len_str - 1:
            next_char = mol_str[i + 1]
            # Append symbols and subscripts to list based on rules.
            if curr_char.isupper() and next_char.isupper():
                atom_syms.append(curr_char)
                subscripts.append(float(1))
            elif curr_char.isupper() and next_char.islower():
                atom_syms.append(curr_char + next_char)
            elif curr_char.islower() and next_char.isupper():
                subscripts.append(float(1))
            elif curr_char.isupper() and next_char.isdigit():
                atom_syms.append(curr_char)
            # If char. is digit, check following chars. for add'l digits.
            elif curr_char.isdigit() and mol_str[i - 1].isalpha():
                subscript = ""
                for j in range(i, len_str):
                    addl_sub = mol_str[j]
                    if addl_sub.isdigit() or addl_sub == ".":
                        subscript += addl_sub
                        if j == len_str - 1:
                            subscripts.append(float(subscript))
                        else:
                            continue
                    else:
                        subscripts.append(float(subscript))
                        break
            else:
                continue
        # Check final string char. when reached in for loop.
        else:
            prev_char = mol_str[i - 1]
            if curr_char.islower():
                subscripts.append(float(1))
            elif curr_char.isupper():
                atom_syms.append(curr_char)
                subscripts.append(float(1))
            elif curr_char.isdigit() and prev_char.isalpha():
                subscripts.append(float(curr_char))
    return atom_syms, subscripts


def parse_atom_data(return_type):
    ''' Returns a dictionary of atomic weights, or atomic numbers.
    Reads a text file containing atom symbols, weights, and numbers, and
    returns a dictionary of weights or numbers based on return_type arg,
    "weight" or "number".
    '''
    # Read atomic weights/number data from text file.
    import os
    file_dir = os.path.dirname(__file__)
    rel_path = 'ref_tables/atomic_weights.txt'
    abs_path = os.path.join(file_dir, rel_path)
    f = open(abs_path, "r")
    data = f.readlines()
    # Parse data and create key/values based on return_type.
    atom_dict = {}
    for info in data:
        indices = []
        for i, j in enumerate(info):
            if j == "\t":
                indices.append(i)
            elif j == "\n":
                indices.append(i)
            else:
                continue
        index1 = indices[0]
        index2 = indices[1]
        index3 = indices[2]
        atom_name = info[:index1]
        if return_type == "weight":
            atom_weight = info[index1 + 1:index2]
            atom_dict[atom_name] = atom_weight
        elif return_type == "number":
            atom_num = info[index2 + 1:index3]
            atom_dict[atom_name] = atom_num
    return atom_dict


def get_atom_info(atom_sym, return_type):
    ''' Takes an atom symbol and return_type, returns atomic weight or number.
    Atom symbol should be standard IUPAC element style (camel-case), and
    return_type is "weight" or "number".
    '''
    atom_dict = parse_atom_data(return_type)
    info = atom_dict.get(atom_sym.title())
    return info


def molecule_wt(mol_str):
    '''Calculates total molecular weight of a molecule.
    Argument must be a single string, no spaces, with elements represented in
    standard (camel-case) IUPAC format. A coefficient may precede the molecule.
    Examples of appropriate formatting:
        "2C6H12O6" (i.e., 2 glucose molecules)
        "CH1.08O0.5N0.2" (gen. formula for dry cell biomass, float-type subs.)
    '''
    # Get atomic weight dictionary.
    atom_wt = parse_atom_data("weight")
    # Find number of characters in mol_str
    len_str = len(mol_str)
    # Check for presence of coefficient and differentiate coeff. from molecule.
    if mol_str[0].isalpha():
        single_molecule = mol_str
        coeff = 1
    else:
        coeff = ""
        for i in range(0, len_str):
            char = mol_str[i]
            if char.isdigit() or char == ".":
                coeff += char
            else:
                coeff = float(coeff)
                single_molecule = mol_str[i:]
                break
    # Use parse_mol func. to put molecule into sub-parts (elems and subscrpts).
    mol_data = parse_mol(single_molecule)
    atom_sym = mol_data[0]
    atom_subs = mol_data[1]
    num_elems = len(atom_sym)
    # Calculate total molecular weight.
    tot_wt = 0
    for i in range(0, num_elems):
        tot_wt += float(atom_wt.get(atom_sym[i])) * atom_subs[i]
    tot_wt *= coeff
    return tot_wt
