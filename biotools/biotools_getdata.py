def open_read(filename):
    ''' Reads a text file and stores each line as an item in a list. '''
    with open(filename, 'r') as f:
        data = f.read().splitlines()
    return data


def open_txt_1str(filename):
    ''' Reads text file, returns one continuous string without linebreaks. '''
    with open(filename, 'r') as f:
        tot_str = f.read().replace("\n", "")
    return tot_str


def open_fasta(filename, return_type):
    ''' Reads a text file in FASTA format and returns FASTA IDs and/or sequences.
    Returned values are based on return_type argument: return_type is either
    "fasta", "seqs", "both".  Returned values are lists, or a tuple of lists.
    '''
    fasta_ids = []
    seq_strs = []
    seq_str = ""
    fasta_count = 0
    seq_count = 0
    with open(filename, 'r') as f:
        for line in f:
            cur_line = line[:-1]
            if "\n" in cur_line:
                cur_line = cur_line[:-1]
            else:
                if ">" in cur_line:
                    fasta_ids.append(cur_line[1:])
                    fasta_count += 1
                else:
                    if fasta_count == seq_count + 1:
                        seq_str += cur_line
                    else:
                        seq_strs.append(seq_str)
                        seq_count += 1
                        seq_str = ""
                        seq_str += cur_line
        seq_strs.append(seq_str)
    if return_type.lower() == "fasta":
        return fasta_ids
    elif return_type.lower() == "seq":
        return seq_strs
    elif return_type.lower() == "both":
        return fasta_ids, seq_strs
    else:
        return "invalid return type"


def parse_fasta(data, return_type):
    ''' Takes a string in FASTA format and returns FASTA IDs and/or sequences.
    Returned values are based on return_type argument: return_type is either
    "fasta", "seqs", "both".  Returned values are lists, or a tuple of lists.
    '''
    fasta_ids = []
    seq_strs = []
    seq_str = ""
    fasta_count = 0
    seq_count = 0
    data_by_line = data.splitlines()
    for line in data_by_line:
        cur_line = line
        if "\n" in cur_line:
            cur_line = cur_line[:-1]
        else:
            if ">" in cur_line:
                fasta_ids.append(cur_line[1:])
                fasta_count += 1
            else:
                if fasta_count == seq_count + 1:
                    seq_str += cur_line
                else:
                    seq_strs.append(seq_str)
                    seq_count += 1
                    seq_str = ""
                    seq_str += cur_line
    seq_strs.append(seq_str)
    if return_type.lower() == "fasta":
        return fasta_ids
    elif return_type.lower() == "seq":
        return seq_strs
    elif return_type.lower() == "both":
        return fasta_ids, seq_strs
    else:
        return "invalid return type"


def uniprot_get(uni_id):
    ''' Gets data on proteins from the uniprot.org website.
    Given a uni_id, returns in FASTA format.
    '''
    import requests
    data = requests.get('http://uniprot.org/uniprot/' + uni_id + '.fasta')
    fasta_data = data.text
    return fasta_data


def check_seq(seq, str_type):
    ''' Checks if all characters are valid DNA or RNA bases in string.
    Will check for DNA if seq_type="DNA" or RNA if seq_type="RNA".
    '''
    # Get relevant alphabet based on seq_type.
    if str_type.upper() == "DNA":
        check = set("ATGC")
    elif str_type.upper() == "RNA":
        check = set("AUGC")
    else:
        raise ValueError('the seq_type arg. must be either "RNA" or "DNA".')
    # Check seq against check alphabet.
    leftover = set(seq.upper()) - check
    return not leftover


def walk_str(str, block_size, offset):
    ''' Walks a string, returns a list of substrings.
    Size of substring and offset of walk must be specified as arguments.
    '''
    len_str = len(str)
    offset_end = len_str - block_size
    str_blocks = []
    if block_size > len_str:
        raise ValueError("length of string must be greater than block size.")
    elif block_size == len_str:
        str_blocks.append(str)
        return str_blocks
    else:
        for i in range(0, offset_end + 1, offset):
            mini_str = str[i: i + block_size]
            str_blocks.append(mini_str)
    return str_blocks


def list_to_matr(whol_list, num_rows):
    ''' Takes a list of integers and returns a matrix (lists within a list).
    The number of rows must be specified as an argument, and the length of the
    whole list modulo the number of rows must be 0.
    '''
    len_list = len(whol_list)
    nums_in_row = len_list // num_rows
    matr_list = []
    if num_rows < 1:
        raise ValueError("the number of rows must be greater than 0.")
    elif num_rows == 1:
        return whol_list
    elif len_list % num_rows != 0:
        raise ValueError("length of list modulo number of rows must be 0.")
    else:
        for i in range(0, num_rows):
            sub_list = "sub_" + str(i)
            matr_list.append(sub_list)
            start_ind = int(nums_in_row * i)
            end_ind = int(nums_in_row * i + nums_in_row)
            matr_list[i] = whol_list[start_ind:end_ind]
            print(matr_list[i])
    return matr_list
