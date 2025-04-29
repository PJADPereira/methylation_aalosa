#!/usr/bin/env python3

import sys
import re
import gzip
import argparse


class Enzyme(object):
    '''
    Object to represent an enzyme and its function
    On init requires a name (name of the enzyme) and a sequence (recognition site with cut site, provide the 5'-3' strand)
    Accepts IUPAC nucleotide codes in the recognition site (except gaps for obvious reasons)
    The insilico digestion is not sensitive to any of the methylation statuss
    '''
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.offset = self.seq.find("^")
        self.site = self.seq.replace("^", "")
        self.sequence_code = {"B": "[CGT]",
                              "D": "[AGT]",
                              "H": "[ACT]",
                              "K": "[GT]",
                              "M": "[AC]",
                              "N": "[ACTG]",
                              "R": "[AG]",
                              "S": "[CG]",
                              "V": "[ACG]",
                              "W": "[AT]",
                              "Y": "[CT]",
                              "A": "A",
                              "C": "C",
                              "G": "G",
                              "T": "T"
                              }
        self.regex = "".join([self.sequence_code[x] for x in self.site])

    def digest(self, scaffold_name, scaffold_seq):
        '''
        Performs the insilico genome digestion of a given string of DNA by finding the recognition site and offsetting
        by the distance between start of recognition site with cut position
        '''
        return [x.start()+self.offset for x in re.finditer(self.regex, scaffold_seq)]

    def __str__(self):
        '''
        String representation of the enzyme
        '''
        return 'Enzyme: {}\tCut site: {}'.format(self.name, self.site)

    def __repr__(self):
        return self.__str__()


def digest_scaffold(output_file_object, scaffold_name, scaffold_seq, enzymes, retain_between=None):
    '''
    This function digests any given scaffold with any number of enzyme objects, it is made to be used inside the main function of the project
    but can be run independently, for that it requires:
    output_file_object: a file object to the output file (e.g. in ```f= open('output.txt','w')``` f is a file object)
    scaffold_name: name of the sequence to be digested
    scaffold_seq: Sequence to be digested
    enzymes: A list of Enzyme objects to be used to digest the genome
    retain_between: A list with two elements, start and stop of range, this will ensure that only fragments that fit inside the range are recorded 
    '''
    enzymes_cuts = [x.digest(scaffold_name, scaffold_seq) for x in enzymes]
    enzymes_dict = {}
    for index, cuts in enumerate(enzymes_cuts):
        for pos in cuts:
            enzymes_dict[pos] = enzymes[index].name
    cuts_sorted = sorted(list(enzymes_dict.keys()))
    # the first fragment will start on the beggining of the scaffold
    start = [0, "scaffold_start"]
    retained_length = 0
    for i in cuts_sorted:

        if retain_between is None:
            retained_length += i-start[0]
            output_file_object.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                scaffold_name, start[0], i, start[1], enzymes_dict[i], i-start[0]))
        else:
            if retain_between[0] <= i-start[0] <= retain_between[1]:
                retained_length += i-start[0]
                output_file_object.write(("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    scaffold_name, start[0]+1, i, start[1], enzymes_dict[i], i-start[0])))
        start = [i, enzymes_dict[i]]
    print(f"Finished digesting {scaffold_name} - retained {retained_length} bps of {len(scaffold_seq)} bps ({(retained_length*100)/len(scaffold_seq)} %)")


def read_fasta_and_digest(fasta_object, output_file_object, enzymes, retain_between=None):
    ''' 
    Function responsible for reading in the fasta file object and passing the individual scaffolds to the digest_scaffold function.
    Takes as arguments:
    fasta_object: A file object to the fasta file.
    output_file_object: A file object to the output file to be passed to digest_scaffold function
    enzymes: A list of Enzyme objects to be used to digest the genome, to be passed to digest_scaffold
    retain_between: A list with two elements, start and stop of range, this will ensure that only fragments that fit inside the range are recorded, to be passed to digest_scaffold 
    '''
    name = None
    seq = ""
    for line in fasta_object:
        data = line.strip()
        if data.startswith('>') and name is None:
            name = data.replace('>', '')
        elif data.startswith('>') and name is not None:
            digest_scaffold(output_file_object, name, seq.upper(), enzymes, retain_between)
            
            name = data.replace('>', '')
            seq = ""

        else:
            seq += data
    
    digest_scaffold(output_file_object, name, seq.upper(), enzymes, retain_between)
            

def main(sizes_to_retain,enzyme_list,genome_path,output_path):
    '''
    main function of insilico_digestion module. The function is responsible for:
    1) Build the Enzyme objects and keep them in a list
    2) Create the fasta file object from either a text fasta file or a gziped fasta  file.
    3) Create the output file object to write the results to
    4) Pass the fasta file to the read_fasta_and_digest function
    '''
    enzymes_object = []
    
    for s_enzyme in enzyme_list:
        enzyme_name, enzyme_cut = s_enzyme.strip().split(":")
        enzyme_cut = enzyme_cut.upper()
        enzymes_object.append(Enzyme(enzyme_name, enzyme_cut))

    if genome_path.endswith('gz') or genome_path.endswith('gzip'):
        with gzip.open(genome_path, 'rt') as fasta, \
                open(output_path, 'w') as output:
            read_fasta_and_digest(fasta, output, enzymes_object,
                    retain_between=sizes_to_retain)
    else:
        with open(genome_path) as fasta, \
                open(output_path, 'w') as output:
            read_fasta_and_digest(fasta, output, enzymes_object,
                    retain_between=sizes_to_retain)


if __name__ == '__main__':
    '''
    usage: insilico_digestion.py [-h] -e enzymes -f fasta -o output [-r range_sizes]

    Insilico digestion

    Required arguments:
    -e enzymes      A comma separated list of enzymes, e.g.: MspI:C^CGG,ApeKI:G^CWGC
    -f fasta        Path to the genome fasta file or fasta.gz
    -o output       Path to the output file

    Optional arguments:
    -r range_sizes  Range to retain fragments, default is none, range should be in the format minimum:maximum, e.g.: 60:300. Default = None
    '''

    parser = argparse.ArgumentParser(description="Insilico digestion")
    parser._action_groups = []  # remove all argument groups from the list
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-e', metavar='enzymes', type=str, required=True,
                          help='A comma separated list of enzymes, e.g.: MspI:C^CGG,ApeKI:G^CWGC')
    required.add_argument('-f', metavar='fasta',  type=str, required=True,
                          help="Path to the genome fasta file or fasta.gz")
    required.add_argument('-o', metavar='output', type=str, required=True,
                          help="Path to the output file")
    optional.add_argument('-r', metavar='range_sizes', type=str, default=None,
                          help="Range to retain fragments, default is none, range should be in the format\
                            minimum:maximum, e.g.: 60:300. Default = None")

    args = parser.parse_args()

    # This assertion should never fail since enzymes is a required argument
    enzyme_valid_symbols = {"B", "D", "H", "K", "M", "N",
                            "R", "S", "V", "W", "Y", "A", "C", "G", "T", "^"}
    assert len(args.e.split(",")) >= 1, "Something went wrong on the inputs"

    # Make sure that all entries on the enzymes have a column separating potential name and cutsite
    assert sum([1 if ":" in x else 0 for x in args.e.split(",")]) == len(
        args.e.split(",")), "Enzymes are not in the format of name:cutsite"

    # Check if any invalid characters are in the enzyme cutsites
    assert sum([1 if set(x.split(":")[1].upper()).issubset(enzyme_valid_symbols)
                else 0 for x in args.e.split(",")]) == len(args.e.split(',')), \
        "Invalid character in one of the cutsites, valid characters are {}".format(
            enzyme_valid_symbols)

    # Check if all enzymes have an actual cut position
    assert sum([1 if "^" in x else 0 for x in args.e.split(",")]) == len(
        args.e.split(",")), "One of the enzymes doesn't have a cutsite specified, cutsite is defined by a carrot symbol (^)"


    sizes_to_retain = list(map(int, args.r.split(
        ':'))) if args.r != None else None
    enzyme_list = args.e.split(",")
    genome_path = args.f
    output_path = args.o
    main(sizes_to_retain=sizes_to_retain,enzyme_list=enzyme_list,genome_path=genome_path,output_path=output_path)
