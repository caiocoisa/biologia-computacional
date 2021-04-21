import os
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description = 'test')

    parser.add_argument('-d', '--dataset', help = 'Dataset file' )
    parser.add_argument('-i', '--id_problem', help = 'ID of Rosalind problem' )

    return parser.parse_args()

def switch(id_problem, dataset):
    switcher = {
        'DNA': dna(dataset),
        'RNA': rna("caminho"),
    }
    return switcher.get(id_problem)

def dna(dataset):
    dna = open (dataset, 'r')
    a, c, g, t = 0, 0, 0, 0
    line = dna.read()

    for l in line:
        if l == 'A':
            a = a+1
        if l == 'C':
            c = c+1
        if l == 'G':
            g = g+1
        if l == 'T':
            t = t+1

    print("{} {} {} {}".format(a, c, g, t))

def replace_t(code):
        return re.sub('T', 'U', code)

def rna(dataset):
    dna = open(dataset, 'r')
    line = dna.read()
    result = map(replace_t, line)
    print("".join(list(result)))

def complement_symbol(symbol):
    if symbol == 'A':
        return 'T'
    elif symbol == 'T':
        return 'A'
    elif symbol == 'C':
        return 'G'
    elif symbol == 'G':
        return 'C'

def complement_dna(dna):
    return map(complement_symbol, dna)

def revc(dataset):
    dna = open(dataset, 'r')
    line = dna.read()
    lista = list(line)
    lista.reverse()
    line = "".join(lista)
    complement_line = complement_dna(line)
    x = list(complement_line) 
    s = "".join(x[1:])
    print(s)

if __name__ == '__main__':
    args = parse_args()
    id_lower = args.id_problem.lower()

    if id_lower == 'dna':
        dna(args.dataset)
    elif id_lower == 'rna':
        rna(args.dataset)
    elif id_lower == 'revc':
        revc(args.dataset)