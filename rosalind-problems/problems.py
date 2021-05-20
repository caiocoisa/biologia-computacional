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

def fibo_rabbit(n, k):
    month1 = 1
    month2 = 1
    for i in range(2, n):
        population = month1 + (month2 * k)
        month2 = month1
        month1 = population
    return population

def fib(dataset):
    dna = open(dataset, 'r')
    line = dna.read()
    lista = line.split(' ')
    a = fibo_rabbit(int(lista[0]), int(lista[1]))
    print(a)

#
def gc(dataset):
    fasta = open(dataset, 'r')
    fasta_stripped = (line.strip() for line in fasta)
    fasta_splitted = (line.split(">") for line in fasta_stripped if line)

    dic = {}
    id_dna = ''
    dna = ''
    for info in fasta_splitted:
        
        if info[0] == '':
            if id_dna != '' : dic[id_dna] = dna
            id_dna = info[1]
            dna = ''
        
        else:
            dna = ''.join([dna, info[0]])

    output = {}
    for key, value in dic.items():
        
        qtd=0
        for char in value:
            if char in ['C', 'G']:
                qtd += 1
        perc = 100 * ( qtd/len(value))
        output[key] = perc
    print(output)

if __name__ == '__main__':
    args = parse_args()
    id_lower = args.id_problem.lower()

    if id_lower == 'dna':
        dna(args.dataset)
    elif id_lower == 'rna':
        rna(args.dataset)
    elif id_lower == 'revc':
        revc(args.dataset)
    elif id_lower == 'fib':
        fib(args.dataset)
    elif id_lower == 'gc':
        gc(args.dataset)
        