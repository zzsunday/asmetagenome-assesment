from utilities.read_barrnap import read_barrnap
import pandas as pd
import logging
from pathlib import Path
import argparse
import os

FORMAT = "%(asctime)s\t%(message)s"
logging.basicConfig(level = logging.INFO, format=FORMAT, datefmt="%Y%m%d %H:%M:%S")

def extract_gene(gff_path,fasta_path,copy_number):
    """ extract single copy or multiple copy of 16s gene from 
    barrnap 
    gff_path : gff format of barrnap result
    fasta_path : fasta format of barrnap result
    
    """

    asm_id = []
    barrnap_dat = read_barrnap(gff_path)
    logging.info('the working dir is {}'.format(Path.cwd()))
    cwd = str(Path.cwd())
    if copy_number ==1:
        for x in list(barrnap_dat.groupby('bin_id')): 
            count = 0
            for i in x[1]['Ribosome']: # for each contig to find 16s gene
                if 'Name=16S_rRNA' in i:
                    asm_id.append(x[0])
        logging.info('the number of contig contain 16s gene is {}'.format(len(asm_id)))
        os.system('mkdir -p {}'.format(cwd + '/ctg_16_single'))
        for ctg in asm_id:
            f = open(fasta_path)
            for line in f: # extract contig ,which contains 16s gene
                if line.startswith('>16S_rRNA') and ctg in line:
                    output = cwd + '/ctg_16_single/'+ ctg + '_16s.fa'
                    out = open(output, 'a')
                    out.write(line)
                    seq = f.readline()
                    out.write(seq)
                    out.close()
                else:
                    next(f)
            f.close()
    else:
        for x in list(barrnap_dat.groupby('bin_id')): 
            count = 0
            for i in x[1]['Ribosome']: # for each contig to find 16s gene
                if 'Name=16S_rRNA' in i:
                        count+=1
            if count >=2:
                asm_id.append(x[0])
        logging.info('the number of contig contain more than one copy of 16s geneis {}'.format(len(asm_id)))
        os.system('mkdir -p {}'.format(cwd + '/ctg_16_multiple'))

        for ctg in asm_id:
            f = open(fasta_path)
            for line in f: # extract contig , whose copies of 16s gene are more than once 
                if line.startswith('>16S_rRNA') and ctg in line:
                    output = cwd + '/ctg_16_multiple/'+ ctg + '_16s.fa'
                    out = open(output, 'a')
                    out.write(line)
                    seq = f.readline()
                    out.write(seq)
                    out.close()
                else:
                    next(f)
            f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, help = 'barrnap gff file')
    parser.add_argument("-a", required=True,  help = "output")
    parser.add_argument("-n", type=int, nargs="?", default=2)
    arg = parser.parse_args()
    extract_gene(arg.i, arg.a,arg.n)