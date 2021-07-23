import argparse
import os
from pathlib import Path
import logging


description = """description: assement.py combine checkm barnnap prodgial to assess metagenomics assembly"""
usage = """ python assement.py -i input -o output -b bins
        eg:python assesment.py -i test_data/test.fa -o test -b test_data/bins/ 
===============================================================================
        """

FORMAT = "%(asctime)s\t%(message)s"
logging.basicConfig(level = logging.INFO, format=FORMAT, datefmt="[%Y-%m-%d %H:%M:%S]")

cwd = str(Path.cwd())

def Parsearg():

    """ get args form argparse """
    parser = argparse.ArgumentParser(description=description, usage=usage,
                            formatter_class=argparse.RawTextHelpFormatter        ) # -t: the number of threads 
    parser.add_argument('-i', '--input', required=True, help='input file')
    parser.add_argument('-o','--output' , required=True, help = 'output')
    parser.add_argument("--checkm", nargs='?',default='checkm', help = 'checkm path')
    parser.add_argument("--barrnap", nargs='?',default='barrnap', help = 'barrnap path')
    parser.add_argument('--prodigal', nargs="?", default='prodigal', help = 'prodigal path')
    parser.add_argument("-b", required=True, help = "the split dir of assemble fasta dir")
    return parser.parse_args()

def run_checkm(input, path='checkm',threads=8, prefix='fa'): # need to change 
    logging.info('The working dir is {}'.format(cwd))
    checkm_dir = cwd + '/' + 'checkm'
    print(checkm_dir)
    try: 
        if not Path(checkm_dir).is_dir():
            os.system('mkdir {}'.format(checkm_dir))
        else: 
            raise FileExistsError
    except FileExistsError:
        logging.info("\033[1;31merror: \033[0m checkm dir exists, plz remove")
        exit()

    cmd = '{} lineage_wf -t {} -x {} --nt --tab_table -f bins_stats.txt {} {}'.format(path, threads, prefix, input, checkm_dir)
    logging.info('cmd: {}'.format(cmd))
    os.system(cmd)
    stat_dir = cwd + '/' + 'bins_stats.txt'
    cmd1 = 'mv {} {}'.format(stat_dir, checkm_dir)
    os.system(cmd1)

def run_barrnap(input, output, path='barrnap', threads=16):
    barrnap_dir = cwd + '/' + 'barrnap'
    output = barrnap_dir + '/' + output

    try: 
        if not Path(barrnap_dir).is_dir():
            os.system('mkdir {}'.format(barrnap_dir))
        else: 
            raise FileExistsError
    except FileExistsError:
        logging.info("\033[1;31merror: \033[0m barrnap dir exists, plz remove")
        exit()
    cmd = '{} --thread {} -o {}_barrnap.fa < {} 1> {}_barrnap.gff 2>{}_barrnap_err.log'.format(path,threads,output,input,output,output)
    logging.info('cmd: {}'.format(cmd))
    os.system(cmd)

def run_prodigal(input,output, path='prodigal'):
    prodigal_dir = cwd + '/' + 'prodigal'
    output = prodigal_dir + '/' + output
    try: 
        if not Path(prodigal_dir).is_dir():
            os.system('mkdir {}'.format(prodigal_dir))
        else: 
            raise FileExistsError
    except FileExistsError:
        logging.info("\033[1;31merror: \033[0m prodigal dir exists, plz remove")
        exit()
    cmd = '{} -i {} -o {}.gbk -a {}.faa -a {},ffn -p meta'.format(path,input,output,output,output)
    logging.info('cmd: {}'.format(cmd))
    os.system(cmd)
if __name__ == "__main__":
    arg = Parsearg()
    run_checkm(arg.b)
    logging.info('checkm is finished')
    logging.info('barrnap starts')
    run_barrnap(arg.i,arg.o)
    logging.info('barrnap is finished')
    logging.info('prodigal starts')
    run_prodigal(arg.i,arg.o)
    logging.info('Job is finished')
