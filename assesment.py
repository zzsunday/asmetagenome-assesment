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

def Get_arg():

    """ get args form argparse """
    parser = argparse.ArgumentParser(description=description, usage=usage,
                            formatter_class=argparse.RawTextHelpFormatter        ) # -t: the number of threads 
    parser.add_argument("-i" ,required=True, help="input file")
    parser.add_argument("-o", nargs='?', default='.', help = "output")
    parser.add_argument("--checkm", nargs="?",default="checkm", help = "checkm path")
    parser.add_argument("--barrnap", nargs="?",default="barrnap", help = "barrnap path")
    parser.add_argument("--prodigal", nargs="?", default="prodigal", help = "prodigal path")
    parser.add_argument('--cdhit', nargs="?", default='cd-hit', help = 'cd-hit path')
    parser.add_argument("-b", required=True, help = "the split dir of assemble fasta dir")
    parser.add_argument("-c", nargs="?",default=0.99, help = " seqence identity")
    parser.add_argument('-t', nargs="?", default=16, help = 'number of threads')
    return parser.parse_args()

def run_checkm(input, path="checkm",threads=8, prefix="fa"): # need to change 
    logging.info("The working dir is {}".format(cwd))
    checkm_dir = cwd + "/" + "checkm"
    print(checkm_dir)
    try: 
        if not Path(checkm_dir).is_dir():
            os.system("mkdir {}".format(checkm_dir))
        else: 
            raise FileExistsError
    except FileExistsError:
        logging.info("\033[1;31merror: \033[0m checkm dir exists, plz remove")
        exit()

    cmd = "{} lineage_wf -t {} -x {} --nt --tab_table -f bins_stats.txt {} {}".format(path, threads, prefix, input, checkm_dir)
    logging.info("cmd: {}".format(cmd))
    os.system(cmd)
    stat_dir = cwd + "/" + "bins_stats.txt"
    cmd1 = "mv {} {}".format(stat_dir, checkm_dir)
    os.system(cmd1)

def run_barrnap(input, output, path="barrnap", threads=16):
    barrnap_dir = cwd + "/" + "barrnap"
    output = barrnap_dir + "/" + output

    try: 
        if not Path(barrnap_dir).is_dir():
            os.system("mkdir {}".format(barrnap_dir))
        else: 
            raise FileExistsError
    except FileExistsError:
        logging.info("\033[1;31merror: \033[0m barrnap dir exists, plz remove")
        exit()
    cmd = "{} --thread {} -o {}_barrnap.fa < {} 1> {}_barrnap.gff 2>{}_barrnap_err.log".format(path,threads,output,input,output,output)
    logging.info("cmd: {}".format(cmd))
    os.system(cmd)

def run_prodigal(input,output, path="prodigal"):
    prodigal_dir = cwd + "/prodigal"
    output = prodigal_dir + "/" + output
    try: 
        if not Path(prodigal_dir).is_dir():
            os.system("mkdir {}".format(prodigal_dir))
        else: 
            raise FileExistsError
    except FileExistsError:
        logging.info("\033[1;31merror: \033[0m prodigal dir exists, plz remove")
        exit()
    cmd = "{} -i {} -o {}.gbk -a {}.faa -p meta".format(path,input,output,output)
    logging.info("cmd: {}".format(cmd))
    os.system(cmd)

def run_cdhit(input,output,similarity,path='cd-hit'):
    cdhit_dir = cwd + "/cd-hit"
    output = cdhit_dir + "/" + output
    try: 
        if not Path(cdhit_dir).is_dir():
            os.system("mkdir {}".format(cdhit_dir))
        else: 
            raise FileExistsError
    except FileExistsError:
        logging.info("\033[1;31merror: \033[0m prodigal dir exists, plz remove")
        exit()
    cmd = "{} -i {} -o {}_{}.faa -aS 0.9 -c {} -G 0 -M 0 -T 9 -n 5 -d 0 -g 1".format(path,input,output,similarity,similarity)
    logging.info("cmd: {}",format(cmd))
    os.system(cmd)

if __name__ == "__main__":
    arg = Get_arg()
    run_checkm(arg.b)
    logging.info("checkm is finished")
    logging.info("barrnap starts")
    run_barrnap(arg.i,arg.o)
    logging.info("barrnap is finished")
    logging.info("prodigal starts")
    run_prodigal(arg.i,arg.o)
    logging.info("prodigal is finished")
    logging.info('cd-hit starts')
    run_cdhit(arg.i,arg.o,arg.c,arg.cdhit)
    print(arg.i,arg.o,arg.c,arg.cdhit)