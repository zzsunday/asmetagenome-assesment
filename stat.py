import sys
import glob
import os
import argparse
from utilities.read_barrnap import read_barrnap
from utilities.checkm2csv import checkm2csv
from utilities.extract_genes import extract_gene
from utilities.read_prodigal import read_prodigal
from pathlib import Path
import logging

FORMAT = "%(asctime)s\t%(message)s"
logging.basicConfig(level = logging.INFO, format=FORMAT, datefmt="[%Y-%m-%d %H:%M:%S]")

#files = glob.glob('/data/sunda/04.Project/meta_hifi_hic/result/05.barrnap/canu/ctg_16/tmp/*.fasta')
#print(files)exit

description = "analyse the result of checkm barrnap prodigal"

cwd = str(Path.cwd())
#print(cwd)
stat_dir = cwd + '/asm_stats'
if Path(stat_dir).exists():
    pass
else:
    os.system("mkdir -p {}".format(stat_dir)) #create dir

out = open('{}/stat.csv'.format(stat_dir),'a') #open statistic file



def Get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('--checkm_dir', required=True, help='checkm result dir')
    parser.add_argument('--barrnap_dir', required=True, help = 'barrnap resultls dir')
    parser.add_argument('--vesearch', required=True, help = 'vesearch path')
    parser.add_argument('--prodigal_dir', required=True, help = 'prodigal dir')
    parser.add_argument('--cdhit_dir', required=True, help = 'cdhit dir')

    return parser.parse_args()


def Parse_checkm(input):
    """input : bin_stats_ext.tsv 
    """
    input_ = input + '/storage/bin_stats_ext.tsv'
    output_ = stat_dir + '/checkm.csv'
    dat = checkm2csv(input_,output_)
    length_10kbp = int(dat[dat['Contig length'] >10000].sum()[1]/1000000)
    length_100kbp = int(dat[dat['Contig length'] >1000000].sum()[1]/1000000)
    length_1Mbp =  int(dat[dat['Contig length'] >100000].sum()[1]/1000000)
    complete_90 = len(dat[(dat['Completeness']> 90 ) & (dat['Contamination'] < 5)])
    complete_25 = len(dat[dat['Completeness']> 25])
    out.write('Length in contigs>10kbp' + '\t' + str(length_10kbp) + 'Mbp' + '\n')
    out.write('Length in contigs>100kbp' + '\t' + str(length_100kbp) + 'Mbp' + '\n')
    out.write('Length in contigs>1Mbp' + '\t' + str(length_1Mbp) + 'Mbp' + '\n')
    out.write('Checkm>90/% complete' + '\t' + str(complete_90) +'\n')
    out.write('Checkm>25/% complete' + '\t' + str(complete_25) +'\n')






def Parse_barrnap(barrnap_dir, vesearch_path='vesearch',threads=16):
    #print("{}*barrnap.gff".format(barrnap_dir))
    file_gff = glob.glob("{}/*barrnap.gff".format(barrnap_dir))[0]
    logging.info('gff path is {}'.format(file_gff)) #
    file_fa =  glob.glob("{}/*barrnap.fa".format(barrnap_dir))[0] # get barrnap fa path
    logging.info('fasta path is {}'.format(file_fa))
    extract_gene(file_gff,file_fa) # extract contig, which contain more than one copy of 16s RNA
    extract_gene(file_gff,file_fa,copy_number=1) # extract contig contain genes
    #处理提取的多拷贝16s RNA
    multi_dir = stat_dir + '/ctg_16s_multiple' # multiple copy dir
    single_dir = stat_dir + '/ctg_16s_single'  #single copy dir
    try: 
        if Path(multi_dir).is_dir():
            pass    # cluter multiple_dir
        elif Path(single_dir).is_dir():
            pass
        else: 
            raise FileExistsError
    except FileExistsError:
        logging.info("\033[1;31merror: \033[0m 16s dir does not exists, plz check barrnap result")
        exit()

    
    #16s RNA clusters(95%)
    os.system('cat {}/*.fa > {}/16s_RNA.fa'.format(single_dir,single_dir))
    total_genes = cmd_count="grep '>16S' {}/16s_RNA.fa |wc -l".format(single_dir)
    total_genes_out = int(os.popen(cmd_count).read().strip())
    logging.info('the number of 16S RNA genes are {}'.format(total_genes_out))
    out.write("16s rRNA genes" + '\t' + str(total_genes_out) + '\n')
    
    # cluster at 95% similarity
    vesearch_cmd = '{} --cluster_fast {}/16s_RNA.fa --id 0.95 --threads {} --sizeout --centroids {}/16s_RNA_0.95.fa'.format(vesearch_path,single_dir,threads,single_dir)
    logging.info('CMD of cluster at 95\% similarity: {}'.format(vesearch_cmd))
    os.system(vesearch_cmd)
    cmd_count="grep '>16S' {}/16s_RNA_0.95.fa | wc -l".format(single_dir)
    #print(cmd_count)
    cmd_count_single = "grep -w 'size=1' {}/16s_RNA_0.95.fa |wc -l".format(single_dir)
    #print(cmd_count_single)
    #print(cmd_count)
    cmd_res = os.popen(cmd_count)
    cmd_res_single = os.popen(cmd_count_single)
    #print(cmd_count_single)
    a = int(cmd_res.read().strip())
    b = int(cmd_res_single.read().strip())
    logging.info('16s RNA clusters are {}'.format(a-b))
    out.write('16s rRNA clusters(95%)' + '\t' + str(a-b) +'\n')



    # # multiple copies of 16s genes
    # contigs with matching 16s (97% similarity)
    files = glob.glob("{}/*fa".format(multi_dir))
    chirimic_count = []
    tmp_dir = multi_dir + '/tmp'
    print(tmp_dir)
    os.system("mkdir -p {}".format(tmp_dir))
    for file in files:
        prefix = file.split('/')[-1].split('.')[0]
        vesearch_cmd = "{} --cluster_fast {} --id 0.97 --threads {} --centroids {}/{}_0.97.fa".format(vesearch_path,file,threads,tmp_dir,prefix)
        os.system(vesearch_cmd) 
        cmd_count = "grep '>16S' {}/{}_0.97.fa | wc -l".format(tmp_dir,prefix)
        cmd_res = os.popen(cmd_count)
        a = int(cmd_res.read().strip())
        chirimic_count.append(a)
    logging.info('Conitgs with matching 16s are {}'.format(chirimic_count.count(1)))
    out.write('Contigs with matching 16s' + '\t' + str(chirimic_count.count(1)) + '/' + str(len(files)) + '\n')
        

def Parse_cdhit(prodigal_dir, cdhit_dir):
    prodigal_faa = glob.glob("{}/*faa".format(prodigal_dir))[0]
    cdhit_faa = glob.glob("{}/*faa".format(cdhit_dir))[0]
    ORF = read_prodigal(prodigal_faa)
    ORF_99 = read_prodigal(cdhit_faa)
    out.write('Full length ORFS' + '\t' + str(ORF) + '\n')
    out.write('ORFS clusters(99%)' + '\t' + str(ORF_99) + '\n')

    









if __name__ == "__main__":
    arg = Get_arg()
    if arg.checkm_dir.endswith('/'):
        arg.checkm_dir = arg.checkm_dir.strip('/')
    Parse_checkm(arg.checkm_dir)
    logging.info('checkm analysis finished')
    if arg.barrnap_dir.endswith('/'):
        arg.barrnap_dir = arg.barrnap_dir.strip('/')  
    Parse_barrnap(arg.barrnap_dir,arg.vesearch)
    logging.info('barrnap analysis finished')
    if arg.cdhit_dir.endswith('/'):
        arg.cdhit_dir = arg.cdhit_dir.strip('/')
    Parse_cdhit(arg.prodigal_dir,arg.cdhit_dir)
    logging.info('analysis finished')
    out.close()