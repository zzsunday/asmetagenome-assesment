import pandas as pd
import argparse 

__descirption__ = "read barrnap results to dataFrame"
__author__ = "sunda"
__date__  = "2021.7.19"


def read_barrnap(file_path): 
    """ file path : gff format of barrnap results"""
    dat = pd.read_csv(file_path, skiprows=1, sep='\t', header=None)
    barrnap_dat = dat[[0,2,3,4,6,8]]
    barrnap_dat.columns = ['bin_id', 'Gene type', 'start','end', 'strand', 'Ribosome']
    return barrnap_dat


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True,help = 'input file')
    parser.add_argument("-o", nargs="?", default='barrnap.csv', help = "output")
    arg = parser.parse_args()
    barrnap_dat = read_barrnap(arg.i)
    barrnap_dat.to_csv(arg.o, index=False)