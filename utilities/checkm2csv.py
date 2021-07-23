import panda as pd
from pathlib import Path




cwd = str(Path.cwd())

def check2csv(file_path:str, csv_path):
    with open(file_path, 'r') as f:
        load = {}
        for line in f:
        #print(line)
            newline = line.split('\t')
            bin_id = newline[0]
            #print(bin_id)
            new = newline[1].strip().replace('\'','\"')
            #print(new)
            dat_dict = json.loads(new)
            load[bin_id] = dat_dict
    csv_path = cwd + '/' + csv_path
    with open(csv_path, 'w+') as output:
        output.write('Bin_id\tContig length\tCompleteness\tContamination\n')
        for key in load:
            output.write(key + '\t')
            output.write(str(load[key]['Genome size']) + '\t')
            output.write(str(load[key]['Completeness']) + '\t')
            output.write(str(load[key]['Contamination']) + '\t')
            output.write('\n')
    return pd.read_csv(csv_path, sep='\t', index_col=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help = 'input file')
    parser.add_argument('-o',required=True, help = 'output file' )