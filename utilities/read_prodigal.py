

__description__ = """ to extract full length of ORF from prodigal """



def read_prodigal(file_path:str): # gbk format
    with open(file_path, 'r') as f:
        count = 0
        for line in f:
            if line.startswith(">") and "partial=00" in line:
                count +=1
    return count
 
