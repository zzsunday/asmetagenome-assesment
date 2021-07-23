import sys
import glob
import os
##  未完成 进行结果统计
files = glob.glob('/data/sunda/04.Project/meta_hifi_hic/result/05.barrnap/canu/ctg_16/tmp/*.fasta')
#print(files)

out = open('stat.csv','a')
for file in files:
    cmd = "grep '>' {}| wc -l".format(file)
    p = os.popen(cmd) 
    out.write(file+'\t'+str(p.read().strip())+'\n')
    
out.close()