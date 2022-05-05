import sys
import glob
import subprocess as sp

folder = sys.argv[1]

prev = ''
count = 0
with open(folder+'/clustering.tsv') as tsv:
    for line in tsv:
        repres = line.rstrip().split()[0]
        member = line.rstrip().split()[1]
        if repres != prev: count += 1
        rename = str(count)+'_'+member

        sp.run('mv {}/models_trimmed/{} {}/models_trimmed/{}'\
                .format(folder, member, folder, rename), shell=True)
        prev = repres



