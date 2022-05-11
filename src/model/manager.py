import os
import glob
import shutil
import subprocess as sp

class FastaManager():
    '''
    Organize a [new] folder in electron-chain/models/ in order to run MSAs
    and AF2 over a set of selected protein sequences.
    
    Initialize with absolute path of fasta file. Filename, trimmed of
    format, will be the new folder name.

    '''

    def __init__(self, fasta_path, models_folder):
        self.fasta_path = fasta_path
        
        self.models_folder = models_folder
        if not os.path.exists(self.models_folder):
            os.mkdir(self.models_folder)
        
        newset = fasta_path.split('/')[-1].split('.')[0]
        self.models_folder = models_folder+'/'+newset
        if not os.path.exists(self.models_folder):
            os.mkdir(self.models_folder)

        self.cluster_path = self.models_folder+'/clusters.tsv'

        self.single_fasta_folder = self.models_folder+'/models/'
        if not os.path.exists(self.single_fasta_folder):
            os.mkdir(self.single_fasta_folder)

    def run_clustering(self, binary, thr='0.9'):
        command = [binary,
                   'createdb',
                   self.fasta_path,
                   self.models_folder+'/seqdb']
        sp.run(command)
        
        if not os.path.exists(self.models_folder+'/tmp'):
            os.mkdir(self.models_folder+'/tmp')
        
        command = [binary,
                   'cluster',
                   self.models_folder+'/seqdb',
                   self.models_folder+'/clustering',
                   self.models_folder+'/tmp',
                   '--cluster-reassign',
                   '1',
                   '--min-seq-id',
                   thr]
        sp.run(command)

        command = [binary,
                   'createtsv',
                   self.models_folder+'/seqdb',
                   self.models_folder+'/seqdb',
                   self.models_folder+'/clustering',
                   self.cluster_path]
        sp.run(command)

        for path in glob.glob(self.models_folder+'/clustering.*'):
            os.remove(path)
        shutil.rmtree(self.models_folder+'/tmp')

    def pick_fasta(self):
        with open(self.fasta_path) as fasta_file:
            fasta = ''.join([line for line in fasta_file])

        fasta = fasta.split('>')[1:]
        fasta = {'>'+entry.split('\n')[0].split()[0]:entry.split('\n')[1] \
                 for entry in fasta}
 
        representatives = []
        with open(self.cluster_path) as clusters:
            for line in clusters:
                repres = line.rstrip().split()[0]
                if repres not in representatives:
                    representatives.append(repres)

        with open(self.models_folder+'/id_list', 'w') as outlist:
            for code, seq in fasta.items():
                filename = code.split('|')[0].strip('>')
                single_fasta_path = self.single_fasta_folder+filename+'.fasta'
                if code.strip('>') in representatives:
                    outlist.write(filename+'\n')
                    with open(single_fasta_path, 'w') as out:
                        out.write(code+'\n'+seq+'\n')


class ModelManager():
    
    '''
    Manages AF2 output models in a folder contained in electron-chain/models/
    whose substructure is organized in the following way:

    folder/
        af_out/
            target1/
            target2/
            ...
    
    Initialize with the absolute folder path.

    '''
    
    def __init__(self, folder):
        self.folder = folder
        self.af_targets = glob.glob(folder+'/af_out/*')

    def rename(self):
        renamed_folder = self.folder+'/models'
        if not os.path.exists(renamed_folder):
            os.mkdir(renamed_folder)

        for target_folder in self.af_targets:
            target_name = target_folder.split('/')[-1]
            if target_name == '': continue

            for model in glob.glob(target_folder+'/ranked_*.pdb'):
                model_suffix = model.split('_')[-1]
                os.rename(
                    model,
                    f'{renamed_folder}/{target_name}_{model_suffix}')

    def trim(self, thr=75):
        trimmed_folder = self.folder+'/trimmed_models'
        if not os.path.exists(trimmed_folder): 
            os.mkdir(trimmed_folder)

        for model in glob.glob(self.folder+'/models/*'):
            modelname = model.split('/')[-1]
            with open(f'{trimmed_folder}/{modelname}', 'w') as out:
                
                for line in open(model):
                    if line.startswith('ATOM'):
                        if float(line[60:66]) > 70: out.write(line)
                    else: out.write(line)

    def run_foldseek_clustering(self):
        print ('To implement!')

    def parse_clustering(self):
        print ('To implement!')

    def plot_clustering(self):
        print ('To implement!')
