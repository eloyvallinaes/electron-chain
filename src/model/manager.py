import os
import glob
import subprocess as sp

class FastaManager():
    '''
    Organize a [new] folder in electron-chain/models/ in order to run MSAs
    and AF2 over a set of selected protein sequences.
    
    Initialize with absolute path of fasta file. Filename, trimmed of
    format, will be the new folder name
    '''

    def __init__(self, fasta_path, models_folder):
        self.fasta_path = fasta_path
        self.models_folder = models_folder
 
        # check cdhit exe

    def run_cdhit_clustering(thr=0.9):
        command = 'cd-hit'
        sp.run()

    def split_fasta():



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
