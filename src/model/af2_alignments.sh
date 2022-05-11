#!/bin/bash -x
#SBATCH -A SNIC2021-5-297
#SBATCH --output=/proj/nobackup/snic2019-35-62/gpozzati/job_out/hh%j.out
#SBATCH --error=/proj/nobackup/snic2019-35-62/gpozzati/job_err/hh%j.err
#SBATCH --array=1-1000
#SBATCH -c 8
#SBATCH -t 04:00:00

LIST=$1
OFFSET=$2
FOLDER=$3

##### CONFIG #####
PROJ="/proj/nobackup/snic2019-35-62"	# common path to everything is used in this script
HOME="$PROJ/gpozzati/electron-chain/models/"	# folder structure with data and scripts

DB="$PROJ/Database"	# common folder containing all databases
UNI="$DB/uniref90/uniref90.fasta"	# uniref90 path
MGN="$DB/mgnify/mgy_clusters.fa"	# mgnify path
BFD="$DB/small_bfd/bfd-first_non_consensus_sequences.fasta"	# small bfd path
PDB70="$DB/pdb70_220313/pdb70"	# pdb70 path
UNICL30="$DB/uniclust30_2018_08/uniclust30_2018_08"  # uniclust30 path
SEQRES="$DB/pdb_seqres/pdb_seqres.txt"
MMCIF="$DB/pdb_mmcif/pdb_mmcif/"
OBS="$DB/pdb_mmcif/obsolete.dat"

AF2="$PROJ/gpozzati/af2-v2.2.0"	# local installation of AF2
AFDATA="$AF2/AF_data_v220"	# data to run AF2
SING="$AFDATA/alphafold_v220.sif"	# singularity image
##################

pos=$(($SLURM_ARRAY_TASK_ID + $OFFSET))
ID=`tail -n+$pos $LIST | head -n 1`

echo 'Processing ' $ID '...'

ml singularity/3.8.2

singularity exec --bind $PROJ:$PROJ $SING \
	python3 $AF2/alphafold/run_alphafold.py \
		--fasta_paths=$HOME/$FOLDER/models/${ID}.fasta \
		--model_preset=monomer \
		--output_dir=$HOME/$FOLDER/models/ \
		--data_dir=$AFDATA \
		--run_msas_only=True \
		--use_templates=True \
		--db_preset=reduced_dbs \
		--uniref90_database_path=$UNI \
		--mgnify_database_path=$MGN \
		--small_bfd_database_path=$BFD \
		--pdb70_database_path=$PDB70 \
                --template_mmcif_dir=$MMCIF \
                --obsolete_pdbs_path=$OBS \
		--max_template_date=9999-01-01
