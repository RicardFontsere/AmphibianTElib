#!/bin/bash

#SBATCH --output=out/%x.%j.out
#SBATCH --error=out/%x.%j.err

# VUB-HPC specific customizations - do not modify!

OS_SUBMIT=$VSC_OS_LOCAL
ARCH_SUBMIT=$VSC_ARCH_LOCAL
SUFFIX_SUBMIT=$VSC_ARCH_SUFFIX


# workaround to make --export=ALL work in Hydra on different partitions
unset VSCPROFILELOADED VSC_OS_LOCAL VSC_ARCH_LOCAL VSC_ARCH_SUFFIX
source /etc/profile.d/vsc.sh
module update

# workaround to make sure we use the correct python executable in compute nodes
envpython=${EBROOTPYTHON/$OS_SUBMIT\/$ARCH_SUBMIT$SUFFIX_SUBMIT/$VSC_OS_LOCAL\/$VSC_ARCH_LOCAL$VSC_ARCH_SUFFIX}/bin/python
eval "$envpython() { python \"\$@\"; }"


#########################################################################################

# properties = {"type": "single", "rule": "edta", "local": false, "input": ["/user/brussel/109/vsc10945/home/scratch/AmphibianTELibrary/Genomes/Eleutherodactylus_coqui/Eleutherodactylus_coqui_headers.fna"], "output": ["/user/brussel/109/vsc10945/home/scratch/AmphibianTELibrary/Genomes/Eleutherodactylus_coqui/Eleutherodactylus_coqui.denovo"], "wildcards": {"species": "Eleutherodactylus_coqui"}, "params": {"species_dir": "/user/brussel/109/vsc10945/home/scratch/AmphibianTELibrary/Genomes/Eleutherodactylus_coqui", "curatedlib": "/user/brussel/109/vsc10945/home/scratch/TE/Libraries/repbase/31/RepBase31.03.fasta/vrtrep.ref"}, "log": ["/user/brussel/109/vsc10945/home/scratch/AmphibianTELibrary/LOGS/edta/edta_Eleutherodactylus_coqui.log"], "threads": 1, "resources": {"mem_mb": 6511, "mem_mib": 6210, "disk_mb": 100000, "disk_mib": 95368, "tmpdir": "<TBD>", "nodes": 1, "tasks": 1, "cpus_per_task": 16, "mem_mb_per_cpu": 25000, "runtime": 7200, "partition": "zen4,zen5_mpi"}, "jobid": 85}
cd /rhea/scratch/brussel/vo/000/bvo00034/vsc10945/Snakemake/AmphibianTElib && /apps/brussel/RL9/zen2-ib/software/Python/3.12.3-GCCcore-13.3.0/bin/python -m snakemake --snakefile '/rhea/scratch/brussel/vo/000/bvo00034/vsc10945/Snakemake/AmphibianTElib/Snakefile' --target-jobs 'edta:species=Eleutherodactylus_coqui' --allowed-rules edta --cores 24 --attempt 1 --force-use-threads  --resources 'mem_mb=6511' 'mem_mib=6210' 'disk_mb=100000' 'disk_mib=95368' 'nodes=1' 'tasks=1' 'cpus_per_task=16' 'mem_mb_per_cpu=25000' --wait-for-files '/rhea/scratch/brussel/vo/000/bvo00034/vsc10945/Snakemake/AmphibianTElib/.snakemake/tmp.n1eow8ha' '/user/brussel/109/vsc10945/home/scratch/AmphibianTELibrary/Genomes/Eleutherodactylus_coqui/Eleutherodactylus_coqui_headers.fna' --force --target-files-omit-workdir-adjustment --keep-storage-local-copies --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers input code software-env params mtime --deployment-method env-modules --conda-frontend 'conda' --shared-fs-usage software-deployment source-cache sources persistence input-output storage-local-copies --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 60 --scheduler 'ilp' --local-storage-prefix .snakemake/storage --scheduler-solver-path '/apps/brussel/RL9/zen2-ib/software/Python/3.12.3-GCCcore-13.3.0/bin' --configfiles /rhea/scratch/brussel/vo/000/bvo00034/vsc10945/Snakemake/AmphibianTElib/master/config/config.yaml --default-resources base64//bWVtX21iPW1pbihtYXgoMippbnB1dC5zaXplX21iLCAxMDAwKSwgODAwMCk= base64//ZGlza19tYj0xMDAwMDA= base64//dG1wZGlyPXN5c3RlbV90bXBkaXI= base64//bm9kZXM9MQ== base64//dGFza3M9MQ== base64//Y3B1c19wZXJfdGFzaz02 base64//bWVtX21iX3Blcl9jcHU9NDAwMA== base64//cnVudGltZT00MzIw base64//cGFydGl0aW9uPXplbjQsemVuNV9tcGk= --mode 'remote' && touch '/rhea/scratch/brussel/vo/000/bvo00034/vsc10945/Snakemake/AmphibianTElib/.snakemake/tmp.n1eow8ha/85.jobfinished' || (touch '/rhea/scratch/brussel/vo/000/bvo00034/vsc10945/Snakemake/AmphibianTElib/.snakemake/tmp.n1eow8ha/85.jobfailed'; exit 1)



