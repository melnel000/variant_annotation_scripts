#!/bin/bash
cat << EOF
#!/bin/bash
#SBATCH -D ${1}
#SBATCH -J vep_${2}
#SBATCH -o ${1}/vep_${2}.log
#SBATCH --no-requeue
#SBATCH --mem 16G
#SBATCH --time 48:00:00
#SBATCH -p Main

shopt -s expand_aliases
module add vep/singularity

vep -i $3 -o $2_vep.vcf --cache --offline --fork 16 --buffer_size 1000000 \
--sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical \
--protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_esp --af_gnomad --max_af \
--pubmed --var_synonyms --variant_class --mane \
--flag_pick --pick_order mane,canonical,appris,tsl,biotype,ccds,rank,length \
--nearest symbol --overlaps $4 --force_overwrite --vcf \
--fasta /cbio/dbs/gatk/hg38/Homo_sapiens_assembly38.fasta \
--custom $5,gnomADv3genomes,vcf,exact,0,AF_amr,AF_nfe,AF_eas,AF_fin,AF_asj,AF_sas,AF_afr,AF_oth,AF_ami \
--plugin SpliceAI,snv=/cbio/projects/003/soft/SpliceAI/spliceai/annotations/spliceai_scores.raw.snv.hg38.vcf.gz,\
indel=/cbio/projects/003/soft/SpliceAI/spliceai/annotations/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
--plugin SpliceRegion \
--plugin GeneSplicer,/cbio/projects/003/soft/GeneSplicer/sources/genesplicer,/cbio/projects/003/soft/GeneSplicer/training_data_sets/Human \
--plugin dbscSNV,/cbio/dbs/snpeff/human/b38/dbscSNV1.1_GRCh38.txt.gz \
--plugin MaxEntScan,/software/data/bio/vep/Plugins/files/MaxEntScan/fordownload,SWA,NCSS
filter_vep -i $2_vep.vcf --format vcf -o $2_vep_filtered.txt -filter "PICK"

EOF
