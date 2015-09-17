### Specify where the main pipeline file is (in general use full paths)
pipeline=/home/ucbtsmo/Projects/viral_denovo_pipeline/scripts/Morfopoulou_pipeline_viral_LEGION_v5.sh


### Specify your working directory
results=/home/ucbtsmo/Scratch


## Name of your sample
sample=10042

if [ ! -e ${results} ]; then mkdir ${results}; fi

outDir=${results}/${sample}                                                                                     ####output directory      


##location of input files
inputFiles="2 /home/ucbtsmo/Projects/viral_denovo_pipeline/example_quick/rawData/UCLGNS1108-EBV-10042_*R1*.fastq.gz /home/ucbtsmo/Projects/viral_denovo_pipeline/example_quick/rawData/UCLGNS1108-EBV-10042_*R2*.fastq.gz"



scriptStep1=${outDir}/cluster/submission/step.sh
 dbViral=/home/ucbtsmo/Projects/viral_denovo_pipeline/blast_databases/EBV_blastdb/EBV_124.fasta   
 dbViral_single=/home/ucbtsmo/Projects/viral_denovo_pipeline/blast_databases/EBV_single_db/EBV.fasta


genome=continuous #### set value to continuous for all pathogens that have non-segmented genome (i.e this flag should now be used for EBV as well) or  fluA (value needs to be FLUA)
instrument=NextSeq    ### can change to HiSeq or MiSeq 
nthreads=4               ### number of processors
collapseLev=72

bash ${pipeline} --script ${scriptStep1} --inputFiles ${inputFiles}  --outDir ${outDir} --sample ${sample} --db ${dbViral} --dbSingle ${dbViral_single}  --type ${instrument}  --genome ${genome} --nthreads ${nthreads} --collapseLev ${collapseLev}    ### collapsing, merging, trimming,blastn, assembly, pseudogenome, variants                                                                                                                                      



