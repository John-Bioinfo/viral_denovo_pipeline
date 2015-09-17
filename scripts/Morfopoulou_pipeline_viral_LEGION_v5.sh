#!/bin/bash                                                                                                                                                                                               
########### Here are the 2 lines that need to be modified                                                                                                                                                 
Illumeta=/home/ucbtsmo/Projects/viral_denovo_pipeline
R=/home/ucbtsmo/R/x86_64-unknown-linux-gnu-library/3.1/


########## everything else below should be automated                       
Software=${Illumeta}/exec
python=${Software}/python-2.7.6/bin/python
fastqCollapse2=${Software}/fastqCollapse2/fastqCollapse
blastn=${Software}/ncbi-blast-2.2.29+/bin/blastn
novoalign=${Software}/novoalign
novoindex=${Software}/novoindex
pear=${Software}/pear-0.9.0-src/src/pear
trimG=${Software}/trim_galore_0.3.7/trim_galore
seqtk=${Software}/seqtk-master/seqtk
spades=${Software}/SPAdes-3.5.0-Linux/bin/spades.py
quast=${Software}/quast-2.3/quast.py
samtools=${Software}/samtools-0.1.19/samtools
bcftools=${Software}/samtools-0.1.19/bcftools/bcftools
vcfutils=${Software}/samtools-0.1.19/bcftools/vcfutils.pl
varscan=${Software}/VarScan.v2.3.7.jar
extractLargeContigs=${Illumeta}/scripts/extract_largeContigs_new.pl
stichContigs=${Illumeta}/scripts/stichContigs.R  ##general                                                                                                                                                 
stichContigs_segm=${Illumeta}/scripts/stichContigs_flu.R ##for flu                                                                                                                                         
consensus=${Software}/QUASR_v7.03/pileupConsensus.jar
picard=${Software}/picard-tools-1.138/picard.jar


  
############ default values
name=contigs
cutoff=200


until [ -z "$1" ]; do
	# use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--inputFiles)
	    shift
	    i=0
	    for fileloc in $@; do 
		inputFiles[ $i ]=$fileloc
		((i=i+1))
	    done;;
	--script)
	    shift
	    script=$1;;
	 --db)
            shift
            db=$1;;
	 --dbSingle)
            shift
            dbSingle=$1;;
	 --genome)
            shift
            genome=$1;;
         --sample )
           shift
           sample=$1;;
         --type )
           shift
           type=$1;;
	--step1 )
	    shift
	    step1=$1;;
	--nthreads )
	    shift
	    nthreads=$1;;
	--outDir)
	    shift
	    output=$1;;
	--collapseLev)
	    shift
	    collapseLev=$1;;


	-* )
	    echo "Unrecognized option: $1"
	    exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done 


################ creating all the output folders
echo -e "Output folder: $output\n"


clusterDir=${output}/cluster
cluster_out=${clusterDir}/out
cluster_error=${clusterDir}/error
cluster_submission=${clusterDir}/submission
output_qc=${output}/QC
output_collapsed=${output}/collapsed
output_merged=${output}/merged
output_blastn=${output}/blastn
output_assembly=${output}/assembly
output_variants=${output}/variants

myFolders="$output $clusterDir $cluster_out $cluster_error $cluster_submission $output_assembly  $output_merged $output_qc $output_collapsed  $output_blastn $output_variants"


 for folder in $myFolders; do
    if [ ! -e $folder ]; then 
	echo "Creating $folder"
	mkdir $folder
    fi
done


 nfiles=${inputFiles[0]}

 if [[ "$nfiles" != "2" ]]; then
     echo "You specified $nfiles input files".
     echo "Error: currently the input data MUST be paired end."
     exit;
 fi


 seq1=${inputFiles[ 1 ]}
 seq2=${inputFiles[ 2 ]}


###############  check that raw data files  exist                                                                                                                                                            
 for file in $seq1 $seq2; do
     ls -lh $file
     if [ ! -e "$file" ]; then
         echo "Error, file $file does not exist"
         exit
     fi
 done


printf "#!/bin/bash -l \n#$ -S /bin/bash                                                                                                                                                                    
" >> $script

echo "Output script:  $script"




######################### Write the header of the script
###   Time/Memory required are estimates. If you know your datasets require more time/memory, modify as suitable.

######## MiSeq
if [[ "$type" == "MiSeq" ]]; then
     printf "#$ -l h_rt=02:30:00                                                                                                                                                                        
#$ -l thr=${nthreads}                                                                                                                               
#$ -l mem=1.9G                                                                                                                                                                                            
#$ -N viral_denovo_pipeline                                                                                                                                                                      
#$ -wd  ${output}                                                                                                                                                                               

cd \$TMPDIR                                                                                                                                                                                                 

module unload compilers/intel/11.1/072                                                                                                                                                                    
module unload mpi/qlogic/1.2.7/intel                                                                                                                                                                     
module unload mkl/10.2.5/035                                                                                                                                                                     
module load recommended/r                                                                                                                                                                                

export R_LIBS=${R}:$R_LIBS                                                                                                                                                                 
" >> $script

  	fi


########## NextSeq
if [[ "$type" == "NextSeq" ]]; then
     printf "#$ -l h_rt=08:30:00                                                                                                                                                                        
#$ -l thr=${nthreads}                                                                                                                               
#$ -l mem=1.9G                                                                                                                                                                                            
#$ -N viral_denovo_pipeline                                                                                                                                                                      
#$ -wd  ${output}                                                                                                                                                                               

cd \$TMPDIR                                                                                                                                                                                                 

module unload compilers/intel/11.1/072                                                                                                                                                                    
module unload mpi/qlogic/1.2.7/intel                                                                                                                                                                     
module unload mkl/10.2.5/035                                                                                                                                                                     
module load recommended/r                                                                                                                                                                                

export R_LIBS=${R}:$R_LIBS                                                                                                                                                                 
" >> $script

  	fi

############### HiSeq
if [[ "$type" == "HiSeq" ]]; then
     printf "#$ -l h_rt=12:30:00                                                                                                                                                                        
#$ -l thr=${nthreads}                                                                                                                               
#$ -l mem=3.9G                                                                                                                                                                                            
#$ -N viral_denovo_pipeline                                                                                                                                                                      
#$ -wd  ${output}                                                                                                                                                                               

cd \$TMPDIR                                                                                                                                                                                                 

module unload compilers/intel/11.1/072                                                                                                                                                                    
module unload mpi/qlogic/1.2.7/intel                                                                                                                                                                     
module unload mkl/10.2.5/035                                                                                                                                                                     
module load recommended/r                                                                                                                                                                                

export R_LIBS=${R}:$R_LIBS                                                                                                                                                                 
" >> $script

  	fi


# #################################################################  Collapsing                                                                                                                  
 echo "1) Collapse paired reads"
	
 summaryfile=${sample}_summary.txt


     echo "                                                                                                                                                       
${fastqCollapse2}  -i $seq1 $seq2 -o ${output_collapsed}/${sample} -sigLength ${collapseLev} -summary ${output_collapsed}/${summaryfile}                                                           
" >> $script



	



# # # # # ##########################################################################Merging                                                                                                                       
         echo "2) Merge Reads"
	
	echo "
${pear} -f ${output_collapsed}/${sample}_1.fq  -r ${output_collapsed}/${sample}_2.fq  -y 1G -j ${nthreads}  -o ${output_merged}/${sample} 

mv ${output_merged}/${sample}.unassembled.forward.fastq ${output_merged}/${sample}.1.fq
mv ${output_merged}/${sample}.unassembled.reverse.fastq ${output_merged}/${sample}.2.fq

" >> $script

	
		       
# # # ############################################### Trim adapters and low quality ends


  	echo "3) Trim adapters from reads"

	echo "
${trimG} -q 20 --paired  ${output_merged}/${sample}.1.fq ${output_merged}/${sample}.2.fq -o ${output_qc} 
${trimG} -q 20  ${output_merged}/${sample}.assembled.fastq  -o ${output_qc}
                                                                                                                                                                                                             
" >> $script


	echo "
${seqtk} seq -a ${output_qc}/${sample}.1_val_1.fq > ${output_qc}/${sample}_1.fasta
${seqtk} seq -a ${output_qc}/${sample}.2_val_2.fq > ${output_qc}/${sample}_2.fasta
${seqtk} seq -a ${output_qc}/${sample}.assembled_trimmed.fq > ${output_qc}/${sample}_assembled.fasta


" >> $script


# ###################################### blastn against pathogens datbase                                                        

    
      if [ ! -e $dbnr ]; then echo "File $dbnr does not exist."; exit; fi
                                                                                                                                                                                     
            echo -e "4) Submit blastn job against pathogen reference"



	  echo "                                                                                                                                                                                            

blastn_input_files=\"${output_qc}/${sample}_1.fasta  ${output_qc}/${sample}_2.fasta  ${output_qc}/${sample}_assembled.fasta\"

            echo "Input files for blastn are '$blastn_input_files'"


        for i in \$blastn_input_files 
         do
         echo "Print '$i'"

            if [ ! -e "'$i'" ]; then
                echo "Error, input file for blastn '$i' does not exist"    ##########Check input for blastn (created in previous step) exists                                                    
                exit



            else

                filename=\`basename \$i .fasta\`
                outputblastn=${output_blastn}/\${filename}_blastn.tab
 
                echo "Filename is '$filename'"
                echo "Output file for blastn is '$outputblastn'"

$blastn -db $db -query  \$i  -outfmt 6  -num_alignments 1   -culling_limit 1 -num_threads ${nthreads}  > \$outputblastn                  

         fi

	done



" >> $script


          echo -e "5)  keep  viral reads"

	  echo "
	 awk '{printf substr(\$0,1,length);getline;printf \"\t\"\$0;getline;getline;print \"\t\"\$0}' ${output_qc}/${sample}.1_val_1.fq |  awk 'BEGIN{OFS=FS=\"\t\"} {print \$1, \$2}' > ${output_qc}/${sample}_1.tab   
	 awk '{printf substr(\$0,1,length);getline;printf \"\t\"\$0;getline;getline;print \"\t\"\$0}' ${output_qc}/${sample}.2_val_2.fq |  awk 'BEGIN{OFS=FS=\"\t\"} {print \$1, \$2}' > ${output_qc}/${sample}_2.tab   
	 awk '{printf substr(\$0,1,length);getline;printf \"\t\"\$0;getline;getline;print \"\t\"\$0}' ${output_qc}/${sample}.assembled_trimmed.fq |  awk 'BEGIN{OFS=FS=\"\t\"} {print \$1, \$2}' > ${output_qc}/${sample}_assembled.tab   
	 
	 blastnfiles="${output_blastn}/*_blastn.tab"                                                                                                                                                    
	 for i in \$blastnfiles                                                                                                                                                                         
	   do                                                                                                                                                                                          
	   filename=\`basename \$i _blastn.tab\`                                                                                                                                                       
	   echo \$i                                                                                                                                                                                
	   echo \$filename                                               

	   awk -F\"\t\" '{print \"@\"\$1}' \$i > ${output_blastn}/\${filename}_readID.tab
	   awk  'BEGIN {OFS=\"\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]==\$1{print \$0}' ${output_blastn}/\${filename}_readID.tab ${output_qc}/\${filename}.tab | awk  'BEGIN{OFS=\"\t\"}{print \">\"\$1\" \"\$2\"\n\"\$3}' > ${output_blastn}/\${filename}.fasta


	 done

	 cat ${output_blastn}/${sample}_1_readID.tab ${output_blastn}/${sample}_2_readID.tab | sort | uniq > ${output_blastn}/${sample}_combined_readID.tab

	   awk  'BEGIN {OFS=\"\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]==\$1{print \$0}' ${output_blastn}/${sample}_combined_readID.tab ${output_qc}/${sample}_1.tab | awk  'BEGIN{OFS=\"\t\"}{print \">\"\$1\" \"\$2\"\n\"\$3}' > ${output_blastn}/${sample}_1.fasta
	   awk  'BEGIN {OFS=\"\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]==\$1{print \$0}' ${output_blastn}/${sample}_combined_readID.tab ${output_qc}/${sample}_2.tab | awk  'BEGIN{OFS=\"\t\"}{print \">\"\$1\" \"\$2\"\n\"\$3}' > ${output_blastn}/${sample}_2.fasta

rm ${output_blastn}/${sample}_1_readID.tab
rm ${output_blastn}/${sample}_2_readID.tab

" >> $script



 ######################### ASSEMBLY (with SPades)

     echo -e "6) Assembly with SPades (tries several parameters automatically)."


	echo "

${python} ${spades} -1 ${output_blastn}/${sample}_1.fasta -2 ${output_blastn}/${sample}_2.fasta  -s  ${output_blastn}/${sample}_assembled.fasta --only-assembler  --careful -t ${nthreads} -o ${output_assembly}

" >> $script


	echo "
${python} ${quast} ${output_assembly}/scaffolds.fasta -M 200  -o ${output_assembly}

" >> $script



##################### Pipeline branches out to 2 routines I)   FluA  II) for continuous genomes without looking at more than one strain when ordering contigs (e.g EBV, VZV etc)



 ########################################################################################################## SEGMENTED - FluA
 if [[ "$genome" == "FluA" ]]; then

  ########################################### CREATE SEGMENTS
 echo "
                                                                                                                                                                                                    
	export PERL5LIB=${PERL5LIB}:${Software}/bioperl-live

	inputContigs=${output_assembly}/${name}.fasta     


	perl $extractLargeContigs \${inputContigs} ${output_assembly} ${cutoff} ${name}

" >> $script


echo  "
	awk '/^>/ {printf(\"\n%s\n\",\$0);next; } { printf(\"%s\",\$0);}  END {printf(\"\n\");}' ${output_assembly}/${name}_grt${cutoff}.fa > ${output_assembly}/${name}1Line_grt${cutoff}.fasta

	awk 'BEGIN{RS=\">\"}NR>1{sub(\"\n\",\"\t\"); gsub(\"\n\",\"\"); print RS\$0}' ${output_assembly}/${name}1Line_grt${cutoff}.fasta | sed 's/>//g' > ${output_assembly}/${name}_grt${cutoff}.tab

" >> $script




echo " 


 for fluA_type in H1N1 H3N2; do

	$blastn -db  ~/scratchDir2/sofia/\${fluA_type}_db/\${fluA_type}_8segments.fasta -query ${output_assembly}/${name}_grt${cutoff}.fa -outfmt '6 qacc qstart qend pident sstart send evalue bitscore sacc' -evalue 1e-30 -culling_limit 1 -num_threads ${nthreads} > ${output_assembly}/out.blastn

 	echo \"done\"


	awk  'BEGIN {FS=OFS = \"\t\"}                                                                                                                                                                       
	{  	    max[\$1] = !(\$1 in max) ? \$8 : (\$8 > max[\$1]) ? \$8 : max[\$1]                                                                                                    
	    }                                                                                                                 
	    END {              
		for (i in max){            
			print i,max[i]                                                                                                                                             
		    }                                                                                                                                      
		}' ${output_assembly}/out.blastn > ${output_assembly}/max.blastn

awk 'BEGIN {FS=OFS=\"\t\" } {print \$1\"_\"\$2, \$0}' ${output_assembly}/max.blastn > ${output_assembly}/max.blastn.tmp
awk 'BEGIN {FS=OFS=\"\t\" } {print \$1\"_\"\$8, \$0}' ${output_assembly}/out.blastn > ${output_assembly}/out.blastn.tmp
awk 'BEGIN {FS=OFS=\"\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]==\$1{print \$0}' ${output_assembly}/max.blastn.tmp ${output_assembly}/out.blastn.tmp |  sort | uniq | awk  'BEGIN {FS=OFS=\"\t\"} {print \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10}'  > ${output_assembly}/out.blastn.max

awk 'BEGIN{FS=OFS=\"\t\"} NR==FNR{a[\$1]=\$2;next} (\$1) in a{print \$0, a[\$1]}' ${output_assembly}/${name}_grt${cutoff}.tab ${output_assembly}/out.blastn.max > ${output_assembly}/out
awk 'BEGIN{FS=OFS=\"\t\"}{print \$1,  \$5, \$6,  substr(\$10, \$2, (\$3-\$2+1)), \$9}' ${output_assembly}/out  > ${output_assembly}/out_seq_correct
awk 'BEGIN{OFS=FS=\"\t\"}{if (\$2<\$3) print \$0}' ${output_assembly}/out_seq_correct > ${output_assembly}/out_positive
awk 'BEGIN{OFS=FS=\"\t\"}{if (\$2>\$3) print \$0}' ${output_assembly}/out_seq_correct > ${output_assembly}/out_negative
awk 'BEGIN{OFS=FS=\"\t\"}{if (\$2>\$3) print \">\"\$1\"\n\"\$4}' ${output_assembly}/out_seq_correct > ${output_assembly}/negative.fa

awk 'BEGIN {                                                                                                                                                                                                
    j = n = split(\"A C G T\", t)                                                                                                                                                                             
    for (i = 0; ++i <= n;)                                                                                                                                                                                
    map[t[i]] = t[j--]                                                                                                                                                                                     
		}                                                                                                                                                                                                        
{                                                                                                                                                                                                          
    if (/^>/) print                                                                                                                                                                                       
    else {                                                                                                                                                                                                 
	    for (i = length; i; i--)                                                                                                                                                                     
	    printf \"%s\", map[substr(\$0, i, 1)]                                                                                                                                                           
	    print x                                                                                                                                                                                       
  	}                                                                                                                                                                                           
}'  ${output_assembly}/negative.fa  > ${output_assembly}/reverse.fasta


awk 'BEGIN{RS=\">\"}NR>1{sub(\"\n\",\"\t\"); gsub(\"\n\",\"\"); print RS\$0}' ${output_assembly}/reverse.fasta | sed 's/>//g' > ${output_assembly}/reverse.tab

awk 'BEGIN{FS=OFS=\"\t\"} NR==FNR{a[\$1]=\$2;next} (\$1) in a{print \$1,\$3,\$2, a[\$1], \$5}' ${output_assembly}/reverse.tab ${output_assembly}/out_negative > ${output_assembly}/out_negative_reversed


cat ${output_assembly}/out_positive ${output_assembly}/out_negative_reversed > ${output_assembly}/correct_contigs

sort -k5,5 -k2,3n  ${output_assembly}/correct_contigs  | awk 'BEGIN{OFS=OS=\"\t\"}{print \$1, \$2, \$3, \$4, \$5}' > ${output_assembly}/\${fluA_type}_ordered_contigs

 done



" >> $script





    echo "

R CMD BATCH  --no-save --restore --arg1=${output_assembly}/H1N1_ordered_contigs  --arg2=${output_assembly}/H3N2_ordered_contigs --arg3=${output_assembly}  ${stichContigs_segm} ${output_assembly}/stichContigs.out

#rm ${output_assembly}/*max*
rm ${output_assembly}/out*
rm ${output_assembly}/reverse*
rm ${output_assembly}/negative.fa
rm ${output_assembly}/correct_contigs
rm ${output_assembly}/ordered_contigs
rm ${output_assembly}/*grt250*



 " >> $script


#### Align and call variants against segments ###################################
         echo "7) Align against segments, call consnsus and call variants"
	
 	echo "

cat ${output_assembly}/H*segment1.fasta  ${output_assembly}/H*segment2.fasta  ${output_assembly}/H*segment3.fasta  ${output_assembly}/H*segment4.fasta  ${output_assembly}/H*segment5.fasta  ${output_assembly}/H*segment6.fasta  ${output_assembly}/H*segment7.fasta  ${output_assembly}/H*segment8.fasta >  ${output_assembly}/segments.fasta



###make index	
${novoindex} ${output_assembly}/${sample}.segments.novoindex ${output_assembly}/segments.fasta

" >> $script                                                                                                                                                                                               


#########remove duplicates and extract consensus by aligning reads against pseudogenome
	echo "
                                                                                                                                                                                               
############# Hardclipping and quality calibration                                                                                                                                             
${novoalign} -c ${nthreads} -t250 -H -r Random -a -k -e 500 -o SAM -F STDFQ -f ${output_qc}/${sample}.1_val_1.fq ${output_qc}/${sample}.2_val_2.fq -d ${output_assembly}/${sample}.segments.novoindex > ${output_assembly}/${sample}.sam

${novoalign} -c ${nthreads} -t180  -H -r Random -a -e 500  -o SAM -F STDFQ -f ${output_qc}/${sample}.assembled_trimmed.fq -d ${output_assembly}/${sample}.segments.novoindex > ${output_assembly}/${sample}_assembled.sam                                                                                                                                                                                                                                             
" >> $script                                                                                                                                                                                               


echo "


cat  ${output_assembly}/${sample}_assembled.sam  <(grep -v '^@'  ${output_assembly}/${sample}.sam)  > ${output_assembly}/${sample}_combined_tmp.sam

awk 'BEGIN{OFS=FS=\"\t\"}{ if (\$11 == \")\") {\$11 = \"*\"}; print }' ${output_assembly}/${sample}_combined_tmp.sam > ${output_assembly}/${sample}_combined.sam

${samtools} view -bS  -o ${output_assembly}/${sample}.bam ${output_assembly}/${sample}_combined.sam  


${samtools} sort  ${output_assembly}/${sample}.bam ${output_assembly}/${sample}_sorted                                                                                                  

${samtools} index  ${output_assembly}/${sample}_sorted.bam                                                                                                          

java -Xmx4g -jar ${picard} MarkDuplicates I=${output_assembly}/${sample}_sorted.bam O=${output_assembly}/${sample}_sorted_rmDUP.bam M=${output_assembly}/picard.dup.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT 

java -Xmx4g -jar ${picard}  BuildBamIndex I=${output_assembly}/${sample}_sorted_rmDUP.bam VALIDATION_STRINGENCY=SILENT

${samtools} depth ${output_assembly}/${sample}_sorted_rmDUP.bam > ${output_assembly}/${sample}.depth
                                                                                                          
java -Xmx4g -jar ${picard} SamToFastq INPUT=${output_assembly}/${sample}_sorted_rmDUP.bam FASTQ=${output_assembly}/${sample}_1.fq SECOND_END_FASTQ=${output_assembly}/${sample}_2.fq UNPAIRED_FASTQ=${output_assembly}/${sample}_unpaired.fq VALIDATION_STRINGENCY=SILENT                                                                                                                                                                                                             
${samtools} mpileup -A -f  ${output_assembly}/segments.fasta   ${output_assembly}/${sample}_sorted_rmDUP.bam   > ${output_assembly}/${sample}.pileup   

java -Xmx1000m -jar ${consensus} -i ${output_assembly}/${sample}.pileup  -r ${output_assembly}/segments.fasta -o ${output_assembly}/${sample}


###make index	
${novoindex} ${output_assembly}/${sample}.consensus.novoindex ${output_assembly}/${sample}.consensus.fa


" >> $script



#####call variants by aligning against consensus

	echo "
                                                                                                                                                                                               
############# Hardclipping and quality calibration                                                                                                                                             
${novoalign} -c ${nthreads} -t250 -r Random -H -a -k -e 500 -o SAM -F STDFQ -f ${output_assembly}/${sample}_1.fq ${output_assembly}/${sample}_2.fq -d ${output_assembly}/${sample}.consensus.novoindex > ${output_variants}/${sample}.sam

${novoalign} -c ${nthreads} -t180 -H -r Random -a -e 500  -o SAM -F STDFQ -f ${output_assembly}/${sample}_unpaired.fq -d ${output_assembly}/${sample}.consensus.novoindex > ${output_variants}/${sample}_assembled.sam                                                                                                                                                                                                                                             
" >> $script                                                                                                                                                                                             



 echo "

cat  ${output_variants}/${sample}_assembled.sam  <(grep -v '^@'  ${output_variants}/${sample}.sam)  > ${output_variants}/${sample}_combined_tmp.sam

awk 'BEGIN{OFS=FS=\"\t\"}{ if (\$11 == \")\") {\$11 = \"*\"}; print }' ${output_variants}/${sample}_combined_tmp.sam > ${output_variants}/${sample}_combined.sam


${samtools} view -bS -o ${output_variants}/${sample}.bam ${output_variants}/${sample}_combined.sam  

${samtools} sort ${output_variants}/${sample}.bam ${output_variants}/${sample}_sorted                                                                                                  

${samtools} index  ${output_variants}/${sample}_sorted.bam                                                                                                          

${samtools} mpileup -A -f  ${output_assembly}/${sample}.consensus.fa   ${output_variants}/${sample}_sorted.bam   > ${output_variants}/${sample}.pileup

java -Xmx1000m -jar ${varscan} pileup2cns ${output_variants}/${sample}.pileup  --min-reads2 5 --min-avg-qual 30 --min-var-freq 0.02 --variants 1 --p-value 0.01  > ${output_variants}/variants.tab

${samtools} depth ${output_variants}/${sample}_sorted.bam > ${output_variants}/${sample}.depth 

${samtools} depth ${output_variants}/${sample}_sorted.bam  |   awk '{sum+=\$3;count++}END{print \"Average = \"sum/count, \"\nTotal = \"sum}' > ${output_variants}/coverageSAMTOOLS.stats


" >> $script



fi





      if [[ "$genome" == "continuous" ]]; then


################## CREATE PSEUDOGENOME FROM CONTIGS
  echo "
                                                                                                                                                                                                    
 	export PERL5LIB=${PERL5LIB}:${Software}/bioperl-live

 	inputContigs=${output_assembly}/${name}.fasta     


 	perl $extractLargeContigs \${inputContigs} ${output_assembly} ${cutoff} ${name}

 " >> $script


echo  "
	awk '/^>/ {printf(\"\n%s\n\",\$0);next; } { printf(\"%s\",\$0);}  END {printf(\"\n\");}' ${output_assembly}/${name}_grt${cutoff}.fa > ${output_assembly}/${name}1Line_grt${cutoff}.fasta

	awk 'BEGIN{RS=\">\"}NR>1{sub(\"\n\",\"\t\"); gsub(\"\n\",\"\"); print RS\$0}' ${output_assembly}/${name}1Line_grt${cutoff}.fasta | sed 's/>//g' > ${output_assembly}/${name}_grt${cutoff}.tab

	$blastn -db ${dbSingle} -query ${output_assembly}/${name}_grt${cutoff}.fa -outfmt '6 qacc qstart qend pident sstart send evalue bitscore' -evalue 1e-30  -num_threads ${nthreads} > ${output_assembly}/out.blastn

 	echo \"done\"

" >> $script


echo "                                                                                                                                                                                                    

awk 'BEGIN{FS=OFS=\"\t\"} NR==FNR{a[\$1]=\$2;next} (\$1) in a{print \$0, a[\$1]}' ${output_assembly}/${name}_grt${cutoff}.tab ${output_assembly}/out.blastn > ${output_assembly}/out
awk 'BEGIN{FS=OFS=\"\t\"}{print \$1,  \$5, \$6,  substr(\$9, \$2, (\$3-\$2+1))}' ${output_assembly}/out  > ${output_assembly}/out_seq_correct
awk 'BEGIN{OFS=FS=\"\t\"}{if (\$2<\$3) print \$0}' ${output_assembly}/out_seq_correct > ${output_assembly}/out_positive
awk 'BEGIN{OFS=FS=\"\t\"}{if (\$2>\$3) print \$0}' ${output_assembly}/out_seq_correct > ${output_assembly}/out_negative

awk 'BEGIN {                                                                                                                                                                                                
    j = n = split(\"A C G T\", t)                                                                                                                                                                             
    for (i = 0; ++i <= n;)                                                                                                                                                                                
    map[t[i]] = t[j--]                                                                                                                                                                                     
		}                                                                                                                                                                                                        
{                                                                                                                                                                                                          
    if (/^>/) print                                                                                                                                                                                       
    else {                                                                                                                                                                                                 
	    for (i = length; i; i--)                                                                                                                                                                     
	    printf \"%s\", map[substr(\$0, i, 1)]                                                                                                                                                           
	    print \"\t\"\$1\"\t\"\$2\"\t\"\$3                                                                                                                                                                                        
  	}                                                                                                                                                                                           
}'  ${output_assembly}/out_negative  | awk 'BEGIN{OFS=FS=\"\t\"}{print \$2, \$4, \$3, \$1}'  > ${output_assembly}/out_reverse

cat ${output_assembly}/out_positive ${output_assembly}/out_reverse > ${output_assembly}/correct_contigs

sort -k2,3 -n ${output_assembly}/correct_contigs > ${output_assembly}/ordered_contigs

" >> $script



     echo "

R CMD BATCH  --no-save --restore --arg1=${output_assembly}/ordered_contigs --arg2=${output_assembly} --arg3=${sample}  ${stichContigs} ${output_assembly}/stichContigs.out


rm ${output_assembly}/out*
rm ${output_assembly}/reverse*
rm ${output_assembly}/negative.fa
rm ${output_assembly}/correct_contigs
rm ${output_assembly}/ordered_contigs
rm ${output_assembly}/*grt250*


  " >> $script








# #####Align, exctract consensus and call variants ###################################
          echo "7) Align against pseudogenome, extract consensus and call variants"
	
 	echo "

###make index	
${novoindex} ${output_assembly}/${sample}.pseudogenome.novoindex ${output_assembly}/pseudogenome.fasta

" >> $script                                                                                                                                                                                               


############ remove duplicates and extract consensus by aligning reads against pseudogenome
 	echo "
                                                                                                                                                                                               
# ############# Hardclipping and quality calibration                                                                                                                                             
 ${novoalign} -c ${nthreads} -t250 -r Random  -H -a -k -e 500 -o SAM -F STDFQ -f ${output_qc}/${sample}.1_val_1.fq ${output_qc}/${sample}.2_val_2.fq -d ${output_assembly}/${sample}.pseudogenome.novoindex > ${output_assembly}/${sample}.sam

 ${novoalign} -c ${nthreads} -t180 -r Random -H -a -e 500  -o SAM -F STDFQ -f ${output_qc}/${sample}.assembled_trimmed.fq -d ${output_assembly}/${sample}.pseudogenome.novoindex > ${output_assembly}/${sample}_assembled.sam  

 " >> $script                                                                                                                                                                                               




  echo "
cat  ${output_assembly}/${sample}_assembled.sam  <(grep -v '^@'  ${output_assembly}/${sample}.sam)  > ${output_assembly}/${sample}_combined_tmp.sam


awk 'BEGIN{OFS=FS=\"\t\"}{ if (\$11 == \")\") {\$11 = \"*\"}; print }' ${output_assembly}/${sample}_combined_tmp.sam > ${output_assembly}/${sample}_combined.sam


 ${samtools} view -bS -o ${output_assembly}/${sample}.bam ${output_assembly}/${sample}_combined.sam  


 ${samtools} sort  ${output_assembly}/${sample}.bam ${output_assembly}/${sample}_sorted                                                                                                  

 ${samtools} index  ${output_assembly}/${sample}_sorted.bam


java -Xmx4g -jar ${picard} MarkDuplicates I=${output_assembly}/${sample}_sorted.bam O=${output_assembly}/${sample}_sorted_rmDUP.bam M=${output_assembly}/picard.dup.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT 

java -Xmx4g -jar ${picard}  BuildBamIndex I=${output_assembly}/${sample}_sorted_rmDUP.bam  VALIDATION_STRINGENCY=SILENT

${samtools} depth ${output_assembly}/${sample}_sorted_rmDUP.bam > ${output_assembly}/${sample}.depth 
                                                                                                         
java -Xmx4g -jar ${picard} SamToFastq INPUT=${output_assembly}/${sample}_sorted_rmDUP.bam FASTQ=${output_assembly}/${sample}_1.fq SECOND_END_FASTQ=${output_assembly}/${sample}_2.fq UNPAIRED_FASTQ=${output_assembly}/${sample}_unpaired.fq VALIDATION_STRINGENCY=SILENT                                                                                                                                                                                                             
${samtools} mpileup -A -f  ${output_assembly}/pseudogenome.fasta   ${output_assembly}/${sample}_sorted_rmDUP.bam   > ${output_assembly}/${sample}.pileup   

java -Xmx1000m -jar ${consensus} -i ${output_assembly}/${sample}.pileup  -r ${output_assembly}/pseudogenome.fasta -o ${output_assembly}/${sample}


###make index	
${novoindex} ${output_assembly}/${sample}.consensus.novoindex ${output_assembly}/${sample}.consensus.fa

" >> $script



####call variants by aligning against consensus
 	echo "
                                                                                                                                                                                               
# ############# Hardclipping and quality calibration                                                                                                                                             
 ${novoalign} -c ${nthreads} -t250 -r Random -H -a -k -e 500 -o SAM -F STDFQ -f ${output_assembly}/${sample}_1.fq ${output_assembly}/${sample}_2.fq -d ${output_assembly}/${sample}.consensus.novoindex > ${output_variants}/${sample}.sam

 ${novoalign} -c ${nthreads} -t180 -H -r Random -a -e 500  -o SAM -F STDFQ -f ${output_assembly}/${sample}_unpaired.fq -d ${output_assembly}/${sample}.consensus.novoindex > ${output_variants}/${sample}_assembled.sam                                                                                                                                                                                                                                             
 " >> $script                                                                                                                                                                                             



  echo "

cat  ${output_variants}/${sample}_assembled.sam  <(grep -v '^@'  ${output_variants}/${sample}.sam)  > ${output_variants}/${sample}_combined_tmp.sam

awk 'BEGIN{OFS=FS=\"\t\"}{ if (\$11 == \")\") {\$11 = \"*\"}; print }' ${output_variants}/${sample}_combined_tmp.sam > ${output_variants}/${sample}_combined.sam


${samtools} view -bS -o ${output_variants}/${sample}.bam ${output_variants}/${sample}_combined.sam  

${samtools} sort ${output_variants}/${sample}.bam ${output_variants}/${sample}_sorted                                                                                                  

${samtools} index  ${output_variants}/${sample}_sorted.bam                                                                                                          

${samtools} mpileup -A -f  ${output_assembly}/${sample}.consensus.fa   ${output_variants}/${sample}_sorted.bam   > ${output_variants}/${sample}.pileup

java -Xmx1000m -jar ${varscan} pileup2cns ${output_variants}/${sample}.pileup  --min-reads2 5 --min-avg-qual 30 --min-var-freq 0.02 --variants 1 --p-value 0.01  > ${output_variants}/variants.tab

${samtools} depth ${output_variants}/${sample}_sorted.bam > ${output_variants}/${sample}.depth 


${samtools} depth ${output_variants}/${sample}_sorted.bam  |   awk '{sum+=\$3;count++}END{print \"Average = \"sum/count, \"\nTotal = \"sum}' > ${output_variants}/coverageSAMTOOLS.stats

 rm ${output_assembly}/*tmp*
 rm ${output_variants}/*tmp*
 rm ${output_variants}/${sample}.sam
 rm ${output_variants}/${sample}_assembled.sam
 rm ${output_assembly}/${sample}.sam
 rm ${output_assembly}/${sample}_assembled.sam


" >> $script



 fi



        qsub  -o $cluster_out -e $cluster_error   $script


