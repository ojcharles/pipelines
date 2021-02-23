#!/bin/bash
# HCMV reference based assembly pipeline
# Used to process ancient DNA SRA entries 
# Keen fastq, final bam, and output files
# v1 simple trim and bbmap - next rem human? decrease quality?
# OC 2021

#refs
# sra toolkit now supported here https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration

### pipeline variables
outdir=/archive/21_2_ancient_dna # no /
ref_seq=/archive/21_2_ancient_dna/ref/NC_006273.2.fasta # merlin
sra_list=${outdir}/sra_list3.tab # file with per line accession numbers
#sratoolkit=/SAN/breuerlab/pathseq1/oc/tools/sratoolkit.2.10.9-ubuntu64
sratoolkit=/usr/local/ncbi/sra-tools
in_qual=20 # quality cutoff
###


### begin
#---------------------------- activate conda env
source /opt/anaconda/bin/activate CMV_refbased
    
for sample in `cat $sra_list`; do
    echo $sample


    #---------------------------- check if already processed
    if [ -e ${sample}/refbased/${sample}.pileup ]
    then
        echo "${sample} sample already processed skipping ..."
        echo "${sample} sample already processed skipping ..." >> simple.log
        continue      # Skip rest of this particular loop iteration.
    fi

    ### end variables
    mkdir ${outdir}/${sample}/
    mkdir ${outdir}/${sample}/QC
    mkdir ${outdir}/${sample}/refbased




    #---------------------------- SRA -> fastq_1/_2 if not exist
    if [ -e ${sample}/${sample}_1.fastq ] # this will identify split or single fastq
    then
        echo "${sample} fastq present"
        echo "${sample} fastq present ..." >> simple.log
    else
        echo "${sample} fastq missing - downloading..."
        echo "${sample} fastq missing - downloading..." >> simple.log
        cd ${outdir}/${sample}
        ${sratoolkit}/bin/prefetch ${sample} -O ${outdir} -f yes
        #${sratoolkit}/bin/fastq-dump ${sample} --split-files --skip-technical # download directly
        ${sratoolkit}/bin/fastq-dump ${outdir}/${sample}/${sample}.sra --split-files --skip-technical
    fi
    # convert fastq to fasta for blast
    #paste - - - - < *_1.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > 1.fa
    cd ${outdir}/${sample}
    

    
    

    #---------------------------- fastq -> QC -> BAM handle split or single reads
    if [ -e ${sample}_2.fastq ] # if split else condier alone only
    then
        echo "Split read files"
        #---------------------------- QC
        trim_galore -q 20 -j 20 --paired ${outdir}/${sample}/${sample}*_1.fastq ${outdir}/${sample}/${sample}*_2.fastq  -o ${outdir}/${sample}/QC
        qc_file1=`find ${outdir}/${sample}/QC/${sample}_1*.fq`
        qc_file2=`find ${outdir}/${sample}/QC/${sample}_2*.fq`

        #---------------------------- mapping
        ###### reads -> bam
        bbmap.sh -Xmx8g in1=${qc_file1} in2=${qc_file2} minid=0.8 covstats=${outdir}/${sample}/refbased/${sample}.constats.txt out=${outdir}/${sample}/refbased/${sample}.sam ref=${ref_seq}  nodisk

    else
        echo "single read file"
        #---------------------------- QC
        trim_galore -q 20 -j 20 ${outdir}/${sample}/${sample}*_1.fastq -o ${outdir}/${sample}/QC
        qc_file1=`find ${outdir}/${sample}/QC/${sample}_1*.fq`

        #---------------------------- mapping
        ###### reads -> bam
        bbmap.sh -Xmx8g in=${qc_file1} minid=0.8 covstats=${outdir}/${sample}/refbased/${sample}.constats.txt out=${outdir}/${sample}/refbased/${sample}.sam ref=${ref_seq}  nodisk
    fi






    #---------------------------- assembly
    ###### sam -> bam       
    samtools view -bS  -o ${outdir}/${sample}/refbased/${sample}.bam ${outdir}/${sample}/refbased/${sample}.sam  


    ###### bam -> bam_sort
    samtools sort -m 1000000000 ${outdir}/${sample}/refbased/${sample}.bam -o ${outdir}/${sample}/refbased/${sample}_sorted.bam                                                                                                  

    ###### bam_sort index
    samtools index  ${outdir}/${sample}/refbased/${sample}_sorted.bam         

    ###### remove duplicates
    picard MarkDuplicates INPUT=${outdir}/${sample}/refbased/${sample}_sorted.bam  OUTPUT=${outdir}/${sample}/refbased/${sample}_accepted_hits.nodups.bam METRICS_FILE=${outdir}/${sample}/refbased/dupsv.txt REMOVE_DUPLICATES=TRUE AS=true VALIDATION_STRINGENCY=LENIENT
    
    ###### sort
    samtools sort -m 1000000000 ${outdir}/${sample}/refbased/${sample}_accepted_hits.nodups.bam -o ${outdir}/${sample}/refbased/${sample}_nodups_sorted.bam
    
    ###### index
    samtools index  ${outdir}/${sample}/refbased/${sample}_nodups_sorted.bam

    ###### bam_nodup_sort -> pileup
    samtools mpileup -A -f ${ref_seq} ${outdir}/${sample}/refbased/${sample}_nodups_sorted.bam >  ${outdir}/${sample}/refbased/${sample}.pileup




    #---------------------------- pileup -> consensus + variant outputs
    # generate consensus sequence
    python3 /home/oscar/tools/QUASR/extras/pileup_consensus.py -f ${outdir}/${sample}/refbased/${sample}.pileup -r ${ref_seq} -a 0.5 -o ${outdir}/${sample}/refbased/${sample} -d -l 1 # set the read depth to min of 0 required
    
    ### all positions with variants
    varscan mpileup2cns ${outdir}/${sample}/refbased/${sample}.pileup --min-reads2 1 --min-avg-qual 10 --min-var-freq 0.01  --p-value 0.5  --output-vcf 1 > ${outdir}/${sample}/refbased/${sample}_allpos_variants.vcf


    cd ${outdir} # raedy for next iter

    ### end of script messages
    echo "********************************* ${sample} COMPLETED *********************************"
    echo date
    echo "********************************* ${sample} COMPLETED *********************************"  >> simple.log
    
    
    
    # remove many intermediate files
    find ${outdir}/${sample}/QC -name "*.fq" -type f -delete
    find ${outdir}/${sample}/refbased/ -name "*.sam" -type f -delete
    find ${outdir}/${sample}/refbased/ -name "${sample}.bam" -type f -delete
    find ${outdir}/${sample}/refbased/ -name "${sample}_sorted.bam" -type f -delete
    find ${outdir}/${sample}/refbased/ -name "${sample}_accepted_hits.nodups.bam" -type f -delete
    #keep   ${sample}_nodups_sorted.bam


    
done
