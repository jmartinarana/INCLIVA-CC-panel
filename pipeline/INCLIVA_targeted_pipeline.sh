#!/bin/bash

function checkDirectories(){
        for dir in ${@}
        do
            if [ ! -d "${dir}" ]
            then
                mkdir -p ${dir}
            fi
        done
}

function get_software_version(){
    server_info=`uname -a`
    echo ${server_info} | awk '{print "Server="$0}' >> ${3}
    lscpu >> ${3}
    version_fastqc=`fastqc -v | awk '{print $2}'`
    echo "FASTQC=${version_fastqc}" >> ${1}
    version_inc_UMI=`inc_include_UMI_in_fastq_header.py -v | awk '{print $2}'`
    echo "inc_include_UMI_in_fastq_header.py=${version_inc_UMI}" >> ${1}
    version_prinseq=`prinseq-lite.pl -version`
    echo "PRINSEQ=${version_prinseq}" >> ${1}
    version_bwa=`bwa |& grep Version`
    echo "BWA=${version_bwa}" >> ${1}
    version_samtools=`samtools |& grep Version`
    echo "SAMTOOLS=${version_samtools}" >> ${1}
    if [ "${2}" == "yes" ]
    then
        the_path=`which umi_tools`
        version_umi=`ls -latr ${the_path} | awk '{print $11}' | sed 's/\//\t/gi' | awk '{print $4}'`
        echo "UMI_TOOLS=${version_umi}" >> ${1}
    else
        echo "PICARD=picardv2.18.6" >> ${1}
    fi
    echo "GATK=gatk-4.0.5.0"
    version_bedtools=`bedtools --version`
    echo "BEDTOOLS=${version_bedtools}" >> ${1}
    version_tabix_gzip=`tabix |& grep Version`
    echo "TABIX-BGZIP=${version_tabix_gzip}" >> ${1}
    version_bcftools=`bcftools |& grep Version`
    echo "BCFTOOLS=${version_bcftools}" >> ${1}
    path_vcftools=`which vcftools`
    version_vcftools=`ls -latr ${path_vcftools} | awk '{print $11}' | sed 's/\//\t/gi' | awk '{print $4}'`
    echo "VCFTOOLS=${version_vcftools}" >> ${1}
    path_vardict=`which VarDict`
    version_vardict=`ls -latr ${path_vardict} | awk '{print $11}' | sed 's/\//\t/gi' | awk '{print $4}'`
    echo "VARDICT=${version_vardict}" >> ${1}
    inc_select_variants_version=`inc_select_variants.py &> ${4}/kk && tail -2 ${4}/kk | head -1`
    echo "version_inc_selectvariants=${inc_select_variants_version}" >> ${1}
    path_vcflib=`which vcfuniq`
    version_vcflib=`ls -latr ${path_vcflib} | awk '{print $11}' | sed 's/\//\t/gi' | awk '{print $4}'`
    echo "VCFLIB=${version_vcflib}" >> ${1}
    echo "PICARD=picardv2.18.6" >> ${1}
    inc_iqr=`inc_calculate_IQR_and_coverage_plot.R &> ${4}/kk | head -1 | awk '{print $2}'`
    echo "inc_calculate_IQR_and_coverage_plot.R=${inc_iqr}" >> ${1}
    echo "SNPSIFT=snpSift_v4.3t" >> ${1}
    inc_mv_cosmic_version=`inc_mv_cosmic_id_to_annotation.py &> ${4}/kk && tail -2 ${4}/kk | head -1`
    echo "inc_mv_cosmic_id_to_annotation.py=${inc_mv_cosmic_version}" >> ${1}
    inc_vcftocsv_version=`inc_vcf_to_csv.py &> ${4}/kk && tail -2 ${4}/kk | head -1`
    echo "inc_vcftocsv_version.py=${inc_vcftocsv_version}" >> ${1}
}


function fastq_preprocessing(){
    echo "Running fastq preprocessing"
    date
    echo "fastqc -o ${5}/QC/rawdata -t ${6} ${1} ${2}"
    fastqc -o ${5}/QC/rawdata -t ${6} ${1} ${2}
    date
    echo "Sample has UMIs"
    gunzip --force ${1}
    gunzip --force ${2}
    R1=`echo ${1} | sed 's/.gz//'`
    R2=`echo ${2} | sed 's/.gz//'`
    
    echo "prinseq-lite.pl -fastq ${R1} -fastq2 ${R2} -min_qual_mean 30 -out_good ${4}/QC/filtered_${3} -out_bad ${4}/QC/bad_${3}"
    date
    prinseq-lite.pl -fastq ${R1} -fastq2 ${R2} -min_qual_mean 30 -out_good ${4}/QC/filtered_${3} -out_bad ${4}/QC/bad_${3}
    echo "fastqc -o ${5}/QC/rawdata -t ${6} ${4}/QC/filtered_${3}_1.fastq ${4}/QC/filtered_${3}_2.fastq"
    date
    fastqc -o ${5}/QC/rawdata -t ${6} ${4}/QC/filtered_${3}_1.fastq ${4}/QC/filtered_${3}_2.fastq
    date
    gzip ${R1}
    gzip ${R2}
}

function mapping_and_bampostprocessing(){
    echo "fgbio-1.1.0 FastqToBam --input ${1} ${2} --read-structures +T 12M11S+T -o ${3}/mapping/${6}_tmp.bam --sample ${6} --library Illumina"
    fgbio-1.1.0 FastqToBam --input ${1} ${2} --read-structures +T 12M11S+T -o ${3}/mapping/${6}_tmp.bam --sample ${6} --library Illumina

    echo "samtools fastq -1 ${3}/QC/${6}_fullcrum_R1.fastq -2 ${3}/QC/${6}_fullcrum_R2.fastq -T RX ${3}/mapping/${6}_tmp.bam"
    samtools fastq -1 ${3}/QC/${6}_fullcrum_R1.fastq -2 ${3}/QC/${6}_fullcrum_R2.fastq -T RX ${3}/mapping/${6}_tmp.bam

    echo "bwa mem -T 0 -M -C -R "@RG\tID:AgilentXTHS\tSM:${6}\tPL:Illumina" -t ${5} ${7} ${3}/QC/${6}_fullcrum_R1.fastq ${3}/QC/${6}_fullcrum_R2.fastq | samtools view -S -b -h | samtools sort --threads ${5} > ${3}/mapping/${6}_sorted.bam"
    bwa mem -T 0 -M -C -R "@RG\tID:AgilentXTHS\tSM:${6}\tPL:Illumina" -t ${5} ${7} ${3}/QC/${6}_fullcrum_R1.fastq ${3}/QC/${6}_fullcrum_R2.fastq | samtools view -S -b -h | samtools sort --threads ${5} > ${3}/mapping/${6}_sorted.bam
    samtools index ${3}/mapping/${6}_sorted.bam
    
    samtools view -h -L ${8}/capture.bed -b ${3}/mapping/${6}_sorted.bam > ${3}/mapping/ontarget_${6}_sorted.bam
    samtools index ${3}/mapping/ontarget_${6}_sorted.bam

    echo "picardv2.18.6.sh FixMateInformation I=${3}/mapping/${6}_sorted.bam O=${3}/mapping/${6}_mate.bam"
    picardv2.18.6.sh FixMateInformation I=${3}/mapping/${6}_sorted.bam O=${3}/mapping/${6}_mate.bam

    echo "fgbio-1.1.0 GroupReadsByUmi -i ${3}/mapping/${6}_mate.bam -o ${3}/mapping/${6}_group.bam -s adjacency"
    fgbio-1.1.0 GroupReadsByUmi -i ${3}/mapping/${6}_mate.bam -o ${3}/mapping/${6}_group.bam -s adjacency

    echo "fgbio-1.1.0 CallMolecularConsensusReads -i ${3}/mapping/${6}_group.bam -o ${3}/mapping/${6}_consensus.bam -M 2"
    fgbio-1.1.0 CallMolecularConsensusReads -i ${3}/mapping/${6}_group.bam -o ${3}/mapping/${6}_consensus.bam -M 2

    echo "samtools fastq -1 ${3}/QC/${6}_processed_R1.fastq -2 ${3}/QC/${6}_processed_R2.fastq -T cD,cE,MI,cM,RX ${3}/mapping/${6}_consensus.bam"
    samtools fastq -1 ${3}/QC/${6}_processed_R1.fastq -2 ${3}/QC/${6}_processed_R2.fastq -T cD,cE,MI,cM,RX ${3}/mapping/${6}_consensus.bam

    echo "bwa mem -T 0 -M -C -R "@RG\tID:AgilentXTHS\tSM:${6}\tPL:Illumina" -t ${5} ${7} ${3}/QC/${6}_processed_R1.fastq ${3}/QC/${6}_processed_R2.fastq | samtools view -S -b -h | samtools sort --threads ${5} > ${3}/mapping/${6}_processed.bam"
    bwa mem -T 0 -M -C -R "@RG\tID:AgilentXTHS\tSM:${6}\tPL:Illumina" -t ${5} ${7} ${3}/QC/${6}_processed_R1.fastq ${3}/QC/${6}_processed_R2.fastq | samtools view -S -b -h | samtools sort --threads ${5} > ${3}/mapping/${6}_processed.bam

    echo "fgbio-1.1.0 ClipBam -i ${3}/mapping/${6}_processed.bam -o ${3}/mapping/${6}_clipped.bam -r ${7} --read-one-five-prime=1"
    fgbio-1.1.0 ClipBam -i ${3}/mapping/${6}_processed.bam -o ${3}/mapping/${6}_clipped.bam -r ${7} --read-one-five-prime=1
       
    echo "picardv2.18.6.sh CleanSam I=${3}/mapping/${6}_clipped.bam O=${3}/mapping/${sample}_processed_tmp2.bam"
    picardv2.18.6.sh CleanSam I=${3}/mapping/${6}_clipped.bam O=${3}/mapping/${sample}_processed_tmp2.bam
    echo "samtools view -h ${3}/mapping/${sample}_processed_tmp2.bam | awk 'function abs(v) {return v < 0 ? -v : v} {if ((abs($9) > 100 && length($10) > 100) || $1 ~ /^@/) print $0}' | samtools view -S -b -h > ${3}/mapping/${sample}_processed.bam"
    samtools view -h ${3}/mapping/${sample}_processed_tmp2.bam | awk 'function abs(v) {return v < 0 ? -v : v} {if ((abs($9) > 100 && length($10) > 100) || $1 ~ /^@/) print $0}' | samtools view -S -b -h > ${3}/mapping/${sample}_processed.bam

    #### Deberia anadir aqui FixMateInformation
    samtools index ${3}/mapping/${sample}_processed.bam
    samtools view -b -h -L ${8}/capture.bed ${3}/mapping/${sample}_processed.bam > ${3}/mapping/ontarget_${sample}_processed.bam

    echo "gatk-4.0.5.0 --java-options "-Xmx30g" BaseRecalibrator -I ${3}/mapping/${6}_processed.bam -R ${7} --intervals ${8}/capture.bed --known-sites ${8}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites ${8}/All_20180418.vcf.gz --output ${3}/variant_calling/tmp_${6}.table"
    /home/jmartin/gatk/gatk-4.0.11.0/gatk --java-options "-Xmx30g" BaseRecalibrator -I ${3}/mapping/${6}_processed.bam -R ${7} --intervals ${8}/capture.bed --known-sites ${8}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites ${8}/All_20180418.vcf.gz --output ${3}/variant_calling/tmp_${6}.table
    echo "gatk-4.0.5.0 --java-options "-Xmx30g" ApplyBQSR -R ${7} -I ${3}/mapping/${6}_processed.bam --bqsr-recal-file ${3}/variant_calling/tmp_${6}.table -O ${4}/mapping/${6}.bam"
    /home/jmartin/gatk/gatk-4.0.11.0/gatk --java-options "-Xmx30g" ApplyBQSR -R ${7} -I ${3}/mapping/${6}_processed.bam --bqsr-recal-file ${3}/variant_calling/tmp_${6}.table -O ${4}/mapping/${6}.bam
    date
}

function variant_calling(){
    bedtools sort -i ${7}/target.bed | bedtools merge > ${2}/mapping/target.bed
    samtools depth -d 9898989898 -a -b ${2}/mapping/target.bed ${1} > ${2}/mapping/${5}_depth
    awk '{if ($3 >= '$8') print $1"\t"$2-1"\t"$2}' ${2}/mapping/${5}_depth | bedtools merge > ${3}/QC/mapping/${5}_covered_regions_${8}x.bed
    lofreq call-parallel -f ${6} -l ${3}/QC/mapping/${5}_covered_regions_${8}x.bed --pp-threads 20 --force-overwrite -o ${2}/variant_calling/tmp1_lofreq_${5}.vcf ${1}

    date

    #grep "#" ${2}/variant_calling/tmp1_lofreq_${5}.vcf > ${2}/variant_calling/tmp1_lofreq_${5}_tmp.vcf
    #grep -v "#" ${2}/variant_calling/tmp1_lofreq_${5}.vcf | awk '$6 > 1000 {print $0}' >> ${2}/variant_calling/tmp1_lofreq_${5}_tmp.vcf
    echo "bedtools intersect -header -a ${2}/variant_calling/tmp1_gatk_${5}.vcf -b ${3}/QC/mapping/${5}_covered_regions_${8}x.bed > ${2}/variant_calling/tmp2_gatk_${5}.vcf"
    bedtools intersect -header -a ${2}/variant_calling/tmp1_lofreq_${5}.vcf -b ${3}/QC/mapping/${5}_covered_regions_${8}x.bed > ${2}/variant_calling/tmp2_lofreq_${5}.vcf
    bgzip -f ${2}/variant_calling/tmp2_lofreq_${5}.vcf
    tabix -p vcf -f ${2}/variant_calling/tmp2_lofreq_${5}.vcf.gz
    echo "bcftools norm -m-any ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/gatk_${5}.vcf"
    #bcftools norm -o ${2}/variant_calling/gatk_${5}.vcf -m -both -c w -f ${6} ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz
    bcftools norm -m-any ${2}/variant_calling/tmp2_lofreq_${5}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/lofreq_${5}.vcf
    bgzip -f ${2}/variant_calling/lofreq_${5}.vcf
    tabix -p vcf -f ${2}/variant_calling/lofreq_${5}.vcf.gz
    date
   
    echo "inc_select_variants_lofreq.py -i ${2}/variant_calling/lofreq_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp_lofreq.vcf -v 0.02 -d 10"

    inc_select_variants_lofreq.py -i ${2}/variant_calling/lofreq_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp_lofreq.vcf -v 0.000000000000000000000000002 -d 10
    echo "vcf-sort ${2}/variant_calling/${5}_tmp.vcf | vcfuniq > ${2}/variant_calling/${5}.vcf"
    vcf-sort ${2}/variant_calling/${5}_tmp_lofreq.vcf | vcfuniq > ${2}/variant_calling/${5}.vcf
    #vcfuniq ${2}/variant_calling/${5}_tmp.vcf > ${2}/variant_calling/${5}.vcf
    bgzip -f ${2}/variant_calling/${5}.vcf && tabix -f -p vcf ${2}/variant_calling/${5}.vcf.gz
    python /home/jmartin/scripts/anno_format_lofreq.py ${2}/variant_calling/${5}.vcf.gz ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf
    bgzip -f ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf && tabix -f -p vcf ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf.gz
    python /usr/local/bin/anno_bam_info.py ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf.gz ${2}/variant_calling/${5}_lofreq_anno.vcf Clipped
    bgzip -f ${2}/variant_calling/${5}_lofreq_anno.vcf && tabix -f -p vcf ${2}/variant_calling/${5}_lofreq_anno.vcf.gz
    mv ${2}/variant_calling/${5}.vcf.gz ${2}/variant_calling/${5}_lofreq.vcf.gz
    rm ${2}/variant_calling/${5}.vcf.gz.tbi
    date
}


function variant_calling2(){
   
    lofreq call-parallel -f ${6} -l ${3}/QC/mapping/${10}_covered_regions_${8}x.bed --pp-threads 20 --force-overwrite -o ${2}/variant_calling/tmp1_lofreq_${10}.vcf ${1}

    date
    #grep "#" ${2}/variant_calling/tmp1_lofreq_${10}.vcf > ${2}/variant_calling/tmp1_lofreq_${10}_tmp.vcf
    #grep -v "#" ${2}/variant_calling/tmp1_lofreq_${10}.vcf | awk '$6 > 1000 {print $0}' >> ${2}/variant_calling/tmp1_lofreq_${10}_tmp.vcf
    echo "bedtools intersect -header -a ${2}/variant_calling/tmp1_gatk_${10}.vcf -b ${3}/QC/mapping/${10}_covered_regions_${8}x.bed > ${2}/variant_calling/tmp2_gatk_${10}.vcf"
    bedtools intersect -header -a ${2}/variant_calling/tmp1_lofreq_${10}.vcf -b ${3}/QC/mapping/${10}_covered_regions_${8}x.bed > ${2}/variant_calling/tmp2_lofreq_${10}.vcf
    bgzip -f ${2}/variant_calling/tmp2_lofreq_${10}.vcf
    tabix -p vcf -f ${2}/variant_calling/tmp2_lofreq_${10}.vcf.gz
    echo "bcftools norm -m-any ${2}/variant_calling/tmp2_gatk_${10}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/gatk_${10}.vcf"
    #bcftools norm -o ${2}/variant_calling/gatk_${10}.vcf -m -both -c w -f ${6} ${2}/variant_calling/tmp2_gatk_${10}.vcf.gz
    bcftools norm -m-any ${2}/variant_calling/tmp2_lofreq_${10}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/lofreq_${10}.vcf
    bgzip -f ${2}/variant_calling/lofreq_${10}.vcf
    tabix -p vcf -f ${2}/variant_calling/lofreq_${10}.vcf.gz
    date
    
    echo "inc_select_variants_lofreq.py -i ${2}/variant_calling/lofreq_${10}.vcf.gz -o ${2}/variant_calling/${10}_tmp_lofreq.vcf -v 0.02 -d 10"

    inc_select_variants_lofreq.py -i ${2}/variant_calling/lofreq_${10}.vcf.gz -o ${2}/variant_calling/${10}_tmp_lofreq.vcf -v 0.00000000000000000000000000002 -d 10
    echo "vcf-sort ${2}/variant_calling/${10}_tmp.vcf | vcfuniq > ${2}/variant_calling/${10}.vcf"
    vcf-sort ${2}/variant_calling/${10}_tmp_lofreq.vcf | vcfuniq > ${2}/variant_calling/${10}.vcf
    #vcfuniq ${2}/variant_calling/${10}_tmp.vcf > ${2}/variant_calling/${10}.vcf
    bgzip -f ${2}/variant_calling/${10}.vcf && tabix -f -p vcf ${2}/variant_calling/${10}.vcf.gz
    python /home/jmartin/scripts/anno_format_lofreq.py ${2}/variant_calling/${10}.vcf.gz ${2}/variant_calling/${10}_lofreq_anno_tmp.vcf
    bgzip -f ${2}/variant_calling/${10}_lofreq_anno_tmp.vcf && tabix -f -p vcf ${2}/variant_calling/${10}_lofreq_anno_tmp.vcf.gz
    python /usr/local/bin/anno_bam_info.py ${2}/variant_calling/${10}_lofreq_anno_tmp.vcf.gz ${2}/variant_calling/${10}_lofreq_anno.vcf Clipped
    bgzip -f ${2}/variant_calling/${10}_lofreq_anno.vcf && tabix -f -p vcf ${2}/variant_calling/${10}_lofreq_anno.vcf.gz
    mv ${2}/variant_calling/${10}.vcf.gz ${2}/variant_calling/${10}_lofreq.vcf.gz
    rm ${2}/variant_calling/${10}.vcf.gz.tbi


    #echo "vcf-merge -c none ${2}/variant_calling/${5}_lofreq.vcf.gz ${2}/variant_calling/${10}_lofreq.vcf.gz | bgzip > ${2}/variant_calling/${5}_lofreq_merge.vcf.gz"
    #vcf-merge -c none ${2}/variant_calling/${5}_lofreq_anno.vcf.gz ${2}/variant_calling/${10}_lofreq_anno.vcf.gz | bgzip > ${2}/variant_calling/${5}_lofreq_merge.vcf.gz
    #tabix -p vcf -f ${2}/variant_calling/${5}_lofreq_merge.vcf.gz
    date
}


function variant_calling3(){
    
    echo "gatk-4.0.5.0 --java-options "-Xmx30g" Mutect2 --native-pair-hmm-threads ${4} -R ${6} -I ${1} -L ${3}/QC/mapping/${5}_covered_regions_${8}x.bed -tumor ${5} -O ${2}/variant_calling/tmp1_gatk_${5}.vcf --max-reads-per-alignment-start 0"
    gatk-4.0.5.0 --java-options "-Xmx30g" Mutect2 --native-pair-hmm-threads ${4} -R ${6} -I ${1} -L ${3}/QC/mapping/${5}_covered_regions_${8}x.bed -tumor ${5} -O ${2}/variant_calling/tmp1_gatk_${5}.vcf --max-reads-per-alignment-start 0
    gatk-4.0.5.0 --java-options "-Xmx30g" SelectVariants -R ${6} -V ${2}/variant_calling/tmp1_gatk_${5}.vcf -O ${2}/variant_calling/tmp2_gatk_${5}.vcf --select-type-to-include INDEL

    bgzip -f ${2}/variant_calling/tmp2_gatk_${5}.vcf
    tabix -p vcf -f ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz

    date

    echo "bcftools norm -m-any ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/gatk_${5}.vcf"
    bcftools norm -m-any ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/mutect2_${5}.vcf
    bgzip -f ${2}/variant_calling/mutect2_${5}.vcf
    tabix -p vcf -f ${2}/variant_calling/mutect2_${5}.vcf.gz
    date
    
    echo "inc_select_variants.py -i ${2}/variant_calling/mutect2_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp_mutect2.vcf -v 0.02 -d ${8}"

    inc_select_variants.py -i ${2}/variant_calling/mutect2_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp_mutect2.vcf -v 0.000000000000000000000000002 -d ${8}


    echo "${2}/variant_calling/${5}_tmp_mutect2.vcf | vcfuniq > ${2}/variant_calling/${5}_tmp.vcf"
    vcf-sort ${2}/variant_calling/${5}_tmp_mutect2.vcf | vcfuniq > ${2}/variant_calling/${5}_tmp2_mutect2.vcf

    #vcfuniq ${2}/variant_calling/${5}_tmp.vcf > ${2}/variant_calling/${5}.vcf
    bgzip -f ${2}/variant_calling/${5}_tmp2_mutect2.vcf && tabix -f -p vcf ${2}/variant_calling/${5}_tmp2_mutect2.vcf.gz
    
    python /usr/local/bin/anno_bam_info.py ${2}/variant_calling/${5}_tmp2_mutect2.vcf.gz ${2}/variant_calling/${5}_mutect2.vcf Clipped
    bgzip -f ${2}/variant_calling/${5}_mutect2.vcf && tabix -f -p vcf ${2}/variant_calling/${5}_mutect2.vcf.gz

    date



    echo "gatk-4.0.5.0 --java-options "-Xmx30g" Mutect2 --native-pair-hmm-threads ${4} -R ${6} -I ${9} -L ${3}/QC/mapping/${10}_covered_regions_${8}x.bed -tumor ${10} -O ${2}/variant_calling/tmp1_gatk_${10}.vcf --max-reads-per-alignment-start 0"
    gatk-4.0.5.0 --java-options "-Xmx30g" Mutect2 --native-pair-hmm-threads ${4} -R ${6} -I ${9} -L ${3}/QC/mapping/${10}_covered_regions_${8}x.bed -tumor ${10} -O ${2}/variant_calling/tmp1_gatk_${10}.vcf --max-reads-per-alignment-start 0
    gatk-4.0.5.0 --java-options "-Xmx30g" SelectVariants -R ${6} -V ${2}/variant_calling/tmp1_gatk_${10}.vcf -O ${2}/variant_calling/tmp2_gatk_${10}.vcf --select-type-to-include INDEL

    bgzip -f ${2}/variant_calling/tmp2_gatk_${10}.vcf
    tabix -p vcf -f ${2}/variant_calling/tmp2_gatk_${10}.vcf.gz

    date

    echo "bcftools norm -m-any ${2}/variant_calling/tmp2_gatk_${10}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/gatk_${10}.vcf"
    bcftools norm -m-any ${2}/variant_calling/tmp2_gatk_${10}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/mutect2_${10}.vcf
    bgzip -f ${2}/variant_calling/mutect2_${10}.vcf
    tabix -p vcf -f ${2}/variant_calling/mutect2_${10}.vcf.gz
    date
    
    echo "inc_select_variants.py -i ${2}/variant_calling/mutect2_${10}.vcf.gz -o ${2}/variant_calling/${10}_tmp_mutect2.vcf -v 0.02 -d ${8}"

    inc_select_variants.py -i ${2}/variant_calling/mutect2_${10}.vcf.gz -o ${2}/variant_calling/${10}_tmp_mutect2.vcf -v 0.000000000000000000000000002 -d ${8}


    echo "${2}/variant_calling/${10}_tmp_mutect2.vcf | vcfuniq > ${2}/variant_calling/${10}_tmp.vcf"
    vcf-sort ${2}/variant_calling/${10}_tmp_mutect2.vcf | vcfuniq > ${2}/variant_calling/${10}_tmp2_mutect2.vcf

    #vcfuniq ${2}/variant_calling/${10}_tmp.vcf > ${2}/variant_calling/${10}.vcf
    bgzip -f ${2}/variant_calling/${10}_tmp2_mutect2.vcf && tabix -f -p vcf ${2}/variant_calling/${10}_tmp2_mutect2.vcf.gz
    
    python /usr/local/bin/anno_bam_info.py ${2}/variant_calling/${10}_tmp2_mutect2.vcf.gz ${2}/variant_calling/${10}_mutect2.vcf Clipped
    bgzip -f ${2}/variant_calling/${10}_mutect2.vcf && tabix -f -p vcf ${2}/variant_calling/${10}_mutect2.vcf.gz

    date
}



function qiaseq(){
    mkdir -p ${2}/variant_calling/qiaseq
    cd ${2}/variant_calling/qiaseq
    python /usr/local/src/qiaseq/qiaseq-dna_V2/run_qiaseq_dna_prueba.py /media/scratch3/UMP_FIS_2019/configs/${5}/run_sm_counter_v2.params_${5}.txt v2 tumor-normal ${5} ${10} > ${2}/QC/${5}_qiaseq_run.log
    
    date

    #grep "#" ${2}/variant_calling/qiaseq/${5}.smCounter.cplx.vcf > ${2}/variant_calling/qiaseq/${5}.smCounter.cplx_tmp.vcf
    #grep -v "#" ${2}/variant_calling/qiaseq/${5}.smCounter.cplx.vcf | awk '$6 > 100 {print $0}' >> ${2}/variant_calling/qiaseq/${5}.smCounter.cplx_tmp.vcf
    #grep "#" ${2}/variant_calling/qiaseq/${10}.smCounter.cplx.vcf > ${2}/variant_calling/qiaseq/${10}.smCounter.cplx_tmp.vcf
    #grep -v "#" ${2}/variant_calling/qiaseq/${10}.smCounter.cplx.vcf | awk '$6 > 100 {print $0}' >> ${2}/variant_calling/qiaseq/${10}.smCounter.cplx_tmp.vcf

    echo "bedtools intersect -header -a ${2}/variant_calling/qiaseq/${5}.smCounter.cplx.vcf -b ${3}/QC/mapping/${5}_covered_regions_${8}x.bed > ${2}/variant_calling/qiaseq/${5}.smCounter.tmp.vcf"
    bedtools intersect -header -a ${2}/variant_calling/qiaseq/${5}.smCounter.cplx.vcf -b ${3}/QC/mapping/${5}_covered_regions_${8}x.bed > ${2}/variant_calling/qiaseq/${5}.smCounter.tmp.vcf
    echo "bedtools intersect -header -a ${2}/variant_calling/qiaseq/${10}.smCounter.cplx.vcf -b ${3}/QC/mapping/${10}_covered_regions_${8}x.bed > ${2}/variant_calling/qiaseq/${10}.smCounter.tmp.vcf"
    bedtools intersect -header -a ${2}/variant_calling/qiaseq/${10}.smCounter.cplx.vcf -b ${3}/QC/mapping/${10}_covered_regions_${8}x.bed > ${2}/variant_calling/qiaseq/${10}.smCounter.tmp.vcf

    bgzip -f ${2}/variant_calling/qiaseq/${5}.smCounter.tmp.vcf
    tabix -p vcf -f ${2}/variant_calling/qiaseq/${5}.smCounter.tmp.vcf.gz

    bgzip -f ${2}/variant_calling/qiaseq/${10}.smCounter.tmp.vcf
    tabix -p vcf -f ${2}/variant_calling/qiaseq/${10}.smCounter.tmp.vcf.gz

    echo "bcftools norm -m-any ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/gatk_${5}.vcf"
    #bcftools norm -o ${2}/variant_calling/qiaseq/${5}.smCounter.tmp.vcf.gz -m -both -c w -f ${6} ${2}/variant_calling/qiaseq/tmp2_smcounter_${5}.vcf.gz
    bcftools norm -m-any ${2}/variant_calling/qiaseq/${5}.smCounter.tmp.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/qiaseq/smcounter_${5}.vcf

    #bcftools norm -o ${2}/variant_calling/qiaseq/${10}.smCounter.tmp.vcf.gz -m -both -c w -f ${6} ${2}/variant_calling/qiaseq/tmp2_smcounter_${10}.vcf.gz
    bcftools norm -m-any ${2}/variant_calling/qiaseq/${10}.smCounter.tmp.vcf.gz | bcftools norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/qiaseq/smcounter_${10}.vcf

    bgzip -f ${2}/variant_calling/qiaseq/smcounter_${5}.vcf
    tabix -p vcf -f ${2}/variant_calling/qiaseq/smcounter_${5}.vcf.gz
    bgzip -f ${2}/variant_calling/qiaseq/smcounter_${10}.vcf
    tabix -p vcf -f ${2}/variant_calling/qiaseq/smcounter_${10}.vcf.gz
    date  
    
    echo "inc_select_variants.py -i ${2}/variant_calling/gatk_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp.vcf -v 0.02 -d ${8}"
    inc_select_variants.py -i ${2}/variant_calling/qiaseq/smcounter_${5}.vcf.gz -o ${2}/variant_calling/qiaseq/${5}_smcounter_tmp.vcf -v 0.0000000000000000000000002 -d 20
    inc_select_variants.py -i ${2}/variant_calling/qiaseq/smcounter_${10}.vcf.gz -o ${2}/variant_calling/qiaseq/${10}_smcounter_tmp.vcf -v 0.000000000000000000000002 -d 20


    echo "vcf-sort ${2}/variant_calling/${5}_tmp.vcf | vcfuniq > ${2}/variant_calling/${5}.vcf"
    vcf-sort ${2}/variant_calling/qiaseq/${5}_smcounter_tmp.vcf | vcfuniq > ${2}/variant_calling/qiaseq/${5}_smcounter.vcf
    vcfuniq ${2}/variant_calling/qiaseq/${5}_smcounter.vcf > ${2}/variant_calling/qiaseq/${5}_smcounter_uniq.vcf
    bgzip ${2}/variant_calling/qiaseq/${5}_smcounter_uniq.vcf && tabix -f -p vcf ${2}/variant_calling/qiaseq/${5}_smcounter_uniq.vcf.gz

    vcf-sort ${2}/variant_calling/qiaseq/${10}_smcounter_tmp.vcf | vcfuniq > ${2}/variant_calling/qiaseq/${10}_smcounter.vcf
    vcfuniq ${2}/variant_calling/qiaseq/${10}_smcounter.vcf > ${2}/variant_calling/qiaseq/${10}_smcounter_uniq.vcf
    bgzip ${2}/variant_calling/qiaseq/${10}_smcounter_uniq.vcf && tabix -f -p vcf ${2}/variant_calling/qiaseq/${10}_smcounter_uniq.vcf.gz
    
    echo "python /usr/local/bin/inc_info_caller.py ${2}/variant_calling/qiaseq/${5}_smcounter_uniq.vcf.gz ${2}/variant_calling/qiaseq/${5}_smCounter.vcf smCounter"
    

    echo "vcf-isec -f -n +1 ${2}/variant_calling/qiaseq/${5}_smcounter_uniq.vcf.gz ${2}/variant_calling/${5}_lofreq_anno.vcf.gz ${2}/variant_calling/${5}_mutect2.vcf.gz | bgzip > ${2}/variant_calling/${5}_tumor.vcf.gz"
    vcf-isec -f -n +1 ${2}/variant_calling/qiaseq/${5}_smcounter_uniq.vcf.gz ${2}/variant_calling/${5}_lofreq_anno.vcf.gz ${2}/variant_calling/${5}_mutect2.vcf.gz | bgzip > ${2}/variant_calling/${5}_tumor.vcf.gz
    echo "vcf-isec -f -n +1 ${2}/variant_calling/qiaseq/${10}_smcounter_uniq.vcf.gz ${2}/variant_calling/${10}_lofreq_anno.vcf.gz ${2}/variant_calling/${10}_mutect2.vcf.gz | bgzip > ${2}/variant_calling/${5}_normal.vcf.gz"
    vcf-isec -f -n +1 ${2}/variant_calling/qiaseq/${10}_smcounter_uniq.vcf.gz ${2}/variant_calling/${10}_lofreq_anno.vcf.gz ${2}/variant_calling/${10}_mutect2.vcf.gz | bgzip > ${2}/variant_calling/${5}_normal.vcf.gz

    tabix -p vcf -f ${2}/variant_calling/${5}_tumor.vcf.gz
    tabix -p vcf -f ${2}/variant_calling/${5}_normal.vcf.gz

    echo "vcf-merge -c none ${2}/variant_calling/${5}_tumor.vcf.gz ${2}/variant_calling/${5}_normal.vcf.gz | bgzip > ${2}/variant_calling/${5}_.vcf.gz"
    vcf-merge -c none ${2}/variant_calling/${5}_tumor.vcf.gz ${2}/variant_calling/${5}_normal.vcf.gz | bgzip > ${2}/variant_calling/${5}_.vcf.gz
    python /usr/local/bin/inc_parse_vcf_header.py ${2}/variant_calling/${5}_.vcf.gz ${2}/variant_calling/${5}.vcf
    bgzip -f ${2}/variant_calling/${5}.vcf
    tabix -p vcf -f ${2}/variant_calling/${5}.vcf.gz
    date

}


function statistics(){
    echo "Running statistics for sample ${5}"
    echo "picardv2.18.6.sh CollectAlignmentSummaryMetrics I=${3}/mapping/${5}_sorted.bam O=${3}/mapping/${sample}_alignment_metrics R=${8}"
    picardv2.18.6.sh CollectAlignmentSummaryMetrics I=${3}/mapping/${5}_sorted.bam O=${3}/mapping/${sample}_alignment_metrics R=${8}
    mapping_reads=`head ${3}/mapping/${sample}_alignment_metrics | tail -1 | awk '{print $6}'`
    echo "picardv2.18.6.sh CollectAlignmentSummaryMetrics I=${3}/mapping/ontarget_${5}_sorted.bam O=${3}/mapping/ontarget_${sample}_alignment_metrics R=${8}"
    picardv2.18.6.sh CollectAlignmentSummaryMetrics I=${3}/mapping/ontarget_${5}_sorted.bam O=${3}/mapping/ontarget_${sample}_alignment_metrics R=${8}
    ontarget_mapping_reads=`head ${3}/mapping/ontarget_${sample}_alignment_metrics | tail -1 | awk '{print $6}'`
    echo "picardv2.18.6.sh CollectAlignmentSummaryMetrics I=${3}/mapping/ontarget_${sample}_processed.bam O=${3}/mapping/ontarget_${sample}_processed_alignment_metrics R=${8}"
    picardv2.18.6.sh CollectAlignmentSummaryMetrics I=${3}/mapping/ontarget_${sample}_processed.bam O=${3}/mapping/ontarget_${sample}_processed_alignment_metrics R=${8}
    dedup_reads=`head ${3}/mapping/ontarget_${sample}_processed_alignment_metrics | tail -1 | awk '{print $6}'`
    rawdata=`zcat ${2} | wc -l | awk '{print $1/2}'`
    tmp_duplicates=`echo "${ontarget_mapping_reads} - ${dedup_reads}" | bc`    
    duplicates=`echo "scale=4; ${tmp_duplicates}/${rawdata}" | bc | awk '{print $1*100}'`    
    all_mapped=`echo "scale=4; ${mapping_reads}/${rawdata}" | bc | awk '{print $1*100}'`
    prinseq_reads=`wc -l ${3}/QC/filtered_${sample}_1.fastq | awk '{print $1/2}'`
    prinseq=`echo "scale=4; ${prinseq_reads}/${rawdata}" | bc | awk '{print 100-($1*100)}'`
    bam_target_reads=`echo "scale=4; ${ontarget_mapping_reads}/${rawdata}" | bc | awk '{print $1*100}'`
    bam_target_reads=`head ${3}/mapping/ontarget_${sample}_alignment_metrics | tail -1 | awk '{print $6}'`
    usable_reads=`echo "scale=4; ${dedup_reads}/${rawdata}" | bc | awk '{print $1*100}'`
    specificity=`echo "scale=4; ${bam_target_reads}/${mapping_reads}" | bc | awk '{print $1*100}'`
    echo "picardv2.18.6.sh CollectInsertSizeMetrics I=${3}/mapping/ontarget_${5}_processed.bam O=${3}/QC/${5}_insert H=${3}/QC/${5}_hist"
    picardv2.18.6.sh CollectInsertSizeMetrics I=${3}/mapping/ontarget_${5}_processed.bam O=${3}/QC/${5}_insert H=${3}/QC/${5}_hist
    bam_all_reads=`samtools flagstat ${3}/mapping/${5}_sorted.bam | head -5 | tail -1 | awk '{print $1}'`
    all_vars_lf=`zcat ${3}/variant_calling/${5}_lofreq_anno.vcf.gz | grep -v "#" -c`
    all_vars_smc=`zcat ${3}/variant_calling/qiaseq/${5}_smcounter_uniq.vcf.gz | grep -v "#" -c`
    all_vars_mu=`zcat ${3}/variant_calling/${5}_mutect2.vcf.gz | grep -v "#" -c`
    all_vars=`zcat ${3}/variant_calling/${5}.vcf.gz | grep -v "#" -c`
    mean_insert_tmp=`head -8 ${3}/QC/${5}_insert | tail -1 | awk '{print $6}'`
    mean_insert=`echo "scale=2; ${mean_insert_tmp}/1" | bc`
    deviation_tmp=`head -8 ${3}/QC/${5}_insert | tail -1 | awk '{print $7}'`
    deviation=`echo "scale=2; ${deviation_tmp}/1" | bc`
    design=`awk '{print $3-$2+1}' ${3}/mapping/target.bed | suma`
    echo "bedtools subtract -a ${3}/mapping/target.bed -b ${4}/QC/mapping/${5}_covered_regions_${10}x.bed > ${3}/QC/${5}_non_covered_${10}x.bed"
    bedtools subtract -a ${7}/target.bed -b ${4}/QC/mapping/${5}_covered_regions_${10}x.bed > ${3}/QC/${5}_non_covered_${10}x.bed
    echo "bedtools intersect -a ${3}/QC/${5}_non_covered_${10}x.bed -b ${7}/inc_db.vcf -wao | awk '{if ($5 ~ /chr/) print $12"\t"$1"\t"$2"\t"$3"\t"$4; else print "--\t"$1"\t"$2"\t"$3"\t"$4}' | sed 's/.*CDS=//' | sed 's/;CNT.*\tchr/\tchr/' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' > ${4}/QC/mapping/${5}_non_covered_regions_${10}x.bed"
    bedtools intersect -a ${3}/QC/${5}_non_covered_${10}x.bed -b ${7}/inc_db.vcf -wao | awk '{if ($5 ~ /chr/) print $12"\t"$1"\t"$2"\t"$3"\t"$4; else print "--\t"$1"\t"$2"\t"$3"\t"$4}' | sed 's/.*CDS=//' | sed 's/;CNT.*\tchr/\tchr/' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' > ${4}/QC/mapping/${5}_non_covered_regions_${10}x.bed
    bgzip -f ${3}/QC/mapping/${5}_covered_regions_${10}x.bed
    tabix -p bed -f ${3}/QC/mapping/${5}_covered_regions_${10}x.bed.gz
    sample200=`awk '{if ($3>=250) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | bedtools merge | awk '{print $3-$2+1}' | suma`
    depth200=`echo "scale=4; ${sample200}/${design}" | bc | awk '{print $1*100}'`
    sample100=`awk '{if ($3>=100) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | bedtools merge | awk '{print $3-$2+1}' | suma`
    depth100=`echo "scale=4; ${sample100}/${design}" | bc | awk '{print $1*100}'`
    sample500=`awk '{if ($3>=500) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | bedtools merge | awk '{print $3-$2+1}' | suma`
    depth500=`echo "scale=4; ${sample500}/${design}" | bc | awk '{print $1*100}'`
    inc_calculate_IQR_and_coverage_plot.R ${3}/mapping/${5}_depth ${3}/QC/${5}_R_coverage_summary ${4}/QC/mapping/${5}_coverage_plot.pdf ${5}
    first_quantile=`tail -1 ${3}/QC/${5}_R_coverage_summary | awk '{print $2}'`
    third_quantile=`tail -1 ${3}/QC/${5}_R_coverage_summary | awk '{print $5}'`
    median_coverage=`tail -1 ${3}/QC/${5}_R_coverage_summary | awk '{print $3}'`
    echo "${5},${rawdata},${prinseq},${all_mapped},${duplicates},${usable_reads},${specificity},${first_quantile},${third_quantile},${median_coverage},${depth100},${depth200},${depth500},${mean_insert},${deviation},${all_vars_smc},${all_vars_lf},${all_vars_mu},${all_vars}" >> ${4}/QC/mapping/${6}_global_mapping_stats.csv
    date
}


function check_list_of_transcripts(){
    while IFS='' read -r line || [[ -n "${line}" ]]; do
        transcript=`echo ${line} | awk '{print $3}'`
        gene=`echo ${line} | awk '{print $1}'`
        has_final_info=`grep $transcript ${3}`
        has_intermediate_info=`grep $transcript ${2}`
        if [  -z "${has_final_info}" ]
        then
            if [ ! -z "${has_intermediate_info}" ]
            then
                echo "-------------ERROR----SOMETHING MAY HAPPENS WITH THE LIST OF TRANSCRIPTS-----------------"
                echo "----------------CHECK THAT NO RESULTS ARE OBTAINED FOR ${line}-----------------------"
            fi
        fi
    done < "${1}"
}

function annotate_variants(){
    echo "vep -i ${1} --use_transcript_ref --force_overwrite --cache --merged --dir ${5}/hg38/ --dir_plugins ${5}/hg38/Plugins --fasta ${4} --sift b --polyphen b --gene_phenotype --numbers --vcf_info_field ANN --terms SO --hgvs --shift_hgvs 1 --canonical --biotype --xref_refseq --max_af --af_esp --af_gnomad --pubmed --minimal --force_overwrite --vcf --af --symbol --domains --plugin Conservation,GERP_CONSERVATION_SCORE,mammals --fork ${7} --buffer_size 200 -o ${2}/variant_calling/${6}_annotated_vep.vcf"
    vep -i ${1} --use_transcript_ref --force_overwrite --cache --merged --dir ${5}/hg38/ --dir_plugins ${5}/hg38/Plugins --fasta ${4} --sift b --polyphen b --gene_phenotype --numbers --vcf_info_field ANN --terms SO --hgvs --shift_hgvs 1 --canonical --biotype --xref_refseq --max_af --af_esp --af_gnomad --pubmed --minimal --force_overwrite --vcf --af --symbol --domains --plugin Conservation,database,GERP_CONSERVATION_SCORE,mammals --fork ${7} --buffer_size 200 -o ${2}/variant_calling/${6}_annotated_vep.vcf
    bgzip -f ${2}/variant_calling/${6}_annotated_vep.vcf
    tabix -p vcf -f ${2}/variant_calling/${6}_annotated_vep.vcf.gz
    gunzip -c ${2}/variant_calling/${6}_annotated_vep.vcf.gz > ${2}/variant_calling/${6}_annotated_vep.vcf
    date
    #python /usr/local/bin/inc_parse_vcf_header.py ${2}/variant_calling/${6}_annotated_vep.vcf.gz ${2}/variant_calling/${6}_annotated_vep_tmp.vcf
    echo "bcftools annotate -a ${5}/cancerhospots.bed.gz -c CHROM,FROM,TO,HS -h <(echo '##INFO=<ID=HS,Number=1,Type=String,Description="Tag variants in cancer hotspots cancerhotspot.org">') ${2}/variant_calling/${6}_annotated_vep.vcf > ${2}/variant_calling/${6}_annotated_uniq_h.vcf"
    bcftools annotate -a ${5}/cancerhospots.bed.gz -c CHROM,FROM,TO,HS -h <(echo '##INFO=<ID=HS,Number=1,Type=String,Description="Tag variants in cancer hotspots cancerhotspot.org">') ${2}/variant_calling/${6}_annotated_vep.vcf.gz > ${2}/variant_calling/${6}_annotated_uniq_h.vcf
    echo "snpSift_v4.3t.sh annotate ${5}/AllMutations_COSMIC_v85.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h.vcf > ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf"
    snpSift_v4.3t.sh annotate ${5}/AllMutations_COSMIC_v85.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h.vcf > ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf
    echo "inc_mv_cosmic_id_to_annotation.py -i ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf -o ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf"
    inc_mv_cosmic_id_to_annotation.py -i ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf -o ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf
    echo "snpSift_v4.3t.sh annotate ${5}/inc_db.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf > ${3}/variant_calling/${6}_annotated.vcf"
    snpSift_v4.3t.sh annotate ${5}/inc_db.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf > ${3}/variant_calling/${6}_annotated.vcf
    bgzip -f ${3}/variant_calling/${6}_annotated.vcf
    echo "inc_vcf_to_csv.py -f ${3}/variant_calling/${6}_annotated.vcf.gz -o ${3}/variant_calling/${6}_annotated.csv -m 0.02 -n 1 -r ${5}/list_of_transcripts"
    date
    #python /usr/local/bin/inc_genotype_strelka2.py ${3}/variant_calling/${6}_annotated.vcf.gz ${3}/variant_calling/${6}_annotated2.vcf
    #gzip ${3}/variant_calling/${6}_annotated2.vcf
    #python /usr/local/bin/inc_parse_strelka2_vcf.py ${3}/variant_calling/${6}_annotated.vcf.gz ${3}/variant_calling/${6}_annotated2.vcf
    python /usr/local/bin/inc_vcf_to_csv_lofreq.py -f ${3}/variant_calling/${6}_annotated.vcf.gz -o ${3}/variant_calling/${6}_annotated.csv -m 0.000000000000000000000000002 -n 1 -r ${5}/list_of_transcripts
    python /usr/local/bin/inc_aa_3to1.py ${3}/variant_calling/${6}_annotated.csv
    echo "inc_samples_contamination.R ${3} 0.7"
    inc_samples_contamination.R ${3} 0.7
    echo "inc_non_covered_regions.py ${3}/QC/mapping/${6}_global_mapping_stats.csv ${3} ${2} ${5}/target.bed ${8}"
    python /usr/local/bin/inc_non_covered_regions.py ${3}/QC/mapping/${6}_global_mapping_stats.csv ${3} ${2} ${5}/target.bed ${8}
    echo "inc_gene_cover.py ${3} ${2} ${5} ${8}"
    python /usr/local/bin/inc_gene_cover.py ${3} ${2} ${5} ${8}
    echo "inc_csv_format.R ${3}/variant_calling/${6}_annotated.csv ${3}/variant_calling"
    inc_csv_format.R ${3}/variant_calling/${6}_annotated.csv ${3}/variant_calling
    date
    check_list_of_transcripts ${5}/list_of_transcripts ${2}/variant_calling/${6}_annotated_vep.vcf ${3}/variant_calling/${6}_annotated.csv ${2}/variant_calling
    date
    check_list_of_transcripts ${5}/list_of_transcripts ${2}/variant_calling/${6}_annotated_vep.vcf ${3}/variant_calling/${6}_annotated.csv ${2}/variant_calling
}

function usage() {
    echo "Usage: $0"
    echo "This script runs the pipeline for target seq data"
    echo "Mandatory parameters:"
    echo "-d directory where rawdata is stored"
    echo "-s comma-separated list of samples"
    echo "-u if libraries were built with UMI technology"
    echo "-a the design is amplicon-based"
    echo "-n output file name"
    echo "-t temporal working directory"
    echo "-o output directory"
    echo "-r number of threads"
    echo "-g reference genome"
    echo "-b directory with design file in BED format and other files needed by GATK"
    echo "-c depth of coverage"
    echo "-f minimun value to consider a region covered"
    1>&2; exit 1;
}

while getopts "d:s:n:t:o:r:g:b:c:u:a:z:" opt; do
    case ${opt} in
        d)
                d=${OPTARG} ;;
        s)
                s=${OPTARG} ;;
        n)
                n=${OPTARG} ;;
        t)
                t=${OPTARG} ;;
        o)
                o=${OPTARG} ;;
        r)
                r=${OPTARG} ;;
        g)
                g=${OPTARG} ;;
        b)
                b=${OPTARG} ;;
        c)
                c=${OPTARG} ;;
        z)
                z=${OPTARG} ;;
        u)
                u="yes" ;;
        a)
                a="no";;
        *)
                usage ;;
    esac
done

if [ -z "$d" ] || [ -z "$s" ] || [ -z "$n" ] || [ -z "$c" ] || [ -z "$t" ] || [ -z "$o" ] || [ -z "$r" ] || [ -z "$g" ] || [ -z "$b" ] || [ -z "$z" ]
then
        echo "ERROR: -d -s -n -t -o -r -g -c -b and -m are mandatory parameters"
        usage
        exit
fi

echo "Your command: "$@

checkDirectories ${t}/mapping
checkDirectories ${t}/variant_calling
checkDirectories ${t}/QC
checkDirectories ${o}/mapping
checkDirectories ${o}/variant_calling
checkDirectories ${o}/QC/mapping
checkDirectories ${o}/QC/rawdata
ln -s -f /nfs/home/databases/INCLIVA_somatic/20190805/inc_db.vcf.gz ${b}/
ln -s -f /nfs/home/databases/INCLIVA_somatic/20190805/inc_db.vcf.gz.tbi ${b}/
ln -s -f /nfs/home/databases/ensembl/v96/hg38/ ${b}/
ln -s -f /nfs/home/databases/dbSNP/151/All_20180418.vcf.gz ${b}/
ln -s -f /nfs/home/databases/dbSNP/151/All_20180418.vcf.gz.tbi ${b}/
ln -s -f /nfs/home/databases/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ${b}/
ln -s -f /nfs/home/databases/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi ${b}/
ln -s -f /nfs/home/databases/cosmic/cosmicv85/AllMutations_COSMIC_v85.vcf.gz ${b}/
ln -s -f /nfs/home/databases/cosmic/cosmicv85/AllMutations_COSMIC_v85.vcf.gz.tbi ${b}/
ln -s -f /nfs/home/databases/cancer_hotspots/20180723/cancerhospots.bed.gz ${b}/
ln -s -f /nfs/home/databases/cancer_hotspots/20180723/cancerhospots.bed.gz.tbi ${b}/
ln -s -f /nfs/home/databases/INCLIVA_somatic/20190805/inc_db.vcf ${b}/
bedtools intersect -a /nfs/home/databases/GATK/small_exac_common_3.hg38.vcf.gz -b ${b}/capture.bed -header > ${b}/exac_roche.vcf
bgzip -f ${b}/exac_roche.vcf && tabix -p vcf -f ${b}/exac_roche.vcf.gz


echo "Sample,Rawdata,%LowQReads,%MappedReads,%DuplicateReads,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,Nt_100x,Nt_250x,Nt_500x,Mean_insert,Insert_SD,NumVars(smCounter),NumVars(loFreq),NumIndels(Mutect2),NumVars(Merged)" > ${o}/QC/mapping/${n}_global_mapping_stats.csv
sample_list=`echo ${s} | sed 's/,/ /gi'`
normal_list=`echo ${z} | sed 's/,/ /gi'`
cwd=`pwd`
read -a arr <<< $sample_list

for i in $(seq 1 ${#arr[@]});
do
    
    echo $i
    sample=`echo $sample_list | cut -d" " -f"$i"`
    normal=`echo $normal_list | cut -d" " -f"$i"`
    
    echo "Running analysis pipeline on sample ${sample}"
    R1=`ls ${d}/${sample}_R1.fast*`
    R2=`ls ${d}/${sample}_R2.fast*`
    if [ -z "$u" ]
    then
        echo "------IS NO UMI------"
        has_umi="no"
    else
        echo "------HAS UMI------"
        has_umi="yes"
        umi=`ls ${d}/$sample*_I2*.fast*`
    fi
    if [ -z "$a" ]
    then
        echo "------NO AMPLICONS------"
        amplicon="no"
    else
        echo "------DESIGN IS AMPLICON-BASED------"
        amplicon="yes"
    fi
    get_software_version ${o}/QC/software_versions ${has_umi} ${o}/QC/hardware_info ${t}
    fastq_preprocessing ${R1} ${R2} ${sample} ${t} ${o} ${r} ${has_umi} ${umi}
    mapping_and_bampostprocessing ${t}/QC/filtered_${sample}_1.fastq ${t}/QC/filtered_${sample}_2.fastq ${t} ${o} ${r} ${sample} ${g} ${b} ${has_umi} ${amplicon}
    variant_calling ${o}/mapping/${sample}.bam ${t} ${o} ${r} ${sample} ${g} ${b} ${c} ${o}/mapping/${normal}.bam ${normal}
    variant_calling2 ${o}/mapping/${normal}.bam ${t} ${o} ${r} ${sample} ${g} ${b} ${c} ${o}/mapping/${normal}.bam ${normal}
    variant_calling3 ${o}/mapping/${sample}.bam ${t} ${o} ${r} ${sample} ${g} ${b} ${c} ${o}/mapping/${normal}.bam ${normal}
    qiaseq ${o}/mapping/${sample}.bam ${t} ${o} ${r} ${sample} ${g} ${b} ${c} ${o}/mapping/${normal}.bam ${normal}
    statistics ${o}/mapping/${sample}.bam ${R1} ${t} ${o} ${sample} ${n} ${b} ${g} ${has_umi} ${c} ${amplicon}
    list_files+=" ${t}/variant_calling/${sample}.vcf.gz "


done

echo "Combining files"
echo "vcf-merge -c none $list_files | bgzip > ${o}/variant_calling/${n}.vcf.gz"
vcf-merge -c none $list_files | bgzip > ${t}/variant_calling/${n}.vcf.gz
date
tabix -p vcf -f ${o}/variant_calling/${n}.vcf.gz

echo "annotate_variants ${o}/variant_calling/${n}.vcf.gz ${t} ${o} ${g} ${b} ${n} ${r} ${c}"
annotate_variants ${t}/variant_calling/${n}.vcf.gz ${t} ${o} ${g} ${b} ${n} ${r} ${c}
