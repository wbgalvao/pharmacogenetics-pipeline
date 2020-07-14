#!/usr/bin/env nextflow


/* Command line arguments:
    --samples: Path to directory where the target batch fastq files are stored
    --reference: Path to GRCh37 reference genome (with index and .dict files) directory
    --resources: Path to Broad's Institute b37 resource bundle
    --list: Path to .list file with target intervals (stargazer-target-genes.sorted.list)
    --results: Target directory to publish pipeline outputs
*/


sampleReadsFilesChannel = Channel.fromFilePairs("${params.samples}/*_R{1,2}*.fastq.gz")
referenceGenomeDirectoryChannel = Channel.value(params.reference)
resourceBundleDirectoryChannel = Channel.value(params.resources)
intervalsFilePathChannel = Channel.value(params.list)
resultsDirectory = params.results


targetGenesChannel = Channel.fromList([
    "CACNA1S", "CFTR", "CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2A7", "CYP2A13",
    "CYP2B6", "CYP2B7", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2D7", "CYP2E1",
    "CYP2F1", "CYP2J2", "CYP2R1", "CYP2S1", "CYP2W1", "CYP3A4", "CYP3A5", "CYP3A7",
    "CYP3A43", "CYP4B1", "CYP26A1", "CYP4F2", "CYP19A1", "DPYD", "G6PD", "GSTM1",
    "GSTP1", "GSTT1", "IFNL3", "NAT1", "NAT2", "NUDT15", "POR", "RYR1", "SLC15A2",
    "SLC22A2", "SLCO1B1", "SLCO1B3", "SLCO2B1", "SULT1A1", "TBXAS1", "TPMT", "UGT1A1",
    "UGT1A4", "UGT2B7", "UGT2B15", "UGT2B17", "VKORC1"
])


batchName = file(params.samples).getBaseName()
pipelineOutputPath = "${resultsDirectory}/${batchName}"

// These channels are used later to organize and publish pipeline artifacts per sample
sampleReadsFilesChannel.into {
    samplesFastqsChannel;
    samplesVcfChannel;
    samplesHaplotypesChannel;
    samplesCNVReportsChannel;
}


process alignReadFiles {

    container "bwa:v0.7.17"

    input:
    tuple val(sample), file(fastqs) from samplesFastqsChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel

    output:
    tuple val(sample), file("${sample}.sam") into samChannel

    """
    bwa mem \
        -H \
        -a \
        -d \
        -S \
        -R "@RG\\tID:${fastqs[0]}\\tPL:ILLUMINA\\tSM:${sample}" \
        ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
        ${fastqs} > ${sample}.sam
    """

}


process compressAlignmentFiles {

    container "samtools:1.10"

    input:
    tuple val(sample), file(sam) from samChannel

    output:
    tuple val(sample), file("${sample}.bam") into bamChannel

    """
    samtools view -b ${sam} > ${sample}.bam
    """

}


process sortAlignments {

    container "samtools:1.10"

    input:
    tuple val(sample), file(bam) from bamChannel

    output:
    tuple val(sample), file("${sample}.sorted.bam") into sortedBamChannel

    """
    samtools sort ${bam} > ${sample}.sorted.bam
    """

}


process markDuplicates {

    container "broadinstitute/picard:2.22.8"

    input:
    tuple val(sample), file(sortedBam) from sortedBamChannel

    output:
    tuple val(sample),
        file("${sample}.sorted.duplicate_marked.bam"),
        file("${sample}.sorted.duplicate_marked.bai") into duplicateMarksChannel

    """
    java -jar /usr/picard/picard.jar \
        MarkDuplicates \
        INPUT="${sortedBam}" \
        OUTPUT=${sample}.sorted.duplicate_marked.bam \
        METRICS_FILE=${sample}.duplicate_metrics.tsv \
        VALIDATION_STRINGENCY=SILENT \
        CLEAR_DT="false" \
        ADD_PG_TAG_TO_READS=false \
        && java -jar /usr/picard/picard.jar \
            BuildBamIndex \
            INPUT=${sample}.sorted.duplicate_marked.bam \
            OUTPUT=${sample}.sorted.duplicate_marked.bai
    """

}


duplicateMarksChannel.into { 
    baseRecalibratorChannel;
    printReadsChannel;
}


process createRecalibrationData {

    container "broadinstitute/gatk3:3.8-0"

    input:
    tuple val(sample),
        file(sortedDuplicateMarkedBam),
        file(sortedDuplicateMarkedBamIndex) from baseRecalibratorChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path resourceBundleDirectory from resourceBundleDirectoryChannel

    output:
    file "${sample}.recal_data.table" into recalibrationDataChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type BaseRecalibrator \
        --input_file ${sortedDuplicateMarkedBam} \
        --out ${sample}.recal_data.table \
        --reference_sequence ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY \
        --useOriginalQualities \
        --knownSites ${resourceBundleDirectory}/${params.DBSNP_VCF} \
        --knownSites ${resourceBundleDirectory}/${params.GOLD_STANDART_INDELS_1000G_VCF} \
        --knownSites ${resourceBundleDirectory}/${params.OMNI25_1000G_VCF} \
        --knownSites ${resourceBundleDirectory}/${params.PHASE1_INDELS_1000G_VCF}
    """

}


process applyRecalibration {

    container "broadinstitute/gatk3:3.8-0"
    publishDir "${pipelineOutputPath}/${sample}", mode: "copy"

    input:
    tuple val(sample),
        file(sortedDuplicateMarkedBam),
        file(sortedDuplicateMarkedBamIndex) from printReadsChannel
    file recalibrationDataTable from recalibrationDataChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel

    output:
    tuple val(sample),
        file("${sample}.sorted.duplicate_marked.recalibrated.bam"),
        file("${sample}.sorted.duplicate_marked.recalibrated.bai") into recalibrationChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type PrintReads \
        --input_file ${sortedDuplicateMarkedBam} \
        --out ${sample}.sorted.duplicate_marked.recalibrated.bam \
        --reference_sequence ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
        --BQSR ${recalibrationDataTable} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY \
        --useOriginalQualities
    """

}


recalibrationChannel.into {
    haplotypeCallerChannel;
    depthOfCoverageChannel;
}


process callGerminativeVariants {

    container "broadinstitute/gatk3:3.8-0"

    input:
    tuple val(sample),
        file(sortedDuplicateMarkedRecalibratedBam),
        file(sortedDuplicateMarkedRecalibratedBamIndex) from haplotypeCallerChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path intervalsFilePath from intervalsFilePathChannel
    path resourceBundleDirectory from resourceBundleDirectoryChannel

    output:
    tuple file("${sample}.g.vcf"),
        file("${sample}.g.vcf.idx") into gvcfChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type HaplotypeCaller \
        --input_file ${sortedDuplicateMarkedRecalibratedBam} \
        --out ${sample}.g.vcf \
        --reference_sequence ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
        --intervals ${intervalsFilePath} \
        --emitRefConfidence GVCF \
        --dbsnp ${resourceBundleDirectory}/${params.DBSNP_VCF} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process combineGVCFs {

    container "broadinstitute/gatk3:3.8-0"

    input:
    file gvcf from gvcfChannel.collect()
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path intervalsFilePath from intervalsFilePathChannel
    path resourceBundleDirectory from resourceBundleDirectoryChannel

    output:
    tuple file("${batchName}.g.vcf"),
        file("${batchName}.g.vcf.idx") into combinedGVCFChannel

    """
    ls *.g.vcf > gvcfs.list
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type CombineGVCFs \
        --variant gvcfs.list \
        --out "${batchName}.g.vcf" \
        --reference_sequence ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
        --intervals ${intervalsFilePath} \
        --dbsnp ${resourceBundleDirectory}/${params.DBSNP_VCF} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process genotypeGVCF {

    container "broadinstitute/gatk3:3.8-0"

    input:
    file gvcf from combinedGVCFChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path intervalsFilePath from intervalsFilePathChannel
    path resourceBundleDirectory from resourceBundleDirectoryChannel

    output:
    file "${batchName}.vcf" into vcfChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type GenotypeGVCFs \
        --variant ${gvcf[0]} \
        --out "${batchName}.vcf" \
        --reference_sequence ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
        --intervals ${intervalsFilePath} \
        --dbsnp ${resourceBundleDirectory}/${params.DBSNP_VCF} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process filterVariants {

    container "broadinstitute/gatk3:3.8-0"

    input:
    file vcf from vcfChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path intervalsFilePath from intervalsFilePathChannel

    output:
    file "${batchName}.filtered.vcf" into filteredVcfChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
    --analysis_type VariantFiltration \
    --variant ${vcf} \
    --out "${batchName}.filtered.vcf" \
    --reference_sequence ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
    --intervals ${intervalsFilePath} \
    --filterExpression "QUAL <= 50.0" \
    --filterName QUALFilter \
    --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


filteredVcfChannel.into {
    stargazerVcfChannel;
    selectVariantsVcfChannel;
}


process assessAlignmentCoverage {

    container "broadinstitute/gatk3:3.8-0"

    input:
    file resultAlignment from depthOfCoverageChannel.collect()
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path intervalsFilePath from intervalsFilePathChannel

    output:
    file "${batchName}.table" into depthOfCoverageResultChannel

    """
    ls *.bam > bam.list
    java -jar /usr/GenomeAnalysisTK.jar \
    --analysis_type DepthOfCoverage \
    --input_file bam.list \
    --out ${batchName}.table \
    --reference_sequence ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
    --intervals ${intervalsFilePath} \
    --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY \
    --minMappingQuality 1 \
    --omitIntervalStatistics \
    --omitPerSampleStats \
    --omitLocusTable
    """

}


process callHaplotypes {

    container "stargazer:v1.0.8"
    errorStrategy "ignore" // ¯\_(ツ)_/¯

    input:
    file vcf from stargazerVcfChannel
    file gdf from depthOfCoverageResultChannel
    each gene from targetGenesChannel

    output:
    tuple file("${batchName}.${gene}.stargazer-genotype.txt"),
        file("${batchName}.${gene}.stargazer-genotype.log"),
        path("${batchName}.${gene}.stargazer-genotype.project") into stargazerResultsChannel

    """
    mkdir /${batchName}/
    python /usr/Stargazer_v1.0.8/stargazer.py genotype \
    --target_gene ${gene} \
    --control_gene EGFR \
    --data_type ts \
    --gdf ${gdf} \
    --vcf ${vcf} \
    --output_dir . \
    --output_prefix ${batchName}.${gene} \
    """

}


/* The haplotype definition pipeline ends here.
   From now on, all tasks are to publish specific process outputs into
   different sample directories.
   The objectives here are:
    1. Split the result multi sample vcf into different single sample vcfs
    2. Gather all the '.stargazer-genotype.txt' files, group all genes results
       and create and publish a new result tsv per sample.
    3. Apply the same methodoly aforementioned for the cnv reports created by Stargazer
*/


stargazerResultsChannel.into {
    stargazerHaplotypesChannel;
    stargazerCnvReportsChannel;
}


process splitVCFPerSample {

    container "broadinstitute/gatk3:3.8-0"
    publishDir "${pipelineOutputPath}/${sample}", mode: "copy"

    input:
    tuple val(sample), file(_) from samplesVcfChannel
    file vcf from selectVariantsVcfChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel

    output:
    tuple val(sample), file("${sample}.filtered.vcf"), file("${sample}.filtered.vcf.idx")

    """
    java -jar /usr/GenomeAnalysisTK.jar \
    --analysis_type SelectVariants \
    --variant ${vcf} \
    --out ${sample}.filtered.vcf \
    --reference_sequence ${referenceGenomeDirectory}/${params.REFERENCE_GENOME_FASTA} \
    --sample_name ${sample} \
    --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process gatherStargazerResultsPerSample {

    container "pandas:1.0.5"
    publishDir "${pipelineOutputPath}/${sample}", mode: "copy"

    input:
    tuple val(sample), file(_) from samplesHaplotypesChannel
    file stargazerResults from stargazerHaplotypesChannel.collect()

    output:
    file "${sample}.haplotypes.tsv"

    """
    python ${params.PIPELINE_FOLDER}/scripts/merge_stargazer_output_per_sample.py \
        ${sample} \
        ${params.PIPELINE_FOLDER}/intervals/list.json
    """

}


process createSampleCNVReport {

    container "poppler:0.82.0-r1"
    publishDir "${pipelineOutputPath}/${sample}", mode: "copy"

    input:
    tuple val(sample), file(_) from samplesCNVReportsChannel
    path stargazerResults from stargazerCnvReportsChannel.collect()

    output:
    file "${sample}.CNV-report.pdf"

    """
    pdfunite \$(find -L *.project -name ${sample}.pdf) ${sample}.CNV-report.pdf
    """

}
