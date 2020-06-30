#!/usr/bin/env nextflow


sampleReadsFilesChannel = Channel.fromFilePairs("${params.samples}/*_R{1,2}*.fastq.gz")
referenceGenomeDirectoryChannel = Channel.value(params.reference)
resourceBundleDirectoryChannel = Channel.value(params.resources)
intervalsFilePathChannel = Channel.value(params.list)
resultsDirectory = params.results

batchName = file(params.samples).getBaseName()
pipelineOutputPath = "${resultsDirectory}/${batchName}"

REFERENCE_GENOME_FASTA = "Homo_sapiens.GRCh37.dna.primary_assembly.fa"
DBSNP_VCF = "dbsnp_138.b37.vcf"
GOLD_STANDART_INDELS_1000G_VCF = "Mills_and_1000G_gold_standard.indels.b37.vcf"
OMNI25_1000G_VCF = "1000G_omni2.5.b37.vcf"
PHASE1_INDELS_1000G_VCF = "1000G_phase1.indels.b37.vcf"


process alignReadFiles {

    container "gcr.io/pharmacogenetics/bwa:v0.7.17"

    input:
    tuple val(sample), file(fastqs) from sampleReadsFilesChannel
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
        ${referenceGenomeDirectory}/${REFERENCE_GENOME_FASTA} \
        ${fastqs} > ${sample}.sam
    """

}


process compressAlignmentFiles {

    container "gcr.io/pharmacogenetics/samtools:1.10"

    input:
    tuple val(sample), file(sam) from samChannel

    output:
    tuple val(sample), file("${sample}.bam") into bamChannel

    """
    samtools view -b ${sam} > ${sample}.bam
    """

}


process sortAlignments {

    container "gcr.io/pharmacogenetics/samtools:1.10"

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
        --reference_sequence ${referenceGenomeDirectory}/${REFERENCE_GENOME_FASTA} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY \
        --useOriginalQualities \
        --knownSites ${resourceBundleDirectory}/${DBSNP_VCF} \
        --knownSites ${resourceBundleDirectory}/${GOLD_STANDART_INDELS_1000G_VCF} \
        --knownSites ${resourceBundleDirectory}/${OMNI25_1000G_VCF} \
        --knownSites ${resourceBundleDirectory}/${PHASE1_INDELS_1000G_VCF}
    """

}


process applyRecalibration {

    container "broadinstitute/gatk3:3.8-0"
    publishDir "${pipelineOutputPath}/${sample}"

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
        --reference_sequence ${referenceGenomeDirectory}/${REFERENCE_GENOME_FASTA} \
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
    publishDir "${pipelineOutputPath}/${sample}"

    input:
    tuple val(sample),
        file(sortedDuplicateMarkedRecalibratedBam),
        file(sortedDuplicateMarkedRecalibratedBamIndex) from haplotypeCallerChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path intervalsFilePath from intervalsFilePathChannel
    path resourceBundleDirectory from resourceBundleDirectoryChannel

    output:
    tuple file("${sample}.g.vcf"), file("${sample}.g.vcf.idx") into gvcfChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type HaplotypeCaller \
        --input_file ${sortedDuplicateMarkedRecalibratedBam} \
        --out ${sample}.g.vcf \
        --reference_sequence ${referenceGenomeDirectory}/${REFERENCE_GENOME_FASTA} \
        --intervals ${intervalsFilePath} \
        --emitRefConfidence GVCF \
        --dbsnp ${resourceBundleDirectory}/${DBSNP_VCF} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process genotypeGVCF {

    container "broadinstitute/gatk3:3.8-0"

    input:
    file gvcf from gvcfChannel.collect()
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path intervalsFilePath from intervalsFilePathChannel
    path resourceBundleDirectory from resourceBundleDirectoryChannel

    output:
    file "${batchName}.joint.vcf" into vcfChannel

    """
    ls *.g.vcf > gvcf.list
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type GenotypeGVCFs \
        --variant gvcf.list \
        --out "${batchName}.joint.vcf" \
        --reference_sequence ${referenceGenomeDirectory}/${REFERENCE_GENOME_FASTA} \
        --intervals ${intervalsFilePath} \
        --dbsnp ${resourceBundleDirectory}/${DBSNP_VCF} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process filterVariants {

    container "broadinstitute/gatk3:3.8-0"
    publishDir "${pipelineOutputPath}"

    input:
    file vcf from vcfChannel
    path referenceGenomeDirectory from referenceGenomeDirectoryChannel
    path intervalsFilePath from intervalsFilePathChannel

    output:
    file "${batchName}.joint.filtered.vcf" into filteredVcfChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
    --analysis_type VariantFiltration \
    --variant ${vcf} \
    --out "${batchName}.joint.filtered.vcf" \
    --reference_sequence ${referenceGenomeDirectory}/${REFERENCE_GENOME_FASTA} \
    --intervals ${intervalsFilePath} \
    --filterExpression "QUAL <= 50.0" \
    --filterName QUALFilter \
    --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process assessAlignmentCoverage {

    container "broadinstitute/gatk3:3.8-0"
    publishDir "${pipelineOutputPath}"

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
    --reference_sequence ${referenceGenomeDirectory}/${REFERENCE_GENOME_FASTA} \
    --intervals ${intervalsFilePath} \
    --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY \
    --minMappingQuality 1 \
    --omitIntervalStatistics \
    --omitPerSampleStats \
    --omitLocusTable
    """

}


process callCYP2D6Alleles {

    container "gcr.io/pharmacogenetics/stargazer:v1.0.8"
    publishDir "${pipelineOutputPath}"

    input:
    file vcf from filteredVcfChannel
    file gdf from depthOfCoverageResultChannel

    output:
    tuple file("${batchName}.stargazer-genotype.txt"),
        file("${batchName}.stargazer-genotype.log"),
        path("${batchName}.stargazer-genotype.project") into stargazerResultsChannel

    """
    mkdir /${batchName}/
    python /usr/Stargazer_v1.0.8/stargazer.py genotype \
    --target_gene cyp2d6 \
    --control_gene egfr \
    --data_type ts \
    --gdf ${gdf} \
    --vcf ${vcf} \
    --output_dir . \
    --output_prefix ${batchName} \
    """

}
