#!/usr/bin/env nextflow


/* Command line arguments:
    --samples: Path to directory where the target batch fastq files are stored
    --sdtype: Sequencing data type of taget samples ('ts' for targeted sequencing or 'wgs' for whole genome sequencing)
    --reference: Path to GRCh37 reference genome (with index and .dict files) directory
    --resources: Path to Broad's Institute b37 resource bundle
    --list: Path to .list file with target intervals (stargazer-target-genes.sorted.list)
    --results: Target directory to publish pipeline outputs
*/


sampleReadsFilesChannel = Channel.fromFilePairs("${params.samples}/*_R{1,2}*.fastq.gz")
sequencingDataTypeChannel = Channel.of(params.sdtype)
referenceGenomeDirectoryChannel = Channel.value(params.reference)
resourceBundleDirectoryChannel = Channel.value(params.resources)
resultsDirectory = params.results


targetGenesChannel = Channel.fromList([
    "CACNA1S", "CFTR", "CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2A13",
    "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1",
    "CYP2F1", "CYP2J2", "CYP2R1", "CYP2S1", "CYP2W1", "CYP3A4", "CYP3A5", "CYP3A7",
    "CYP3A43", "CYP4B1", "CYP26A1", "CYP4F2", "CYP19A1", "DPYD", "G6PD", "GSTM1",
    "GSTP1", "GSTT1", "IFNL3", "NAT1", "NAT2", "NUDT15", "POR", "RYR1", "SLC15A2",
    "SLC22A2", "SLCO1B1", "SLCO1B3", "SLCO2B1", "SULT1A1", "TBXAS1", "TPMT", "UGT1A1",
    "UGT1A4", "UGT2B7", "UGT2B15", "UGT2B17", "VKORC1"
])


batchName = file(params.samples, type: "dir").getBaseName()
pipelineOutputPath = "${resultsDirectory}/${batchName}"

// These channels are used later to organize and publish pipeline artifacts per sample
sampleReadsFilesChannel.into {
    samplesFastqsChannel;
    samplesVcfChannel;
    samplesHaplotypesChannel;
    samplesCNVReportsChannel;
}


// Reference genome fasta and index files channels
referenceGenomeFastaChannel = Channel.value(file("${params.reference}/${params.REFERENCE_GENOME_FASTA}"))
referenceGenomeDictChannel = Channel.value(file("${params.reference}/${params.REFERENCE_GENOME_DICT}"))
referenceGenomeAmbChannel = Channel.value(file("${params.reference}/${params.REFERENCE_GENOME_AMB}"))
referenceGenomeAnnChannel = Channel.value(file("${params.reference}/${params.REFERENCE_GENOME_ANN}"))
referenceGenomeBwtChannel = Channel.value(file("${params.reference}/${params.REFERENCE_GENOME_BWT}"))
referenceGenomeFaiChannel = Channel.value(file("${params.reference}/${params.REFERENCE_GENOME_FAI}"))
referenceGenomePacChannel = Channel.value(file("${params.reference}/${params.REFERENCE_GENOME_PAC}"))
referenceGenomeSaChannel = Channel.value(file("${params.reference}/${params.REFERENCE_GENOME_SA}"))


// Resource bundle vcf files channels
dbsnpVcfChannel = Channel.value(file("${params.resources}/${params.DBSNP_VCF}"))
goldStandardIndelsChannel = Channel.value(file("${params.resources}/${params.GOLD_STANDARD_INDELS_1000G_VCF}"))
omni25VcfChannel = Channel.value(file("${params.resources}/${params.OMNI25_1000G_VCF}"))
phase1IndelsVcfChannel = Channel.value(file("${params.resources}/${params.PHASE1_INDELS_1000G_VCF}"))


// Intervals files channels
intervalsFileChannel = Channel.value(file(params.intervals))
intervalsJsonFileChannel = Channel.fromPath(params.intervalsJson)


// Output organizer script channel
outputOrganizerScriptChannel = Channel.fromPath(params.script)


process alignReadFiles {

    container "751848375488.dkr.ecr.us-east-1.amazonaws.com/alignment:v0.1.1"

    input:
    tuple val(sample), file(fastqs) from samplesFastqsChannel
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel

    output:
    tuple val(sample), file("${sample}.sorted.bam") into sortedBamChannel

    """
    bwa mem \
        -K 100000000 \
        -t 16 \
        -R "@RG\\tID:${fastqs[0]}\\tPL:ILLUMINA\\tSM:${sample}" \
        ${referenceGenomeFasta} \
        ${fastqs} \
        | samtools sort -@ 16 -o ${sample}.sorted.bam
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
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel
    file dbsnpVcf from dbsnpVcfChannel
    file goldStandardIndels from goldStandardIndelsChannel
    file omni25Vcf from omni25VcfChannel
    file phase1IndelsVcf from phase1IndelsVcfChannel

    output:
    file "${sample}.recal_data.table" into recalibrationDataChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type BaseRecalibrator \
        --input_file ${sortedDuplicateMarkedBam} \
        --out ${sample}.recal_data.table \
        --reference_sequence ${referenceGenomeFasta} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY \
        --useOriginalQualities \
        --knownSites ${dbsnpVcf} \
        --knownSites ${goldStandardIndels} \
        --knownSites ${omni25Vcf} \
        --knownSites ${phase1IndelsVcf}
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
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel

    output:
    tuple val(sample),
        file("${sample}.sorted.duplicate_marked.recalibrated.bam"),
        file("${sample}.sorted.duplicate_marked.recalibrated.bai") into recalibrationChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type PrintReads \
        --input_file ${sortedDuplicateMarkedBam} \
        --out ${sample}.sorted.duplicate_marked.recalibrated.bam \
        --reference_sequence ${referenceGenomeFasta} \
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
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel
    file dbsnpVcf from dbsnpVcfChannel
    file intervalsFile from intervalsFileChannel

    output:
    tuple file("${sample}.g.vcf"),
        file("${sample}.g.vcf.idx") into gvcfChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type HaplotypeCaller \
        --input_file ${sortedDuplicateMarkedRecalibratedBam} \
        --out ${sample}.g.vcf \
        --reference_sequence ${referenceGenomeFasta} \
        --intervals ${intervalsFile} \
        --emitRefConfidence GVCF \
        --dbsnp ${dbsnpVcf} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process combineGVCFs {

    container "broadinstitute/gatk3:3.8-0"

    input:
    file gvcf from gvcfChannel.collect()
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel
    file dbsnpVcf from dbsnpVcfChannel
    file intervalsFile from intervalsFileChannel

    output:
    tuple file("${batchName}.g.vcf"),
        file("${batchName}.g.vcf.idx") into combinedGVCFChannel

    """
    ls *.g.vcf > gvcfs.list
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type CombineGVCFs \
        --variant gvcfs.list \
        --out "${batchName}.g.vcf" \
        --reference_sequence ${referenceGenomeFasta} \
        --intervals ${intervalsFile} \
        --dbsnp ${dbsnpVcf} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process genotypeGVCF {

    container "broadinstitute/gatk3:3.8-0"

    input:
    file gvcf from combinedGVCFChannel
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel
    file dbsnpVcf from dbsnpVcfChannel
    file intervalsFile from intervalsFileChannel

    output:
    file "${batchName}.vcf" into vcfChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
        --analysis_type GenotypeGVCFs \
        --variant ${gvcf[0]} \
        --out "${batchName}.vcf" \
        --reference_sequence ${referenceGenomeFasta} \
        --intervals ${intervalsFile} \
        --dbsnp ${dbsnpVcf} \
        --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process filterVariants {

    container "broadinstitute/gatk3:3.8-0"

    input:
    file vcf from vcfChannel
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel
    file intervalsFile from intervalsFileChannel

    output:
    file "${batchName}.filtered.vcf" into filteredVcfChannel

    """
    java -jar /usr/GenomeAnalysisTK.jar \
    --analysis_type VariantFiltration \
    --variant ${vcf} \
    --out "${batchName}.filtered.vcf" \
    --reference_sequence ${referenceGenomeFasta} \
    --intervals ${intervalsFile} \
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
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel
    file intervalsFile from intervalsFileChannel

    output:
    file "${batchName}.table" into depthOfCoverageResultChannel

    """
    ls *.bam > bam.list
    java -jar /usr/GenomeAnalysisTK.jar \
    --analysis_type DepthOfCoverage \
    --input_file bam.list \
    --out ${batchName}.table \
    --reference_sequence ${referenceGenomeFasta} \
    --intervals ${intervalsFile} \
    --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY \
    --minMappingQuality 1 \
    --omitIntervalStatistics \
    --omitPerSampleStats \
    --omitLocusTable
    """

}


process callHaplotypes {

    container "751848375488.dkr.ecr.us-east-1.amazonaws.com/stargazer:debian-v1.0.8-new-beagle"
    errorStrategy "ignore" // ¯\_(ツ)_/¯

    input:
    file vcf from stargazerVcfChannel
    file gdf from depthOfCoverageResultChannel
    val sequencingDataType from sequencingDataTypeChannel
    each gene from targetGenesChannel

    output:
    tuple file("${batchName}.${gene}.stargazer-genotype.txt"),
        file("${batchName}.${gene}.stargazer-genotype.log"),
        path("${batchName}.${gene}.stargazer-genotype.project") into stargazerResultsChannel

    """
    mkdir ${batchName}
    python /usr/Stargazer_v1.0.8/stargazer.py genotype \
    --target_gene ${gene} \
    --control_gene EGFR \
    --data_type ${sequencingDataType} \
    --gdf ${gdf} \
    --vcf ${vcf} \
    --output_dir . \
    --output_prefix ${batchName}.${gene}
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
    file referenceGenomeFasta from referenceGenomeFastaChannel
    file referenceGenomeDict from referenceGenomeDictChannel
    file referenceGenomeAmb from referenceGenomeAmbChannel
    file referenceGenomeAnn from referenceGenomeAnnChannel
    file referenceGenomeBwt from referenceGenomeBwtChannel
    file referenceGenomeFai from referenceGenomeFaiChannel
    file referenceGenomePac from referenceGenomePacChannel
    file referenceGenomeSa from referenceGenomeSaChannel

    output:
    tuple val(sample), file("${sample}.filtered.vcf"), file("${sample}.filtered.vcf.idx")

    """
    java -jar /usr/GenomeAnalysisTK.jar \
    --analysis_type SelectVariants \
    --variant ${vcf} \
    --out ${sample}.filtered.vcf \
    --reference_sequence ${referenceGenomeFasta} \
    --sample_name ${sample} \
    --unsafe ALLOW_SEQ_DICT_INCOMPATIBILITY
    """

}


process gatherStargazerResultsPerSample {

    container "751848375488.dkr.ecr.us-east-1.amazonaws.com/pandas:1.0.5"
    publishDir "${pipelineOutputPath}/${sample}", mode: "copy"

    input:
    tuple val(sample), file(_) from samplesHaplotypesChannel
    file stargazerResults from stargazerHaplotypesChannel.collect()
    file outputOrganizerScript from outputOrganizerScriptChannel
    file intervalsJsonFile from intervalsJsonFileChannel

    output:
    file "${sample}.haplotypes.tsv"

    """
    python ${outputOrganizerScript} ${sample} ${intervalsJsonFile}
    """

}


process createSampleCNVReport {

    container "751848375488.dkr.ecr.us-east-1.amazonaws.com/poppler:v0.85.0-2"
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
