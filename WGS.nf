#!/usr/bin/env nextflow

// Utils modules
include { extractFastqPairFromDir } from './NextflowModules/Utils/fastq.nf'
include { ExportParams as Workflow_ExportParams } from './NextflowModules/Utils/workflow.nf'

// Mapping modules
include { MEM as BWA_MEM } from './NextflowModules/BWAMEM2/2.2.1/MEM.nf' params(genome_fasta: "$params.genome", optional: "-K 100000000 -Y")  //new best practice settings
include { MarkdupMerge as Sambamba_MarkdupMerge } from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

// GATK BaseRecalibrator
include { BaseRecalibratorApplyBQSR as GATK_BaseRecalibratorApplyBQSR } from './NextflowModules/GATK/4.3.0.0/BaseRecalibrator.nf' params(genome_fasta: "$params.genome", optional_br: "$params.gatk_bqsr_options")
include { ViewUnmapped as Sambamba_ViewUnmapped } from './NextflowModules/Sambamba/0.7.0/ViewUnmapped.nf'
include { Merge as Sambamba_Merge } from './NextflowModules/Sambamba/0.7.0/Merge.nf'

// GATK HaplotypeCaller
include { IntervalListTools as PICARD_IntervalListTools } from './NextflowModules/Picard/3.0.0/IntervalListTools.nf' params(scatter_count:'500', optional:'--BREAK_BANDS_AT_MULTIPLES_OF 1000000')
include { HaplotypeCallerGVCF as GATK_HaplotypeCallerGVCF } from './NextflowModules/GATK/4.3.0.0/HaplotypeCaller.nf' params(genome_fasta: "$params.genome", emit_ref_confidence: "GVCF", optional: "$params.gatk_hc_options")
include { CombineGVCFsInterval as GATK_CombineGVCFsInterval } from './NextflowModules/GATK/4.3.0.0/CombineGVCFs.nf' params(genome_fasta: "$params.genome", optional: "")
include { GenotypeGVCFsInterval as GATK_GenotypeGVCFs } from './NextflowModules/GATK/4.3.0.0/GenotypeGVCFs.nf' params(genome_fasta: "$params.genome", optional: "$params.gatk_ggvcf_options")
include { MergeVcfs as PICARD_MergeVcfs } from './NextflowModules/Picard/3.0.0/MergeVcfs.nf'
include { CombineGVCFs as GATK_CombineGVCFs } from './NextflowModules/GATK/4.3.0.0/CombineGVCFs.nf' params(genome_fasta: "$params.genome", optional: "")

// GATK Filter
include { VariantFiltrationSnpIndel as GATK_VariantFiltration } from './NextflowModules/GATK/4.3.0.0/VariantFiltration.nf' params(
    gatk_path: "$params.gatk_path", genome_fasta: "$params.genome", snp_filter: "$params.gatk_snp_filter", snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter"
)

// Fingerprint modules -> no gatk 4 option
include { UnifiedGenotyper as GATK_UnifiedGenotyper_Fingerprint } from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "--intervals $params.dxtracks_path/$params.fingerprint_target --output_mode EMIT_ALL_SITES")

// QC Modules
include { FastQC } from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional: "")
include { CollectMultipleMetrics as PICARD_CollectMultipleMetrics } from './NextflowModules/Picard/3.0.0/CollectMultipleMetrics.nf' params(genome:"$params.genome", optional: "--PROGRAM CollectAlignmentSummaryMetrics --PROGRAM CollectInsertSizeMetrics --PROGRAM CollectGcBiasMetrics --METRIC_ACCUMULATION_LEVEL SAMPLE")
include { Mosdepth } from './NextflowModules/Mosdepth/0.3.3/Mosdepth.nf' params(optional: "-n --fast-mode")
include { MultiQC } from './NextflowModules/MultiQC/1.14/MultiQC.nf' params(optional: "--config $baseDir/assets/multiqc_config.yaml")

// CustomModules
include { VersionLog } from './CustomModules/Utils/VersionLog.nf'

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

// Define chromosomes used to scatter
def chromosomes = Channel.fromPath(params.genome.replace('fasta', 'dict'))
    .splitCsv(sep:'\t', skip:1)
    .map{type, chr, chr_len, md5, file -> [chr.minus('SN:')]}

workflow {
    // Mapping
    BWA_MEM(fastq_files)
    Sambamba_MarkdupMerge(
        BWA_MEM.out.map{
            sample_id, rg_id, bam_file, bai_file -> [sample_id, bam_file]
        }.groupTuple()
    )

    // GATK BaseRecalibrator
    GATK_BaseRecalibratorApplyBQSR(Sambamba_MarkdupMerge.out.bam_file.combine(chromosomes))
    Sambamba_ViewUnmapped(Sambamba_MarkdupMerge.out.bam_file)
    Sambamba_Merge(GATK_BaseRecalibratorApplyBQSR.out.bam_file.mix(Sambamba_ViewUnmapped.out).groupTuple())

    // GATK HaplotypeCaller (GVCF)
    PICARD_IntervalListTools(Channel.fromPath(params.gatk_hc_interval_list))
    GATK_HaplotypeCallerGVCF(Sambamba_Merge.out.combine(PICARD_IntervalListTools.out.flatten()))

    // Create multisample vcf
    GATK_CombineGVCFsInterval(GATK_HaplotypeCallerGVCF.out.map{
        sample_id, gvcf_file, gvcf_idx_file, interval_file ->
        def interval = interval_file.toRealPath().toString().split("/")[-1]
        [sample_id, gvcf_file, gvcf_idx_file, interval_file, interval]
    }.groupTuple(by: 4).map{
        sample_id, gvcf_files, gvcf_idx_files, interval_file, interval -> [analysis_id, gvcf_files, gvcf_idx_files, interval_file[0]]
    })
    GATK_GenotypeGVCFs(GATK_CombineGVCFsInterval.out)
    PICARD_MergeVcfs(GATK_GenotypeGVCFs.out.groupTuple())

    // Create singlessample g.vcf
    GATK_CombineGVCFs(GATK_HaplotypeCallerGVCF.out.map{sample_id, gvcf_file, gvcf_idx_file, interval_file -> [sample_id, gvcf_file, gvcf_idx_file]}.groupTuple())

    // GATK VariantFiltration
    GATK_VariantFiltration(PICARD_MergeVcfs.out)

    // GATK UnifiedGenotyper (fingerprint)
    GATK_UnifiedGenotyper_Fingerprint(Sambamba_Merge.out)

    // ExonCov
    // ExonCov(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]})

    // QC
    FastQC(fastq_files)
    PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
    Mosdepth(Sambamba_Merge.out)

    MultiQC(analysis_id, Channel.empty().mix(
        FastQC.out,
        Sambamba_MarkdupMerge.out.stats_file,
        PICARD_CollectMultipleMetrics.out,
        Mosdepth.out.txt_files.map{sample_id, txt_files -> txt_files}
    ).collect())

    // Create log files: Repository versions and Workflow params
    VersionLog(Channel.of(
        "${workflow.projectDir}/",
        "${params.dxtracks_path}/",
        "${params.exoncov_path}/",
    ).collect())
    Workflow_ExportParams()
}

// Workflow completion notification
workflow.onComplete {
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: analysis_id,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "WGS Workflow Successful: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html, attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html")
    } else {
        def subject = "WGS Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}

// Custom processes
process ExonCov {
    // Custom process to run ExonCov
    tag {"ExonCov ${sample_id}"}
    label 'ExonCov'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(val(analysis_id), sample_id, path(bam_file), path(bai_file))

    script:
        """
        source ${params.exoncov_path}/venv/bin/activate
        python ${params.exoncov_path}/ExonCov.py import_bam --threads ${task.cpus} --overwrite --exon_bed ${params.dxtracks_path}/${params.exoncov_bed} ${analysis_id} WGS ${bam_file}
        """
}
