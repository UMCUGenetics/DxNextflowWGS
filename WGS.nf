#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

// Mapping modules
include MEM as BWA_MEM from './NextflowModules/BWA/0.7.17/MEM.nf' params(genome:"$params.genome", optional: '-c 100 -M')
include ViewSort as Sambamba_ViewSort from './NextflowModules/Sambamba/0.7.0/ViewSort.nf'
include MarkdupMerge as Sambamba_MarkdupMerge from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

// GATK BaseRecalibrator
include BaseRecalibrator as GATK_BaseRecalibrator from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/BaseRecalibrator.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "$params.gatk_bqsr_options")
include ViewUnmapped as Sambamba_ViewUnmapped from './NextflowModules/Sambamba/0.7.0/ViewUnmapped.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'

// GATK HaplotypeCaller
include IntervalListTools as PICARD_IntervalListTools from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(scatter_count:'500', optional: 'BREAK_BANDS_AT_MULTIPLES_OF=1000000')
include HaplotypeCallerGVCF as GATK_HaplotypeCallerGVCF from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/HaplotypeCaller.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "$params.gatk_hc_options")
include CatVariantsGVCF as GATK_CatVariantsGVCF from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/CatVariants.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "")
include GenotypeGVCFs as GATK_GenotypeGVCFs from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/GenotypeGVCFs.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "$params.gatk_ggvcf_options")
include CombineVariants as GATK_CombineVariants from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/CombineVariants.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "--assumeIdenticalSamples")
include VariantFiltrationSnpIndel as GATK_VariantFiltration from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/VariantFiltration.nf' params(
    gatk_path: "$params.gatk_path", genome:"$params.genome", snp_filter: "$params.gatk_snp_filter", snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter"
)

// Fingerprint modules
include UnifiedGenotyper as GATK_UnifiedGenotyper from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "--intervals $params.dxtracks_path/$params.fingerprint_target --output_mode EMIT_ALL_SITES")

// QC Modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include MultiQC from './NextflowModules/MultiQC/1.8/MultiQC.nf' params(optional:"--config $baseDir/assets/multiqc_config.yaml")

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

// Define chromosomes used to scatter
def chromosomes = Channel.fromPath(params.genome.replace('fasta', 'dict'))
    .splitCsv(sep:'\t', skip:1)
    .map{type, chr, chr_len, md5, file -> [chr.minus('SN:')]}

workflow {
    // Mapping
    BWA_MEM(fastq_files)
    Sambamba_ViewSort(BWA_MEM.out)
    Sambamba_MarkdupMerge(
        Sambamba_ViewSort.out.map{
            sample_id, rg_id, bam_file, bai_file -> [sample_id, bam_file]
        }.groupTuple()
    )

    // GATK BaseRecalibrator
    GATK_BaseRecalibrator(Sambamba_MarkdupMerge.out.combine(chromosomes))
    Sambamba_ViewUnmapped(Sambamba_MarkdupMerge.out)
    Sambamba_Merge(GATK_BaseRecalibrator.out.mix(Sambamba_ViewUnmapped.out).groupTuple())

    // GATK HaplotypeCaller (GVCF)
    PICARD_IntervalListTools(Channel.fromPath(params.gatk_hc_interval_list))
    GATK_HaplotypeCallerGVCF(Sambamba_Merge.out.combine(PICARD_IntervalListTools.out.flatten()))
    // Create multisample vcf
    GATK_GenotypeGVCFs(GATK_HaplotypeCallerGVCF.out.map{
        sample_id, gvcf_file, gvcf_idx_file, interval_file -> 
        def interval = interval_file.toRealPath().toString().split("/")[-1]
        [sample_id, gvcf_file, gvcf_idx_file, interval_file, interval]
    }.groupTuple(by: 3).map{
        sample_id, gvcf_files, gvcf_idx_files, interval_file, interval -> [analysis_id, gvcf_file, gvcf_idx_file, interval_file[0]]
    })
    GATK_CombineVariants(GATK_GenotypeGVCFs.out.groupTuple())
    // Create singlessample g.vcf
    GATK_CatVariantsGVCF(GATK_HaplotypeCallerGVCF.out.map{sample_id, gvcf_file, gvcf_idx_file, interval_file -> [sample_id, gvcf_file, gvcf_idx_file]}.groupTuple())

    // GATK VariantFiltration
    GATK_VariantFiltration(GATK_CombineVariants.out)

    // GATK UnifiedGenotyper (fingerprint)
    GATK_UnifiedGenotyper(Sambamba_Merge.out)

    // ExonCov
    // ExonCov(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]})

    // COPY_NUMBER
        // CNV_QDNASEQ
        // CNV_FREEC

    // BAF

    // QC
    FastQC(fastq_files)
        // Poststats
        // multiqc

    // Repository versions
    VersionLog()
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
        sendMail(to: params.email, subject: subject, body: email_html, attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html")
    } else {
        def subject = "WGS Workflow Failed: ${analysis_id}"
        sendMail(to: params.email, subject: subject, body: email_html)
    }
}

// Custom processes
process ExonCov {
    // Custom process to run ExonCov
    tag {"ExonCov ${sample_id}"}
    label 'ExonCov'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, sample_id, path(bam_file), path(bai_file))

    script:
        """
        source ${params.exoncov_path}/venv/bin/activate
        python ${params.exoncov_path}/ExonCov.py import_bam --threads ${task.cpus} --overwrite --exon_bed ${params.dxtracks_path}/${params.exoncov_bed} ${analysis_id} ${bam_file}
        """
}

process GetStatsFromFlagstat {
    // Custom process to run get_stats_from_flagstat.pl
    tag {"GetStatsFromFlagstat"}
    label 'GetStatsFromFlagstat'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(flagstat_files)

    output:
        path('run_stats.txt')

    script:
        """
        python ${baseDir}/assets/get_stats_from_flagstat.py ${flagstat_files} > run_stats.txt
        """
}

process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']

    output:
        path('repository_version.log')

    script:
        """
        echo 'DxNextflowWGS' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        echo 'Dx_tracks' >> repository_version.log
        git --git-dir=${params.dxtracks_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        echo 'ExonCov' >> repository_version.log
        git --git-dir=${params.exoncov_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}
