#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

// Mapping modules
include BWAMapping from './NextflowModules/BWA-Mapping/bwa-0.7.17_samtools-1.9/Mapping.nf' params(genome_fasta: "$params.genome", optional: '-c 100 -M')
include MarkdupMerge as Sambamba_MarkdupMerge from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

// GATK BaseRecalibrator
include BaseRecalibrator as GATK_BaseRecalibrator from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/BaseRecalibrator.nf' params(gatk_path: "$params.gatk_path", genome: "$params.genome", optional_bqsr: "$params.gatk_bqsr_options", optional_pr: "$params.gatk_bqsr_pr_options")
include ViewUnmapped as Sambamba_ViewUnmapped from './NextflowModules/Sambamba/0.7.0/ViewUnmapped.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'

// GATK HaplotypeCaller
include IntervalListTools as PICARD_IntervalListTools from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(scatter_count:'500', optional: 'BREAK_BANDS_AT_MULTIPLES_OF=1000000')
include HaplotypeCallerGVCF as GATK_HaplotypeCallerGVCF from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/HaplotypeCaller.nf' params(gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "$params.gatk_hc_options")
include CatVariantsGVCF as GATK_CatVariantsGVCF from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/CatVariants.nf' params(gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "")
include GenotypeGVCFs as GATK_GenotypeGVCFs from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/GenotypeGVCFs.nf' params(gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "$params.gatk_ggvcf_options")
include CombineVariants as GATK_CombineVariants from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/CombineVariants.nf' params(gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "--assumeIdenticalSamples")
include VariantFiltrationSnpIndel as GATK_VariantFiltration from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/VariantFiltration.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", snp_filter: "$params.gatk_snp_filter", snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter"
)

// Fingerprint modules
include UnifiedGenotyper as GATK_UnifiedGenotyper_Fingerprint from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "--intervals $params.dxtracks_path/$params.fingerprint_target --output_mode EMIT_ALL_SITES")

// CNV modules
include Freec from './NextflowModules/ControlFREEC/11.5/Freec.nf' params(chr_len_file: "$params.freec_chr_len_file", chr_files: "$params.freec_chr_files", gem_mappability_file: "$params.freec_gem_mappability_file", ploidy: "$params.freec_ploidy", window: "$params.freec_window")
include AssessSignificance as Freec_AssessSignificance from './NextflowModules/ControlFREEC/11.5/AssessSignificance.nf'
include MakeGraph as Freec_MakeGraph from './NextflowModules/ControlFREEC/11.5/MakeGraph.nf' params(ploidy:2)

// BAF modules
include UnifiedGenotyper as GATK_UnifiedGenotyper_BAF from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "--intervals $params.baf_snsp_bed --output_mode EMIT_ALL_SITES")

// QC Modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional: "")
include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(genome:"$params.genome", optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=CollectGcBiasMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE")
include EstimateLibraryComplexity as PICARD_EstimateLibraryComplexity from './NextflowModules/Picard/2.22.0/EstimateLibraryComplexity.nf' params(optional:"OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500")
include CollectWgsMetrics as PICARD_CollectWgsMetrics from './NextflowModules/Picard/2.22.0/CollectWgsMetrics.nf' params(genome:"$params.genome", optional: "MINIMUM_MAPPING_QUALITY=20 MINIMUM_BASE_QUALITY=10 COVERAGE_CAP=250")
include MultiQC from './NextflowModules/MultiQC/1.9/MultiQC.nf' params(optional: "--config $baseDir/assets/multiqc_config.yaml")

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

// Define chromosomes used to scatter
def chromosomes = Channel.fromPath(params.genome.replace('fasta', 'dict'))
    .splitCsv(sep:'\t', skip:1)
    .map{type, chr, chr_len, md5, file -> [chr.minus('SN:')]}

workflow {
    // Mapping
    BWAMapping(fastq_files)
    Sambamba_MarkdupMerge(
        BWAMapping.out.map{
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
    }.groupTuple(by: 4).map{
        sample_id, gvcf_files, gvcf_idx_files, interval_file, interval -> [analysis_id, gvcf_files, gvcf_idx_files, interval_file[0]]
    })
    GATK_CombineVariants(GATK_GenotypeGVCFs.out.groupTuple())
    // Create singlessample g.vcf
    GATK_CatVariantsGVCF(GATK_HaplotypeCallerGVCF.out.map{sample_id, gvcf_file, gvcf_idx_file, interval_file -> [sample_id, gvcf_file, gvcf_idx_file]}.groupTuple())

    // GATK VariantFiltration
    GATK_VariantFiltration(GATK_CombineVariants.out)

    // GATK UnifiedGenotyper (fingerprint)
    GATK_UnifiedGenotyper_Fingerprint(Sambamba_Merge.out)

    // ExonCov
    // ExonCov(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]})

    // COPY_NUMBER
    Freec(Sambamba_Merge.out)
    Freec_AssessSignificance(Freec.out.cnv)
    Freec_MakeGraph(Freec.out.cnv)
    QDNAseq(Sambamba_Merge.out)

    // BAF
    GATK_UnifiedGenotyper_BAF(Sambamba_Merge.out)
    BAF(GATK_UnifiedGenotyper_BAF.out)

    // QC
    FastQC(fastq_files)
    PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
    PICARD_EstimateLibraryComplexity(Sambamba_Merge.out)
    PICARD_CollectWgsMetrics(Sambamba_Merge.out)
    
    MultiQC(analysis_id, Channel.empty().mix(
        FastQC.out,
        PICARD_CollectMultipleMetrics.out,
        PICARD_EstimateLibraryComplexity.out,
        PICARD_CollectWgsMetrics.out
    ).collect())

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

process QDNAseq {
    // Custom process to run QDNAseq
    tag {"QDNAseq ${sample_id}"}
    label 'QDNAseq'
    container = '/hpc/diaggen/software/singularity_cache/QDNAseq_v1.9.2-HMF.1.sif'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(sample_id, path(bam_file), path(bai_file))
    
    output:
        tuple(sample_id, path("${sample_id}.vcf"), emit: vcf)
        tuple(sample_id, path("${sample_id}.*.igv"), emit: igv)
        tuple(sample_id, path("${sample_id}.readCountsFiltered.rds"), emit: rds)
    
    script:
        """
        Rscript ${baseDir}/assets/run_QDNAseq.R -s ${sample_id} -b ${bam_file}
        mv readCountsFiltered.rds ${sample_id}.readCountsFiltered.rds
        mv copynumber.igv ${sample_id}.copynumber.igv
        mv segments.igv ${sample_id}.segments.igv
        mv calls.igv ${sample_id}.calls.igv
        """
}

process BAF {
    // Custom process to run BAF analysis
    tag {"BAF ${sample_id}"}
    label 'BAF'
    container = '/hpc/diaggen/software/singularity_cache/baf_nextflow.sif'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(sample_id, path(vcf_file))
    
    output:
        tuple(sample_id, path("${sample_id}_BAF.txt"), path("${sample_id}_BAF.pdf"))
    
    script:
        """
        cat ${vcf_file} | bio-vcf --num-threads ${task.cpus} -i \
        --sfilter '!s.empty? and s.dp>=20' \
        --eval '[r.chrom,r.pos,r.ref+">"+r.alt[0]]' \
        --seval 'tot=s.ad.reduce(:+) ; ((tot-s.ad[0].to_f)/tot).round(2)' \
        > ${sample_id}_BAF.txt

        Rscript ${baseDir}/assets/makeBAFplot.R ./ ${sample_id}_BAF.txt

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
