// Parâmetros de Entrada
params.reads = "dados/*_{1,2}.fastq.gz"
params.ref = "dados/genome.fasta"
params.outdir = "resultados_prova"

process INDEX_REF {
    container 'staphb/bwa'
    
    input:
        path fasta
        
    output:
        path "${fasta}*", emit: indices_bwa
        
    script:
        """
        bwa index ${fasta}
        """
}

process TRIM {
    container 'staphb/fastp'
    publishDir "${params.outdir}/1_trimming", pattern: '*.html', mode: 'copy'

    input:
        tuple val(id), path(reads)

    output:
        tuple val(id), path("*_clean*.fastq.gz"), emit: reads_limpas
        path "*.html", emit: relatorio_html

    script:
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
              -o ${id}_clean_1.fastq.gz -O ${id}_clean_2.fastq.gz \
              --html ${id}_relatorio.html --thread 2
        """
}

process ALIGN_BWA {
    container 'staphb/bwa' 
        publishDir "${params.outdir}/2_alinhamento", mode: 'copy'
    
    input:
        tuple val(id), path(reads)
        path fasta
        path indices 

    output:
        tuple val(id), path("${id}.sam"), emit: sam_bruto

    script:
        """
        bwa mem -t 4 ${fasta} ${reads[0]} ${reads[1]} > ${id}.sam
        """
}

process SORT_BAM {
    container 'staphb/samtools'
    publishDir "${params.outdir}/2_alinhamento", mode: 'copy'

    input:
        tuple val(id), path(sam)

    output:
        tuple val(id), path("${id}_sorted.bam"), emit: bam_alinhado

    script:
        """
        samtools sort -o ${id}_sorted.bam ${sam}
        samtools index ${id}_sorted.bam
        """
}

process VARIANT_CALLING {
    container 'staphb/bcftools'
    publishDir "${params.outdir}/3_variantes", mode: 'copy'

    input:
        tuple val(id), path(bam)
        path fasta

    output:
        tuple val(id), path("${id}_variantes.vcf"), emit: arquivo_vcf

    script:
        """
        bcftools mpileup -Ou -f ${fasta} ${bam} | \
        bcftools call -mv -Ob -o ${id}_variantes.bcf
        
        bcftools view ${id}_variantes.bcf > ${id}_variantes.vcf
        """
}

process CONSENSUS {
    container 'staphb/bcftools'
    publishDir "${params.outdir}/4_consenso", mode: 'copy'

    input:
        tuple val(id), path(vcf)
        path fasta

    output:
        tuple val(id), path("${id}_consenso.fasta"), emit: genoma_consenso

    script:
        """
        bgzip -c ${vcf} > ${vcf}.gz
        bcftools index ${vcf}.gz
        
        bcftools consensus -f ${fasta} ${vcf}.gz > ${id}_consenso.fasta
        """
}

workflow {
    ch_referencia = Channel.fromPath(params.ref).first()
    ch_reads_brutas = Channel.fromFilePairs(params.reads)

    INDEX_REF(ch_referencia)
    TRIM(ch_reads_brutas)
    ALIGN_BWA(TRIM.out.reads_limpas, ch_referencia, INDEX_REF.out.indices_bwa)
    SORT_BAM(ALIGN_BWA.out.sam_bruto)
    VARIANT_CALLING(SORT_BAM.out.bam_alinhado, ch_referencia)
    CONSENSUS(VARIANT_CALLING.out.arquivo_vcf, ch_referencia)
}