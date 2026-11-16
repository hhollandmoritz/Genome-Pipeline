/*
 * NEW default fast-mode tool: SSUextract (https://github.com/NeLLi-team/ssuextract)
 *
 * A full Nextflow pipeline in its own right (pixi-managed, same pattern as
 * GVClass): detects 16S/18S rRNA genes via bundled Infernal covariance
 * models (RF00177 -> 16S index, RF01960 -> 18S index), then assigns
 * taxonomy via BLAST against SILVA+PR2 (curated profile) or SILVA+PR2+IMG
 * (img profile), with an optional tree-based classification step.
 * Archaea are covered through SILVA's 16S taxonomy, not a separate model.
 */

process SSU_EXTRACT {
    tag "$root"
    label 'ssuextract'
    conda "${projectDir}/envs/pixi_tools.yml"
    publishDir { "${params.outdir}/${root}_ssuextract" }, mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly)
    val install_dir

    output:
    tuple val(root), path("results"), emit: results_dir
    tuple val(root), path("results/cmsearch_summary.tsv"), emit: summary, optional: true
    path "results/extracted/*.fna", emit: extracted_fastas, optional: true

    script:
    def tree_flag = params.ssuextract_tree_classification ? '--tree_classification' : ''
    """
    cd "${install_dir}"
    pixi run ssuextract -q "\$OLDPWD/${mylo_assembly}" --outdir "\$OLDPWD/results" ${tree_flag} \\
        --database_path "${params.db_dir}/ssuextract_install/resources/database" \\
        --database_profile ${params.ssuextract_db_profile}
    """
}

workflow SSU_EXTRACT_ALL {
    take:
    assembly   // tuple(root, mylo_assembly)

    main:
    ch_install_dir = Channel.value(file("${params.db_dir}/ssuextract_install", checkIfExists: true).toAbsolutePath().toString())
    SSU_EXTRACT(assembly, ch_install_dir)

    emit:
    results_dir = SSU_EXTRACT.out.results_dir
    summary     = SSU_EXTRACT.out.summary
}

// ---- cmsearch + extraction, once per sample ----
process CMSEARCH {
    tag "$root"
    label 'cmsearch'
    conda "${projectDir}/envs/infernal.yml"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly)
    path cm_file

    output:
    tuple val(root), path("${root}_extracted_ssu.fa"),         emit: ssu_fa
    tuple val(root), path("${root}_cmsearch_summary.tsv"),     emit: summary
    path "${root}_ssu_hits.tblout",                            emit: tblout

    script:
    """
    cmsearch --cpu ${params.threads_cmsearch} --tblout ${root}_ssu_hits.tblout "${cm_file}" "${mylo_assembly}"

    # cmsearch reports reverse-strand hits with seq-from > seq-to (column
    # 10 is the strand itself). samtools faidx requires start <= end and
    # silently returns an empty sequence otherwise - swap coordinates for
    # '-' strand hits before fetching, then reverse-complement the result.
    #
    # Also builds the summary.tsv REPORT expects (contig_id/taxonomy/sequence)
    # in this same loop, since the extraction already walks every hit once.
    # NOTE: a single Rfam covariance model only confirms "this is eukaryotic
    # SSU rRNA" - it has no BLAST-against-reference step like SSUextract, so
    # it genuinely cannot assign a species/taxon name. The taxonomy column
    # records that limitation honestly rather than fabricating a classification.
    : > ${root}_extracted_ssu.fa
    echo -e "contig_name\\ttaxonomy\\tquery_sequence" > ${root}_cmsearch_summary.tsv

    grep -v "^#" ${root}_ssu_hits.tblout | awk '{print \$1, \$8, \$9, \$10}' | \\
    while read -r tgt tfrom tto strand; do
        if [ "\$strand" = "-" ]; then
            region="\${tgt}:\${tto}-\${tfrom}"
            seq=\$(samtools faidx "${mylo_assembly}" "\$region" | seqkit seq -r -p -t dna)
        else
            region="\${tgt}:\${tfrom}-\${tto}"
            seq=\$(samtools faidx "${mylo_assembly}" "\$region")
        fi
        echo "\$seq" >> ${root}_extracted_ssu.fa

        seq_only=\$(echo "\$seq" | grep -v "^>" | tr -d '\\n')
        echo -e "\${tgt}\\tEukaryota (SSU rRNA, RF01960 - no further classification)\\t\${seq_only}" >> ${root}_cmsearch_summary.tsv
    done
    """
}

workflow CMSEARCH_EUK {
    take:
    assembly   // tuple(root, mylo_assembly)

    main:
    // Tells user how to install db. If rfam_cm database isn't detected
    def model_path = "${params.db_dir}/rfam_cm/RF01960.cm"

    if (!file(model_path).exists()) {
        error """
        ERROR: Required Rfam covariance model not found:
        ${model_path}

        Please install the reference databases first:

            nextflow run main.nf --install_databases true
        """
    }

    ch_model = Channel.value(
        file(model_path, checkIfExists: true)
    )
    CMSEARCH(assembly, ch_model)

    emit:
    ssu_fa  = CMSEARCH.out.ssu_fa
    summary = CMSEARCH.out.summary
    tblout  = CMSEARCH.out.tblout
}
