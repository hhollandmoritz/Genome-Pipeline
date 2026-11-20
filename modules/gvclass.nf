
process GVCLASS_SETUP {
    label 'small'
    conda "conda-forge::pixi conda-forge::git"

    output:
    val(install_dir), emit: install_dir

    script:
    install_dir = file("${params.db_dir}/gvclass_install/gvclass").toAbsolutePath().toString()

    """
    INSTALL_DIR="${install_dir}"

    # .setup_complete is created after pixi installation and database setup completes
    # if it already exists, ssuextract will not be reinstalled. If it doesn't
    # exist, creating the folder and cloning the repository will proceed.
    if [ ! -f "\$INSTALL_DIR/.setup_complete" ]; then
        mkdir -p "\$(dirname "\$INSTALL_DIR")"

        if [ ! -d "\$INSTALL_DIR/.git" ]; then
            git clone ${params.gvclass_repo_url} "\$INSTALL_DIR"
        fi

        cd "\$INSTALL_DIR"

        pixi install --frozen
        pixi run setup-db

        touch .setup_complete
    fi
    """
}


/*
 * FNA_PREP prepares drep's output for GVCLASS. It breaks large bins into
 * smaller sub-samples and splits bins into separate folders so that 
 * gvclass can run them on separate nodes, increasing capacity on high
 * compute systems.
 */

process FNA_PREP {
    tag "$root"
    label 'gvclass'
    conda "${projectDir}/envs/genometools.yml"

    input:
    tuple val(root), path(derep_dir)

    output:
    tuple val(root), path("fna_prepped/*"), emit: fna_dirs
    path "subsampled_bins.tsv", emit: subsampled

    script:
    """
    set -euo pipefail
    shopt -s nullglob

    mkdir -p fna_prepped split_work
    : > subsampled_bins.tsv

    max_bins=${params.gvclass_bins_per_job}
    max_size_mb=${params.gvclass_max_bin_size_mb}

    group_number=0
    regular_count=0
    total_bins=0
    split_count=0

    bins=( "${derep_dir}"/*.fa )

    for bin in "\${bins[@]}"; do
        # record information about bins
        bin_name=\$(basename "\$bin")
        bin_stem="\${bin_name%.fa}"
        size_mb=\$(du -m "\$bin" | cut -f1)

        total_bins=\$((total_bins + 1))

        # Split the oversized bins into smaller fragments
        # place in folders called "split_<bin_stem>"
        if (( size_mb > max_size_mb )); then

            # Make directory for split bins
            split_dir="fna_prepped/split_\${bin_stem}"
            mkdir -p "\$split_dir"

            # Make a symlink to bin that we will split (gt will preserve the original bin file)
            local_bin="split_work/$bin_name"
            ln -s "\$(realpath "\$bin")" "\$local_bin"

            gt splitfasta -targetsize "\$max_size_mb" "\$local_bin"

            fragments=( "\$local_bin".* )
            n_parts=\${#fragments[@]}

            # Rename fragments to valid FASTA extensions.
            part=0
            for fragment in "\${fragments[@]}"; do
                part=\$((part + 1))
                printf -v part_name '%s_part_%03d.fa' "\$bin_stem" "\$part"
                mv "\$fragment" "\$split_dir/\$part_name"
            done

            rm -f "\$local_bin"
            split_count=\$((split_count + 1))

            printf '%s\\t%s\\t%s\\t%s\\n' \
                "${root}" "\$bin_name" \
                "\${size_mb}MB" "\${n_parts} parts" \
                >> subsampled_bins.tsv

        # else -> (size_mb <= max_size_mb)
        # Move any smaller bins into groups according to max_bins_per_job parameter
        # subdirectory is called 'group'
        else

            # Regular bins are grouped by maximum bin count.
            if (( regular_count % max_bins == 0 )); then
                group_number=\$((group_number + 1))
                printf -v group_dir \
                    'fna_prepped/group_%03d' "\$group_number"
                mkdir -p "\$group_dir"
            fi

            # we use symlinks rather than copying to save interrum space
            ln -s "\$(realpath "\$bin")" "\$group_dir/\$bin_name"
            regular_count=\$((regular_count + 1))

        fi
    done

    rm -rf split_work

    echo "FNA preparation complete for ${root}"
    echo "Original bins: \$total_bins"
    echo "Oversized bins: \$split_count"
    echo "Regular groups: \$group_number"
    echo "Output directories:"
    find fna_prepped -mindepth 1 -maxdepth 1 -type d | sort
    """
}


process GVCLASS {
    tag "$root"
    label 'gvclass'
    conda "${projectDir}/envs/pixi_tools.yml"
    publishDir { "${params.outdir}/${root}_gvclass" }, mode: params.publish_mode

    input:
    tuple val(root), path(fna_dir)
    val install_dir

    output:
    tuple val(root), path("gvclass_summary.tsv"), emit: summary
    path "gvclass_out", emit: gvclass_dir

     script:
    """
    export GVCLASS_DB="${params.db_dir}/gvclass"
    task_dir="\$(pwd)"
    cd "${install_dir}"
    pixi run ./gvclass -o "\${task_dir}/gvclass_out" -t ${task.cpus} "\${task_dir}/${fna_dir}"
    cd "\${task_dir}"
    cp gvclass_out/gvclass_summary.tsv .
    """
}

process COLLECT_SUBSAMPLE_REPORT {
    label 'small'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path reports, stageAs: 'subsample_report??.tsv'

    output:
    path "gvclass_subsampled_bins_summary.tsv"

    script:
    """
    echo -e "sample\\tbin\\toriginal_size\\tparts" > gvclass_subsampled_bins_summary.tsv
    cat ${reports.join(' ')} >> gvclass_subsampled_bins_summary.tsv 2>/dev/null || true
    """
}

workflow GVCLASS_TAXONOMY {
    take:
    derep_dir   // tuple(root, dir)

    main:
    ch_install_dir = Channel.value(file("${params.software_dir}/gvclass", checkIfExists: true).toAbsolutePath().toString())
    FNA_PREP(derep_dir)
    GVCLASS(FNA_PREP.out.fna_dir, ch_install_dir)
    COLLECT_SUBSAMPLE_REPORT(FNA_PREP.out.subsampled.collect())

    emit:
    gvclass_dir = GVCLASS.out.gvclass_dir
    summary     = GVCLASS.out.summary   // feeds REPORT directly
}
