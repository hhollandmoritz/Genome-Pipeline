
/*
 * Genome-Pipeline installation workflows (Nextflow DSL2)
 *
 * INSTALL_CONTAINERS : pull software images
 * INSTALL_DATABASES  : install reference datasets using existing images
 * INSTALL_ALL        : run both, with image-to-database dependencies
 *
 * No definition files, fakeroot or privileged builds are used.
 *
 * IMPORTANT: Five compound images still need verified published URIs.
 * Set their null entries in the manifest below before attempting a full install.
 * Single-package Quay entries use tag lookup as a TEMPORARY convenience;
 * replace quay-auto with pinned docker://...:tag URIs for reproducibility.
 */
nextflow.enable.dsl = 2

// ============================================================================
// CONTAINER INSTALLATION
// ============================================================================
// This process pulls an image from a list of images in the install_containers workflow. 
// It uses singularity to run and pull it.
process INSTALL_IMAGE {
    tag "$name"
    label 'install'
    storeDir "${params.container_dir}"

    input:
    tuple val(name), val(source)

    output:
    tuple val(name), path("${name}.sif"), emit: image
    tuple val(name), path("${name}.source.txt"), emit: provenance

    script:
    def safeSource = "'" + source.replace("'", "'\\''") + "'"
    """
    set -euo pipefail
    source_uri=${safeSource}

    # A temporary convenience for published BioContainers repositories.
    # This looks up a concrete tag when the task runs. Pin the resolved URI
    # in the manifest for reproducibility after first successful testing.
    if [[ "\$source_uri" == quay-auto:* ]]; then
        repository="\${source_uri#quay-auto:}"
        source_uri=\$(python3 - "\$repository" <<'PY_RESOLVE'
import json, sys, urllib.parse, urllib.request
repo = sys.argv[1]
url = ('https://quay.io/api/v1/repository/biocontainers/'
       + urllib.parse.quote(repo, safe='')
       + '/tag/?limit=100&onlyActiveTags=true')
with urllib.request.urlopen(url, timeout=60) as response:
    obj = json.load(response)
tags = [t for t in obj.get('tags', []) if t.get('name') and t['name'] != 'latest']
if not tags:
    raise SystemExit('No published BioContainers tags found for ' + repo)
tags.sort(key=lambda t: t.get('last_modified') or '', reverse=True)
print('docker://quay.io/biocontainers/' + repo + ':' + tags[0]['name'])
PY_RESOLVE
        )
    fi

    printf '%s\\n' "\$source_uri" > ${name}.source.txt
    singularity pull ${name}.sif "\$source_uri"
    test -s ${name}.sif
    """
}

// ============================================================================
// REFERENCE DATA INSTALLATION
// ============================================================================
// CheckM1 is retained for dRep's optional internal genome-quality filtering.
process INSTALL_CHECKM1_DB {
    label 'install'
    storeDir "${params.db_dir}/checkm1_db"

    input:
    path image

    output:
    path 'checkm_data', emit: db

    script:
    """
    set -euo pipefail
    mkdir -p checkm_data
    wget -O checkm_data.tar.gz '${params.checkm1_data_url}'
    tar -xzf checkm_data.tar.gz -C checkm_data
    rm checkm_data.tar.gz
    """
}

process INSTALL_CHECKM2_DB {
    label 'install'
    storeDir "${params.db_dir}/checkm2_data"

    input:
    path image

    output:
    path 'checkm2_db', emit: db

    script:
    """
    set -euo pipefail
    mkdir -p checkm2_db
    singularity exec --bind "\$(pwd):\$(pwd)" "${image}" \\
        checkm2 database --download --path "\$(pwd)/checkm2_db"
    """
}

process INSTALL_GTDBTK_DB {
    label 'install'
    storeDir "${params.db_dir}/gtdbtk_db"

    input:
    path image

    output:
    path 'gtdbtk_data', emit: db

    script:
    """
    set -euo pipefail
    mkdir -p gtdbtk_data
    wget -O gtdbtk_data.tar.gz '${params.gtdbtk_data_url}'
    tar -xzf gtdbtk_data.tar.gz -C gtdbtk_data --strip-components=1
    rm gtdbtk_data.tar.gz
    """
}

process INSTALL_BUSCO_DB {
    label 'install'
    storeDir "${params.db_dir}/busco_data"

    input:
    path image

    output:
    path 'busco_datasets', emit: db

    script:
    """
    set -euo pipefail
    mkdir -p busco_datasets
    singularity exec --bind "\$(pwd):\$(pwd)" "${image}" \\
        busco --download '${params.busco_lineage ?: 'eukaryota_odb10'}' \\
        --download_path "\$(pwd)/busco_datasets"
    """
}

// NOTE: These images must contain BOTH git and pixi. Their Pixi environments
// are installed to writable work directories, never into a read-only SIF.
process INSTALL_GVCLASS_DB {
    label 'install'
    storeDir "${params.db_dir}/gvclass_install"

    input:
    path image

    output:
    path 'gvclass', emit: db

    script:
    """
    set -euo pipefail
    git clone '${params.gvclass_repo_url}' gvclass
    singularity exec --bind "\$(pwd):\$(pwd)" "${image}" \\
        bash -c 'cd gvclass && pixi install --frozen && pixi run setup-db'
    """
}

process INSTALL_SSUEXTRACT_DB {
    label 'install'
    storeDir "${params.db_dir}"

    input:
    path image

    output:
    path 'ssuextract_install', emit: db

    script:
    """
    set -euo pipefail
    git clone '${params.ssuextract_repo_url}' ssuextract_install
    sed -i '/^[[:space:]]*time[[:space:]]*=/d' ssuextract_install/config/base.config
    singularity exec --bind "\$(pwd):\$(pwd)" "${image}" \\
        bash -c 'cd ssuextract_install && pixi install --frozen && pixi run setup --database_profile ${params.ssuextract_db_profile ?: 'curated'}'
    """
}

process INSTALL_RFAM_DB {
    label 'install'
    storeDir "${params.db_dir}/rfam_cm"

    output:
    path 'RF01960.cm', emit: db

    script:
    """
    set -euo pipefail
    wget -O RF01960.cm '${params.rfam_eukarya_cm_url}'
    test -s RF01960.cm
    """
}

// ============================================================================
// WORKFLOWS
// ============================================================================
workflow INSTALL_CONTAINERS {
    main:
    // Hardcoded image inventory. No params.container_uris map is required.
    // TODO: Pin concrete tags for quay-auto entries after the first test.
    // TODO: Replace null for five compound environments with VERIFIED URIs.
    def sources = [
        myloasm:     'quay-auto:myloasm',
        // contains samtools=1.24 and minimap2=2.31-r1302
        mapping:     'docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:ad4c240457048f1d3cc0bd8c8fa39bf27d222eb7-0', 
        metabat2:    'quay-auto:metabat2',
        semibin2:    'quay-auto:semibin',
        remag:       'quay-auto:remag',
        drep:        'quay-auto:drep',
        checkm2:     'quay-auto:checkm2',
        busco:       'quay-auto:busco',
        gtdbtk:      'quay-auto:gtdbtk',
        gvclass:     null, // NEEDS git + pixi in ONE image
        genometools: 'quay-auto:genometools-genometools',
        ssuextract:  null, // NEEDS git + pixi in ONE image
        cmsearch:    null, // NEEDS infernal + samtools + seqkit in ONE image
        report:      null  // NEEDS pandas + openpyxl in ONE image
    ]

    def missing = sources.findAll { name, uri -> !uri }.keySet()
    if (missing) {
        error("Missing published image URI(s) in modules/install.nf for: ${missing.join(', ')}. " +
              "These multi-tool environments need verified combined images before a full install.")
    }

    ch_specs = Channel.fromList(sources.collect { name, uri -> tuple(name, uri) })
    INSTALL_IMAGE(ch_specs)

    emit:
    images     = INSTALL_IMAGE.out.image
    provenance = INSTALL_IMAGE.out.provenance
}

workflow INSTALL_DATABASES {
    take:
    ch_images

    main:
    ch_drep       = ch_images.filter { name, image -> name == 'drep'       }.map { name, image -> image }.first()
    ch_checkm2    = ch_images.filter { name, image -> name == 'checkm2'    }.map { name, image -> image }.first()
    ch_gtdbtk     = ch_images.filter { name, image -> name == 'gtdbtk'     }.map { name, image -> image }.first()
    ch_busco      = ch_images.filter { name, image -> name == 'busco'      }.map { name, image -> image }.first()
    ch_gvclass    = ch_images.filter { name, image -> name == 'gvclass'    }.map { name, image -> image }.first()
    ch_ssuextract = ch_images.filter { name, image -> name == 'ssuextract' }.map { name, image -> image }.first()

    INSTALL_CHECKM1_DB(ch_drep)
    INSTALL_CHECKM2_DB(ch_checkm2)
    INSTALL_GTDBTK_DB(ch_gtdbtk)
    INSTALL_BUSCO_DB(ch_busco)
    INSTALL_GVCLASS_DB(ch_gvclass)
    INSTALL_SSUEXTRACT_DB(ch_ssuextract)
    INSTALL_RFAM_DB()

    emit:
    checkm1    = INSTALL_CHECKM1_DB.out.db
    checkm2    = INSTALL_CHECKM2_DB.out.db
    gtdbtk     = INSTALL_GTDBTK_DB.out.db
    busco      = INSTALL_BUSCO_DB.out.db
    gvclass    = INSTALL_GVCLASS_DB.out.db
    ssuextract = INSTALL_SSUEXTRACT_DB.out.db
    rfam       = INSTALL_RFAM_DB.out.db
}

workflow INSTALL_ALL {
    main:
    INSTALL_CONTAINERS()
    INSTALL_DATABASES(INSTALL_CONTAINERS.out.images)

    emit:
    images     = INSTALL_CONTAINERS.out.images
    checkm1    = INSTALL_DATABASES.out.checkm1
    checkm2    = INSTALL_DATABASES.out.checkm2
    gtdbtk     = INSTALL_DATABASES.out.gtdbtk
    busco      = INSTALL_DATABASES.out.busco
    gvclass    = INSTALL_DATABASES.out.gvclass
    ssuextract = INSTALL_DATABASES.out.ssuextract
    rfam       = INSTALL_DATABASES.out.rfam
}
