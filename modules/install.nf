/*
 * Genome-Pipeline: Installation processes for the databases and containers.
 *
 * Exports: INSTALL_CONTAINERS(), INSTALL_DATABASES(ch_images), INSTALL_ALL()
 *
 * INSTALL_DATABASES expects tuples (image_name, path_to_sif). For a standalone
 * installation, main.nf should supply existing images from container_dir.
 *
 * Requirements on the installation host: Nextflow, Apptainer, internet access,
 * and permission to build images (root, fakeroot, or an appropriate builder).
 * Install processes must run on the HOST (not inside another container).
 *
 * Warning: Database setup for Pixi-based projects requires networking and can
 * install packages in the writable Nextflow work directory. Pin git revisions
 * and the base image for reproducible production installations.
 */
nextflow.enable.dsl = 2

// All software required across the alternative pipeline branches.
// A single mapping/MetaBAT image can be reused by multiple analysis processes.
process INSTALL_IMAGE {
    tag "$name"
    label 'install'
    storeDir "${params.container_dir}"

    input:
    tuple val(name), val(conda_spec)

    output:
    tuple val(name), path("${name}.sif"), emit: image

    script:
    def packages = conda_spec.replace('bioconda::', '').replace('conda-forge::', '')
    def base_image = params.install_base_image ?: 'continuumio/miniconda3:latest'
    def build_opts = params.apptainer_build_options ?: ''
    """
    cat > ${name}.def <<'APPTAINER_DEFINITION'
Bootstrap: docker
From: ${base_image}

%post
    set -eu
    conda install -y -c conda-forge -c bioconda ${packages} wget
    conda clean -a -y

%environment
    export PATH=/opt/conda/bin:\$PATH
APPTAINER_DEFINITION

    apptainer build ${build_opts} ${name}.sif ${name}.def
    """
}

// CheckM v1 is used indirectly by dRep when genome-quality filtering is on.
process INSTALL_CHECKM1_DB {
    label 'install'
    storeDir "${params.db_dir}/checkm1_db"

    input:
    path image

    output:
    path 'checkm_data', emit: db

    script:
    """
    mkdir -p checkm_data
    apptainer exec "${image}" bash -euc '
        wget -O checkm_data.tar.gz "${params.checkm1_data_url}"
        tar -xzf checkm_data.tar.gz -C checkm_data
        rm checkm_data.tar.gz
    '
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
    mkdir -p checkm2_db
    apptainer exec "${image}" checkm2 database --download --path checkm2_db
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
    mkdir -p gtdbtk_data
    apptainer exec "${image}" bash -euc '
        wget -O gtdbtk_data.tar.gz "${params.gtdbtk_data_url}"
        tar -xzf gtdbtk_data.tar.gz -C gtdbtk_data --strip-components=1
        rm gtdbtk_data.tar.gz
    '
    """
}

// BUSCO is installed as a tool in its own image, but its lineage files
// are also data and must be cached separately for offline HPC use.
process INSTALL_BUSCO_DB {
    label 'install'
    storeDir "${params.db_dir}/busco_data"

    input:
    path image

    output:
    path 'busco_datasets', emit: db

    script:
    """
    mkdir -p busco_datasets
    apptainer exec "${image}" busco --download ${params.busco_lineage ?: 'eukaryota_odb10'} --download_path busco_datasets
    """
}

// The project-specific Pixi environment and DB are installed together,
// matching the original GVCLASS_SETUP implementation.
process INSTALL_GVCLASS_DB {
    label 'install'
    storeDir "${params.db_dir}/gvclass_install"

    input:
    path image

    output:
    path 'gvclass', emit: db

    script:
    """
    apptainer exec "${image}" bash -euc '
        git clone "${params.gvclass_repo_url}" gvclass
        cd gvclass
        pixi install --frozen
        pixi run setup-db
    '
    """
}

// SSUextract is itself a nested Pixi/Nextflow pipeline.
process INSTALL_SSUEXTRACT_DB {
    label 'install'
    storeDir "${params.db_dir}"

    input:
    path image

    output:
    path 'ssuextract_install', emit: db

    script:
    """
    apptainer exec "${image}" bash -euc '
        git clone "${params.ssuextract_repo_url}" ssuextract_install
        cd ssuextract_install
        sed -i "/^[[:space:]]*time[[:space:]]*=/d" config/base.config
        pixi install --frozen
        pixi run setup --database_profile "${params.ssuextract_db_profile ?: 'curated'}"
    '
    """
}

process INSTALL_RFAM_DB {
    label 'install'
    storeDir "${params.db_dir}/rfam_cm"

    output:
    path 'RF01960.cm', emit: db

    script:
    """
    wget -O RF01960.cm "${params.rfam_eukarya_cm_url}"
    """
}

workflow INSTALL_CONTAINERS {
    main:
    ch_specs = Channel.of(
        tuple('myloasm',     'bioconda::myloasm'),
        tuple('mapping',     'bioconda::minimap2 bioconda::samtools'),
        tuple('metabat2',    'bioconda::metabat2'),
        tuple('semibin2',    'bioconda::semibin'),
        tuple('remag',       'bioconda::remag'),
        tuple('drep',        'bioconda::drep conda-forge::python=3.10 conda-forge::setuptools=75.8.2'),
        tuple('checkm2',     'bioconda::checkm2'),
        tuple('busco',       'bioconda::busco'),
        tuple('gtdbtk',      'bioconda::gtdbtk'),
        tuple('gvclass',     'conda-forge::pixi conda-forge::git'),
        tuple('genometools', 'bioconda::genometools-genometools'),
        tuple('ssuextract',  'conda-forge::pixi conda-forge::git'),
        tuple('cmsearch',    'bioconda::infernal bioconda::samtools bioconda::seqkit'),
        tuple('report',      'conda-forge::pandas conda-forge::openpyxl')
    )

    INSTALL_IMAGE(ch_specs)

    emit:
    images = INSTALL_IMAGE.out.image
}

workflow INSTALL_DATABASES {
    take:
    ch_images

    main:
    // These inputs ensure the corresponding image exists before DB setup.
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
