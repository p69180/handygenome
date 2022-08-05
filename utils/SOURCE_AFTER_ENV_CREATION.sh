#211019
#220718 modified

# must be SOURCED after activating conda environment

trap 'echo Error occured. Aborting. ; trap ERR ; return' ERR

PACKAGE_ROOT=$(dirname $(dirname $(realpath -s $BASH_SOURCE)))
MODIFIED_SEQUENZA_SOURCE_DIR=${PACKAGE_ROOT}/sequenza_source_modified

__CONDA_R_ETC=${CONDA_PREFIX}/lib/R/etc
__CONDA_BIN=${CONDA_PREFIX}/bin
RSCRIPT=${__CONDA_BIN}/Rscript
R=${__CONDA_BIN}/R
CRAN_UNIST="https://cran.biodisk.org/"


_make_Renviron_site_()
{
    true
    #echo "RETICULATE_PYTHON=/opt/miniconda_210914/miniconda3/envs/r-reticulate/bin/python3.6" > ${__CONDA_R_ETC}/Renviron.site
    #cp /home/apps/opt/miniconda_210914/Renviron.site ${__CONDA_R_ETC}/Renviron.site
}

_edit_Makeconf_()
{
    cp ${__CONDA_R_ETC}/Makeconf ${__CONDA_R_ETC}/Makeconf.bak
    awk \
        -v pat="([[:blank:]])(x86[^[:blank:]]+)" \
        -v bin="${__CONDA_BIN}" \
        '{print gensub(pat, "\\1" bin "/" "\\2", "g", $0)}' ${__CONDA_R_ETC}/Makeconf > ${__CONDA_R_ETC}/tmp
    mv -f ${__CONDA_R_ETC}/tmp ${__CONDA_R_ETC}/Makeconf
}

_install_packages_()
{
    ### PACKAGES NEEDED TO BE INSTALLED WITHOUT CONDA ###
    
    # hash (required for copynumbersypark)
#    $RSCRIPT -e "install.packages(\"hash\", repos = \"${CRAN_UNIST}\")" 
    
    # sequenza julab
    $R CMD INSTALL ${MODIFIED_SEQUENZA_SOURCE_DIR}/copynumbersypark
    $R CMD INSTALL ${MODIFIED_SEQUENZA_SOURCE_DIR}/sequenza.v2.julab
    $R CMD INSTALL ${MODIFIED_SEQUENZA_SOURCE_DIR}/sequenza.v3.julab
    
#    # update packages (required for installing DUBStepR, RCAv2, and SCopeLoomR)
#    $RSCRIPT -e "devtools::update_packages(upgrade = \"always\")"
#    
#    # DUBStepR
#    $RSCRIPT -e "devtools::install_github(\"prabhakarlab/DUBStepR\")" 
#    
#    # RCAv2
#    $RSCRIPT -e "remotes::install_github(\"prabhakarlab/RCAv2\")"
#    
#    # SCopeLoomR (conda version is too old)
#    $RSCRIPT -e "devtools::install_github(\"aertslab/SCopeLoomR\", build_vignettes = TRUE)"
#    # SCENIC
#    $RSCRIPT -e "devtools::install_github(\"aertslab/SCENIC\")"
}


### MAIN ###

echo "CONDA_PREFIX is : $CONDA_PREFIX"
read -p "Enter y if you want to proceed " __input__

if [[ $__input__ = y ]] ; then
    _make_Renviron_site_
    _edit_Makeconf_
    _install_packages_
fi

unset __input__ _make_Renviron_site_ _edit_Makeconf_ _install_packages_
