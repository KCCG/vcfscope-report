#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail


# TODO in the below: Remove unneeded R packages


main() {
  # Locate the assets bundle, the location of which varies, depending on whether
  # this is an app or an applet.
  if [[ "$DX_RESOURCES_ID" != "" ]]; then
    # This is an app; fetch assets from the app's private asset container
    DX_ASSETS_ID="$DX_RESOURCES_ID"
  else
    # This is an applet; fetch assets from the parent project
    DX_ASSETS_ID="$DX_PROJECT_CONTEXT_ID"
  fi

  # Setup R
  cd ~
  # This R is Aaron's R-3.2.0 with pre-installed packages, with the following 
  # additional packages pre-installed:
  # CRAN: inline, RSQLite, png, gsalib, BH, RcppEigen, gridExtra, StanHeaders, rstan, coda
  # BioC: VariantAnnotation, GenomicRanges, BSgenome
  dx cat "${DX_ASSETS_ID}:/assets/R-3.2.0.compiled.packages_v2.tar.gz" | tar -zxf -
  export PATH="$PWD/bin:$PATH"
  export RHOME=${HOME} # This is needed to make RScript work, since it was compiled in a different dir.

  # Install additional packages needed for rstan
  # rstan itself is as a pre-compiled binary package, as it's slow to compile from source.
  RPACKAGES=(BH_1.58.0-1.tar.gz RcppEigen_0.3.2.5.1.tar.gz gridExtra_2.0.0.tar.gz StanHeaders_2.8.0.tar.gz rstan_2.8.1_R_x86_64-unknown-linux-gnu.tar.gz coda_0.18-1.tar.gz)
  for RPACKAGE in ${RPACKAGES[*]}; do
    dx download "${DX_ASSETS_ID}:/assets/${RPACKAGE}"
    R CMD INSTALL ${RPACKAGE}
  done

  # Install the nlopt library (needed by nloptr, which is needed by lme4)
  dx download "${DX_ASSETS_ID}:/assets/nlopt-2.4.2.tar.gz"
  tar -xvzf nlopt-2.4.2.tar.gz
  cd nlopt-2.4.2
  ./configure --with-pic
  make && sudo make install
  cd ~

  # Install additional packages needed for lme4 (add RcppEigen_0.3.2.5.1.tar.gz also if rstan install code above is omitted)
  RPACKAGES=(minqa_1.2.4.tar.gz nloptr_1.0.4.tar.gz lme4_1.1-10.tar.gz)
  for RPACKAGE in ${RPACKAGES[*]}; do
    dx download "${DX_ASSETS_ID}:/assets/${RPACKAGE}"
    R CMD INSTALL ${RPACKAGE}
  done

  # Fetch inputs
  dx-download-all-inputs --parallel

  # Create a list of all input files
  cat /dev/null > input_files.txt
  for path in ~/in/rdslist/*/*; do
    echo $path >> input_files.txt
  done

  # Run report
  mkdir -p ~/out/pdf/ ~/out/rds/
  ./performance_report.sh input_files.txt ~/out/pdf/performance_report.pdf ~/out/rds/performance_report.rds

  # upload results
  dx-upload-all-outputs
}
