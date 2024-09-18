# anatomist-fmri
Anatomist application for fMRI visualization

A dedicated application to display functional MRI based on Anatomist.

The tool is meant to be used for a specific data layout, a data directory containing:
- a `fMRIprep` directory containg anatomical MRIs, and preprocessed fMRI in a BIDS organization.
- a `fMRIprep/sourcedata/freesurfer` subdirectory containing BIDS-organized anats and meshes from Freesurfer.
- a `first_level` directory containing processed z-maps and contrasts, in a BIDS organization

## Installation

BrainVISA/Anatomist should be installed or part of the build project (see https://brainvisa.info/web/download.html). However as currently it requires new features which will be new in Anatomist >= 5.2 (unreleased yet), a developer install with master branches is needed: see https://brainvisa.info/web/download.html#developer-environment-installation.

## Using it

When things are built or installed, just go (`cd`) to the data directory, and run:

```
ananil.py
```

or to run the "casa-distro container":
```
bv ananil.py
```

## Developing and contributing

For developers, if you are using a Casa-Distro developer container, you can add in the container environment `config/bv_maker.cfg`:

in `[ source $CASA_SRC ]`:
```
  git https://github.com/denisri/anatomist-fmri.git main anatomist/anatomist-fmri/master
```

in `[ build $CASA_BUILD ]`:
```
  + $CASA_SRC/anatomist/anatomist-fmri/master
```
then the `bv_maker` build tool will include it in its update/build process.
