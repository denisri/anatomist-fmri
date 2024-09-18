# anatomist-fmri
Anatomist application for fMRI visualization

A dedicated application to display functional MRI based on Anatomist.

BrainVISA/Anatomist should be installed or part of the build project (see https://brainvisa.info/web/download.html). However as currently it requires new features which will be new in Anatomist >= 5.2 (unreleased yet), a developer install with master branches is needed: see https://brainvisa.info/web/download.html#developer-environment-installation

When built, just go (`cd`) to the data directory, and run:

```
ananil.py
```

or to run the "casa-distro container":
```
bv ananil.py
```

For developers, if you are using a Casa-Distro developer container, you can add in the container encvironment `config/bv_maker.cfg`:

in `[ source $CASA_SRC ]`:
```
  git git@github.com:denisri/anatomist-fmri.git main anatomist/anatomist-fmri/master
```

in `[ build $CASA_BUILD ]`:
```
  + $CASA_SRC/anatomist/anatomist-fmri/master
```
then the `bv_maker` build tool will include it in its update/build process.
