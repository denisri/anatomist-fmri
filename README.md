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

