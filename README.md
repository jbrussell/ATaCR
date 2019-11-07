#### Calculate and remove tilt and compliance from 24 hour long SAC files over the course of a deployment for application to ambient noise.

This branch of ATaCR represents a very early version of the codes, with all of the quality control steps stripped out. Currently this is the only version that has been tested on continuous data for the purposes of removing tilt and compliance before calculating ambient noise empirical Green's functions. In theory, this should also be possible with the current version in `master`, but it has not yet been tested.

---

Main changes from `master`:
- Brute force, no quality checks on spectra
- Expects 24 hour long continuous SAC files
- No frequency limit to tilt and compliance removal