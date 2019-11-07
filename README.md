This branch of ATaCR is built from a pre-release version of the package that is nearly identical to the current `master` branch. Most of the changes made were just for the sake of convenience and ideally will be merged into master eventually.

---

Main changes:
- Specify a station list to loop over rather than running each station individually.
- Rename the files alphanumerically to clarify their order
- Add script to make the list of start times from the event list `a0_make_starttimes.m`
- Add `a1_sac2mat_data.m` and `a3_sac2mat_event.m` to convert local SAC files to the required *.mat format
- A few adjustments to the plotting functions
- Conversion from *.mat files back to SAC at the end `c1_eventmat2sac.m`.