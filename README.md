# TDTR_FDTR_vH4
Same as TDTR_vH3, now with FDTR functionality and tested on Windows (instead of just Mac OS X)

WARNING: This does *NOT* use the latest, most accurate version of the TDTR thermal model. There are subtle inaccuracies or oversimplifications in certain edge cases, which have been addressed in more recent versions on the Cahill group website: http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html

Quoting from Cahill's latest version of the thermal model: "... J. Kimling, January 9, 2017; corrected 1D approximation in earlier versions of bidirectional models; includes temperature measurement at locations other than the heat source and depth dependence of optical absorption; v2 corrects errors in the v1 code from November 2016."

If you need those refinements, you'll have to integrate them into the thermal model that sits at the base of my scripts. Look in the TDTR_REFL and TDTR_TEMP scripts for key differences.

Troubleshooting comments:

0) Automatic fitting is *firmly* discouraged: the superstructure of code I wrote on top of the thermal model is designed to assist with efficient manual fitting and sensitivity analysis, which develops your intuition and avoids several pitfalls regarding low-sensitivity and interdependent parameters.

1) The TDTR thermal model assumes that your system's pump beam is on the delay stage, not the probe beam. If you put the probe beam on the delay stage, the model will not match your data. Your out-of-phase signal will be parabolic instead of linear, and the ratio won't fit. To properly model probe-on-delay-stage systems, go into the REFL_vH4 script, look for the last few lines of code and associated comments, and remove the exponential term that you see there.
