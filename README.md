# clustered-permutation-test
This is a SAS macro to perform a clustered permutation test.  The manuscript that describes this software can be found here:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2674116/

The Simdata_bin.sas simulates a binary dataset to demonstrate the macro.

The permuation_macro.sas is a SAS macro to perform a clustered permuation test

By clustered data, I mean the randomization unit is the cluster. Example clusters include: patients clustered within physicians, students clustered within classrooms. This macro applies to designs where the physician/ classroom is randomized to the treatment of interest and outcomes are measured at the patient/ student level. Please note that this is a different design than randomization within cluster.
