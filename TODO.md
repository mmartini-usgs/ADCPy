Known issues:

Need to Compute depths of cells and add depth as an ordinate variables.
- We will not do this until we convert to EPIC.

It is possible for the 5th beam of a V to have a different number of cells.  This is trapped, but not dealt with, in the code.  In this case 5th beam data will nto be output.

There may be a case where there is only wave data in a Velocity output file.  THe program will probably not handle this well.

Learn how to implment and output log file.

Need process to average very large data sets 6 GB files).  So far working on 300 MB files of some 300k points.