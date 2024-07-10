#Required packages
readxl
changepoint
lubridate
colorRamps

#to install the command is install.packes("package name")
######################################################################

To make a new experiment folder run StartNewExperiment_wSig_S1.R
-This script will copy all of the important functions and the scripts and create all of the empty folders you will need to a directory you specify i.e. newdir
######################################################################
Next in your new directory:
-copy all of your topscan tracking files into the folder rawtrack
-copy all of your AnyMaze TTl files into the folder rawtrack_ttl
-copy all of your fiberphotometry files into the folder rawsig
-copy all of your fiberphotometry ttl files into the folder rawsig_ttl

######################################################################
run MakePointer_wSig_S2.R
-this will output a csv file into pointer_raw folder
-this file will serve as the backbone of the group analysis with the remaining scripts
-you can also make sure that you have all of the ttl and tracking files by visual inspection
- you may want to add additional/missing information into this file manually. Once you have everything you need place it into the pointer folder
######################################################################
run TopScanParse_wSig_S3.R
-this script will parse your raw topscan tracking files, interpolate missing points, and combine the tracking with the ttl pulse files and save the reformatted data in the folder rtrack

######################################################################
run SigParse_wSig_S4
-this script flags all of the behavior from the tracking files in rtrack (the ones generated with TopScanParse_S3) and aligns your neural signals
-Ideally this script should be broken up into several smaller scripts
-it will write csv files for the trial by trial data into csvdir ....
-it will write binary files for all sessions into folder rout. these files will be used to extract summary data with GetSummary_S5.R

######################################################################
GetSummary_wSig_ZS_S5.R
-this script will compile all of the event types and peri-event time histograms based on Z-score and output the data to the folder output_dump
######################################################################
GetSummary_wSig_MM_S6.R
-this script will compile all of the event types and peri-event time histograms based on Min-Max normalization output the data to the folder output_dump
