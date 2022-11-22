This folder contains code for preprocessing CCEP data in BIDS to detect and visually rate ERs (CCEPs).

**Content**

_README.txt 		This file
config_CCEP.m               Matlab script to select datapath. You can run the first part to select the correct datapaths and then use pipeline_CCEP.m to run the rest of this script to configure the right patient settings.
pipeline_CCEP.m             Matlab script to preprocess CCEP data, detect ERs and visually check those detected ERs.
load_ECoGdata.m             Matlab function to load ECoG data from BIDS structure
preprocess_ECoG_spes.m      Matlab function to epoch ECoG data based on events.tsv, avarage the ECoG data. 
detect_n1peak_ECoG_ccep.m   Matlab function to detect the N1 from each CCEP.
visualRating_ccep.m         Matlab function to check the detected CCEPs visually.


