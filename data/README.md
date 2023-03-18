Here we are only able to provide some of the input data needed to run the analyses.
All other files (> 100MB) are generated through the functions.

The raw GBIF download is available at https://doi.org/10.15468/dl.544ewq

Initial bash processing of the GBIF download was performed to reduce file size, but an added column with 'row ID' can be used to trace back occurrence records to the original data as downloaded from GBIF.

### bash commands
~spiritusanti$ unzip 0188599-210914110416597.zip
~spiritusanti$ rename unzipped file to Traqueos_NeoTropics_raw.csv  
~spiritusanti$ awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13"\t"$22"\t"$23"\t"$33}' Traqueos_NeoTropics_raw.csv > Traqueos_NeoTropics_COR.csv  
~spiritusanti$ wc -l Traqueos_NeoTropics_raw.csv  
~spiritusanti$ wc -l Traqueos_NeoTropics_COR.csv  
~spiritusanti$ awk -F"\t" '{print $13"\t"$22"\t"$23"\t"$37}' Traqueos_NeoTropics_raw.csv > Traqueos_NeoTropics_CODES.csv  
