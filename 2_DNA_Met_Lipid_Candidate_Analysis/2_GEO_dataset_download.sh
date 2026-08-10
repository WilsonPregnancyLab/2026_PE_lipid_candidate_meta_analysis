# Screened GEO accession for DNA methylation datasets in preeclamptic patients using Illumina 450K Methylation Array: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

# Scrolled down to Supplmentary File and downloaded the GSE#####_RAW.tar files of each dataset
## Went to custom, selected all and downloaded

# Create raw directory to save these files on the server

# Transfer files from desktop to server 
## go into the folder that the file is saved to, right-click, open terminal, type: 
cp './name_of_file' yourlogin@fhssuperdome.csu.mcmaster.ca/path_to_new_location_you_want_it_saved 

# Make 1 base directory for each GSE#####_RAW.tar. 
# Copy the raw files into the new base directory (while leaving a copy of the raw files alone as a backup)

# Untar and Unzip files in each of the base directories (this allows you to access the files in the zipped folder)
find ./ -type f -name "*.tar" -exec tar -xf {} \;
find ./ -type f -name "*.gz" -exec gunzip {} \;

# Ensure each base directory only has the .idat files. Move all other files into another directory.