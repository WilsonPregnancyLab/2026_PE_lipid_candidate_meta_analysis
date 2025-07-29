#Install STAR

#!/bin/bash

##Go to programs directory
mkdir /workspace/lab/wilsonslab/eyerk/STAR
cd /workspace/lab/wilsonslab/eyerk/STAR

##Install directory is the untarred folder after you download

##Download the STAR
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz

##Extract Program
tar -zxvf 2.7.11b.tar.gz

##Move into extracted folder
cd STAR-2.7.11b

##Configure into installation directory (if binary is not already available)
cd STAR/source
make STAR

##Add the 'binaries (bin)' folder to yout PATH in .bashrc to allow command to run in any directory
echo 'export PATH="/workspace/lab/wilsonslab/eyerk/programs/STAR/STAR-2.7.11b/source:$PATH"' >> ~/.bashrc
##Reload bash
source ~/.bashrc

#Test it's working
which STAR
STAR --version
