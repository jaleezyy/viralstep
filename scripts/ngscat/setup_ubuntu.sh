#!/bin/bash

#################################################################################################################
############# Next-Generation Sequencing Capture Assessment Tool
############# Developed by Genomics and Bioinformatics Platform of Andalusia, Seville, Spain
############# http://www.bioinfomgp.org/ngscat
############# This script downloads and installs ngsCAT and its dependecies for Ubuntu-based Linux distributions. It has been tested for Linux Mint 12, 13 and 15 and Ubuntu 12.04 and 13.10
############# How to use (superuser mode): sudo -E ./setup_ubuntu.sh <targetAbsolutePath> , for example: sudo -E ./setup_ubuntu.sh /home/user/ngsCAT . Close session and log-in again so that ngsCAT is contained in PATH variable
#################################################################################################################


####### Miscellaneous
targetDirectory=$1 # Directory where ngsCAT and pysam will be downloaded and installed
# Create target directory if not exists
if [ ! -d "$targetDirectory" ]; then
	mkdir -p $targetDirectory
fi

# Update repositories
apt-get update


#apt-get install python2.7 # Uncomment this line if python is not installed. ngsCAT has been tested with python 2.6 and 2.7

###################
# Samtools
###################
# Some dependencies
apt-get install zlibc zlib1g zlib1g-dev
# Install samtools
apt-get install samtools


################
# python-numpy and python-scipy
################
# Install numpy 
apt-get install python-numpy
# Install scipy 
apt-get install python-scipy


################
# python-xlwt
################
apt-get install python-xlwt


################
# python-matplotlib
################
# Install package
apt-get install python-matplotlib


################
# python-pysam: latest release (0.7.5)
################
# Install dependencies for pysam
apt-get install gcc
pysamVersion='pysam-0.7.5.tar.gz'
baseNamePysam='pysam-0.7.5'
newOutput=$targetDirectory'/'$pysamVersion
wget -O $newOutput https://pysam.googlecode.com/files/$pysamVersion

# Uncompress downloaded file
tar -xzvf $newOutput -C $targetDirectory'/'
rm $newOutput

# Install package
currentDir=$PWD
cd $targetDirectory'/'$baseNamePysam
python setup.py install
cd $currentDir


#####################
###### Download and install ngsCAT
#####################
currentVersionNGScat='ngscat.v0.1.tar.gz' # Will be changed for new releases
baseNameNGScat='ngscat.v0.1'
newOutput=$targetDirectory'/'$currentVersionNGScat
# Download ngsCAT from GBPA
wget -O $newOutput http://www.bioinfomgp.org/_media/ngscat/download/$currentVersionNGScat

# Uncompress downloaded file
tar -xzvf $newOutput -C $targetDirectory'/'
rm $newOutput
# Change permissions
chmod +x $targetDirectory'/'$baseNameNGScat'/ngscat.py'

################################
# Add ngsCAT to PATH
################################
echo 'PATH=${PATH}:'$targetDirectory'/'$baseNameNGScat'/'  >> ~/.bashrc
echo 'export PATH' >> ~/.bashrc

# ngsCAT can be run. Check some examples at http://www.bioinfomgp.org/ngscat

