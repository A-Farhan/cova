#!/bin/bash

# path for installations
my_path="`pwd`/local"
mkdir -p $my_path
echo "External programs will be installed at $my_path"

# directory of softwares packages
cd softwares

# source directory
sdr="source"
# directory of extracted archives
udr="unpacked"
mkdir -p $udr

### MAFFT ###
echo -e "\n\n`date` -- Proceed with MAFFT installation? [y/n]"
read response

if [ $response == 'y' ]; then
    echo "Okay! Installing MAFFT at $my_path"
    
    if [ ! -s $my_path/bin/mafft ]; then
        
	# untar the source package
        if [ ! -d $udr/mafft-7.453-with-extensions ]; then
            	tar -zvxf $sdr/mafft-7.453-with-extensions-src.tgz -C $udr
	fi
        
	# edit makefile
        cd $udr/mafft-7.453-with-extensions/core/
        echo "Path: `pwd`"
        echo "making a backup copy for the Makefile"
        
	if [ ! -s Makefile.copy ]; then
            cp Makefile Makefile.copy
        else
            echo "a copy of Makefile is already present"
        fi
        
	echo "editing the original Makefile"
        
	function edit_makefile () { 
            sed -i "s+PREFIX = /usr/local+PREFIX = $1+" Makefile; 
        }
        
	edit_makefile $my_path
        # install
        make clean && make && make install
        cd - > /dev/null
    else
        echo -e "\tMAFFT is already installed at the given path\n"
    fi
fi

### FastTree ###
echo -e "\n\n`date` -- Proceed with FastTree installation? [y/n]"
read response
if [ $response == 'y' ]; then
    echo "Okay! Installing FastTree at $my_path"
    if [ ! -s $my_path/bin/FastTree ]; then
        gcc -DUSE_DOUBLE -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o $my_path/bin/FastTree $sdr/FastTree.c -lm
    else
        echo -e "\tFastTree is already installed at the given path\n"
    fi
fi

### HYPHY ###
echo -e "\n\n`date` -- Proceed with Hyphy installation? [y/n]"
read response

if [ $response == 'y' ]; then
    echo "Okay! Installing Hyphy at $my_path"
    
    if [ ! -s $my_path/bin/hyphy ]; then
        
	# untar the source package
        if [ ! -d $udr/hyphy-2.5.8 ]; then
            tar -zvxf $sdr/hyphy-2.5.8.tar.gz -C $udr
        fi
        
	# untar cmake
        if [ ! -d $udr/cmake-3.17.1-Linux-x86_64 ]; then
            tar -xf $sdr/cmake-3.17.1-Linux-x86_64.tar.gz -C $udr
        fi
        
	# install hyphy
        cd $udr/hyphy-2.5.8
	echo "In the hyphy directory. Running cmake"
        ../cmake-3.17.1-Linux-x86_64/bin/cmake -DCMAKE_INSTALL_PREFIX="$my_path" .
        make -j MP
        make install
    else
        echo -e "\tHyphy is already installed at the given path\n"
    fi
fi
exit 0
