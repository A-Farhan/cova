#!/bin/bash

# path for installations
my_path="`pwd`/local"
mkdir -p $my_path
echo "External programs will be installed at $my_path"

cd softwares/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

### MAFFT ###
echo "`date` -- Proceed with MAFFT installation? [y/n]"
read response
if [ $response == 'y' ]; then
    echo "Okay! Installing MAFFT at $my_path"
    if [ ! -s $my_path/bin/mafft ]; then
        # untar the source package
        if [ ! -d mafft-7.453-with-extensions ]; then
            tar xfvz mafft-7.453-with-extensions-src.tgz
        fi
        # edit makefile
        cd mafft-7.453-with-extensions/core/
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
echo "`date` -- Proceed with FastTree installation? [y/n]"
read response
if [ $response == 'y' ]; then
    echo "Okay! Installing FastTree at $my_path"
    if [ ! -s $my_path/bin/FastTree ]; then
        gcc -DUSE_DOUBLE -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o $my_path/bin/FastTree FastTree.c -lm
    else
        echo -e "\tFastTree is already installed at the given path\n"
    fi
fi

### HYPHY ###
echo "`date` -- Proceed with Hyphy installation? [y/n]"
read response
if [ $response == 'y' ]; then
    echo "Okay! Installing Hyphy at $my_path"
    if [ ! -s $my_path/bin/hyphy ]; then
        # untar the source package
        if [ ! -d hyphy-2.5.8 ]; then
            tar xvfz hyphy-2.5.8.tar.gz
        fi
        # untar cmake
        if [ ! -d cmake-3.17.1-Linux-x86_64 ]; then
            tar xf cmake-3.17.1-Linux-x86_64.tar.gz
        fi
        # install hyphy
        cd hyphy-2.5.8
        ../cmake-3.17.1-Linux-x86_64/bin/cmake -DCMAKE_INSTALL_PREFIX="$my_path" .
        make -j MP
        make install
    else
        echo -e "\tHyphy is already installed at the given path\n"
    fi
fi