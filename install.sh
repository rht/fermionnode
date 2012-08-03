#!bash
set -e

echo "This script will install pip, numpy, scipy, matplotlib for you, and the latest version of gfortran if necessary"

# don't use `which`
#http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
hash pip &> /dev/null
if [ $? -eq 1]; then
    easy_install pip
fi

#pip install virtualenv
#pip install virtualenvwrapper
#source /usr/local/share/python/virtualenvwrapper.sh
# we automatically use the newly created virtualenv's environment
#mkvirtualenv ferminode
#virtualenv ferminode

echo "installing numpy..."
pip install --upgrade numpy

# assuming that you have fortran, if not
hash gfortran &> /dev/null
if [ $? -eq 1]; then
    #curl -L -O http://downloads.sourceforge.net/project/hpc/hpc/gcc/gcc-lion.tar.gz
    #brew install gfortran
    echo "You need to install gfortran ... google 'hpc mac gfortran' for precompiled binary or 'homebrew gfortran'"
    exit
fi


echo "installing scipy..."
pip install --upgrade scipy
# pip install -e git+https://github.com/scipy/scipy#egg=scipy-dev


# pip install git+https://github.com/matplotlib/matplotlib.git#matplotlib-dev
echo "installing matplotlib..."
pip install --upgrade matplotlib
