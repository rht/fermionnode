#!bash
set -e

# don't use `which`
#http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
command -v foo >/dev/null 2>&1 || { easy_install pip; }

pip install virtualenv
pip install virtualenvwrapper
source /usr/local/share/python/virtualenvwrapper.sh

# we automatically use the newly created virtualenv's environment
mkvirtualenv ferminode

virtualenv ferminode
pip install numpy

# assuming that you have fortran, if not
# brew install gfortran
pip install scipy
# pip install -e git+https://github.com/scipy/scipy#egg=scipy-dev


# pip install git+https://github.com/matplotlib/matplotlib.git#matplotlib-dev
pip install matplotlib
