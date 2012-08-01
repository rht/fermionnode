set -e

# don't use `which`
#http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
command -v foo >/dev/null 2>&1 || { easy_install pip }

pip install virtualenv
pip install virtualenvwrapper
source /usr/local/share/python/virtualenvwrapper.sh

#we are automatically using the new virtualenv's environment
mkvirtualenv test1

pip install numpy

# assuming that you have fortran, if not
# brew install gfortran
pip install scipy
#pip install -e git+https://github.com/scipy/scipy#egg=scipy-dev


pip install matplotlib
# pip install git+https://github.com/matplotlib/matplotlib.git#matplotlib-dev

