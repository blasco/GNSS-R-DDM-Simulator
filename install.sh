# --- Install pyenv https://github.com/pyenv/pyenv ---
# Installation script: https://github.com/pyenv/pyenv-installer
curl -L https://github.com/pyenv/pyenv-installer/raw/master/bin/pyenv-installer | bash

# --- Install python 3.6.5 ---
# First install dependencies 
sudo apt-get install libbz2-dev
sudo apt-get install libssl-dev
sudo apt-get install libreadline-dev
sudo apt-get install libsqlite3-dev
sudo apt-get install tk-dev 
sudo apt-get install python3-tk
# Then install with pyenv
pyenv install 3.6.5

# --- Create Python Virtual Environment with venv module --- 

python -m venv py_venv_gnssr

# Add the following lines to py_venv_gnssr/bin/activate
    pyenv shell 3.6.5 # Set python version
    export PROJECT_SRC_ROOT='/home/woowapdabug/projects/thesis/py_venv_gnssr/src'
    export TDS_ROOT='/home/woowapdabug/projects/thesis/py_venv_gnssr/src/gnssr/tds'
    export CYGNSS_ROOT='/home/woowapdabug/projects/thesis/py_venv_gnssr/src/gnssr/cygnss'

# Activate enviornment 
source py_venv_gnssr/bin/activate

# --- Install packages ---
pip install -r src/python_packages.txt

# Additional libraries
sudo apt-get install libgeos-dev



