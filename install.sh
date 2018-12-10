# --- Install pyenv https://github.com/pyenv/pyenv ---
# Installation script: https://github.com/pyenv/pyenv-installer
curl -L https://github.com/pyenv/pyenv-installer/raw/master/bin/pyenv-installer | bash

# --- Install python 3.6.5 ---
# First install dependencies 
sudo apt-get install libbz2-dev
sudo apt-get install libssl-dev
sudo apt-get install libsqlite3-dev
sudo apt-get install tk-dev 
sudo apt-get install python3-tk
# Then install with pyenv
pyenv install 3.6.5
# Set python version
pyenv global 3.6.5

# --- Create Python Virtual Environment with venv module --- 
python -m venv py_venv_gnssr

# --- Activate environment ---
source py_venv_gnssr/bin/activate

# --- Install packages ---
pip install -r src/python_packages.txt



