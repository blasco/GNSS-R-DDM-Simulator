# Python version managed with pyenv

Set the global version to 3.6.5
python --version should output 3.6.5

To list the available versions
    pyenv versions
To install a version
    pyenv install 3.6.5
To set the global version
    pyenv global 3.6.5
    
# Automatic version management
I've added at the beginning of the bin/activate script
    pyenv global 3.6.5
and in the deactivate function
    pyenv global system

# Python virtual environment
Environment created with
    python -m venv env_name
To activate the environment:
    source /path/to/venv/bin/activate
To deactivate environment:
    deactivate

Once the right version of python is set and the environment is activated we can install packages with:
    pip install pkg_name

# TDS-1 Data FTP access
ftp://ftp.merrbys.co.uk

Username:     jblasco
Password:      @jYhr=M4rF
