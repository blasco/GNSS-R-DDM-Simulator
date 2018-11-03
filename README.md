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

# Installing the project module
Based on [this answer](https://stackoverflow.com/questions/714063/importing-modules-from-parent-folder)

Install the module in editable state:
    pip install -e .

# TDS-1 Data FTP access
[ftp data](ftp://ftp.merrbys.co.uk)

    Username: jblasco
    Password: @jYhr=M4rF

# CYGNSS Data
[ftp data](ftp://podaac-ftp.jpl.nasa.gov/allData/cygnss/L1/v2.0)


# Sentinel 2 Data
[Sentinel 2 Data explorer](https://apps.sentinel-hub.com/eo-browser/?lat=28.18924&lng=-88.49811&zoom=16&time=2018-06-27&preset=1_TRUE_COLOR&datasource=Sentinel-2%20L1C)

