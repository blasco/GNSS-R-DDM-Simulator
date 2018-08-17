# Python version managed with pyenv
Set the global version to 3.6.5
python --version should output 3.6.5

To list the available versions
    pyenv versions
To install a version
    pyenv install 3.6.5
To set the global version
    pyenv global 3.6.5

# Python virtual environment
Environment created with
    python -m venv env_name
To activate the environment:
    source /path/to/venv/bin/activate
To deactiva the environment:
    deactivate

once the right version of python is set and the environment is activated we can install packages with:
    pip install pkg_name
