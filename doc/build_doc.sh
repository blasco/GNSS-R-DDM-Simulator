rm -r ./source/api_doc
sphinx-apidoc -o ./source/api_doc ./../src/gnssr
make html
chromium build/html/index.html
