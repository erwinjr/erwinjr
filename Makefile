SHELL := /bin/bash -o pipefail po errexit

.PHONY: clean
clean:
	find . -name \*.bak -delete
	find . -name \*.o -delete
	find . -name \*.exe -delete
	find . -name \*.py[cod] -delete
	find . -name __pycache__ -delete
	rm -rf .cache build
	rm -f .coverage .coverage.* junit.xml tmpfile.rc tempfile.rc coverage.xml
	rm -rf .pytest_cache

.PHONY: download-pyqwt
download-pyqwt:
	mkdir -p _src/pyqwt
	curl -L http://prdownloads.sourceforge.net/pyqwt/PyQwt-5.2.0.tar.gz -o PyQwt-5.2.0.tar.gz
	tar -xvf PyQwt-5.2.0.tar.gz --strip-components 1 -C _src/pyqwt
	rm PyQwt-5.2.0.tar.gz
	cd _src/pyqwt/configure
	python configure.py -Q ../qwt-5.2
	make
	make install
