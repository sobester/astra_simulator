Installation
============

Astra is currently closed source, and can be installed from the source code only. The easiest (and generally preferable) way to achieve this is to use `pip`,
	
	cd /path/to/simulator    # the directory containing setup.py

	pip install .

or run python on the setuptools file directly:
	
	python setup.py

Developers should use the editable flag in order to ensure that the version imported by the local python environment is the latest version of the code, i.e.

	cd /path/to/simulator    # the directory containing setup.py

	pip install -e .
