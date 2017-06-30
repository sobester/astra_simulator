# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2016-11-17 17:34:50
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-06-29 12:22:12
from setuptools import setup, find_packages
import subprocess

# Versioning note: I didn't use scm_tools here as I wanted the release to be
# exactly equal to the latest git tag (releasing a new tag updates the version,
# which will in turn be uploaded to anaconda)
git_tag_cmd = "git describe --tags --abbrev=0 | tr -d 'v'"
comm = subprocess.Popen(git_tag_cmd, shell=True, stdout=subprocess.PIPE)
version = comm.communicate()[0].strip().decode("utf-8")


setup(
    name='astra',
    version=version,
    description='ASTRA high altitude balloon flight planner',
    packages=find_packages(),
    include_package_data=True,
    scripts=['bin/astra-sim'],
    # package_data={'': ['*.dat']},
    install_requires=["numpy", "scipy", "grequests", "six", "deap"],
    author='Niccolo Zapponi, Paul Chambers',
    author_email='nz1g10@soton.ac.uk, P.R.Chambers@soton.ac.uk',
)
