#!/usr/bin/python
# -*- coding: utf-8 -*-

myversion = '0.1.6 alpha'

# Without Numpy
# =============
from distutils.core import setup, Extension

setup(
    name = 'QD Interface',
    description = 'A wrapper for the QD library.',
    version = myversion,
    ext_modules = [ Extension('qd', ['qd.c'], libraries=['qd'] ), ],
    url='http://baruchel.hd.free.fr/',
    author='Thomas Baruchel',
    author_email='baruchel@gmx.com',
    )

# With Numpy
# ==========
# from numpy.distutils.core import setup
# 
# def configuration(parent_package = '', top_path=None):
#     from numpy.distutils.misc_util import Configuration
#     config = Configuration('qd', parent_package, top_path)
#     config.add_extension('qd', sources = ['qd.c'], libraries = ['qd'])
#     return config
# setup(configuration=configuration)
