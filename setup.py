#!/usr/bin/python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension

setup(
    name = 'QD Interface',
    description = 'A wrapper for the QD library.',
    version = '0.1.1 alpha',
    ext_modules = [ Extension('qd', ['qd.c'], libraries=['qd'] ), ],
    url='http://baruchel.hd.free.fr/',
    author='Thomas Baruchel',
    author_email='baruchel@gmx.com',
    )
