from distutils.core import setup, Extension

import sys
if len(sys.argv) == 1:
	sys.argv.append('install')


setup (name = 'neuralynxIO',
       version = '0.1',
       author = 'Francesco P. Battaglia',
       description = 'Functions to read Neuralyx Cheetah files, for the moment just TTfiles ',
       packages = ['neuralynx_io'],
       package_dir = {'neuralynx_io': '.'}
       )
