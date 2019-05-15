from distutils.core import setup

import sys
if len(sys.argv) == 1:
	sys.argv.append('install')

setup (name = 'cellStats',
       version = '0.1',
       description = """compute spike features for spike sorting pre-processing""",
       author = "Francesco P. Battaglia",
       author_email = "F.P.Battaglia@uva.nl",
       packages = ['spikeFeatures'],
       package_dir = {'spikeFeatures': '.'}
       )
