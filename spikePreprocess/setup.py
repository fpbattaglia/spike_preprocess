from distutils.core import setup

import sys
if len(sys.argv) == 1:
	sys.argv.append('install')


setup (name = 'spikePreprocess',
       version = '0.1',
       description = """scripts to extract pre-spike sorting information information from spike2 .smr files
		or Neuralynx .ntt files""",
       author = "Francesco P. Battaglia",
       author_email = "F.P.Battaglia@uva.nl",
       packages = ['spikePreprocess'],
       package_dir = {'spikePreprocess': '.'},
       scripts = ['filter_compute_features.py', 'get_av_waveform.py', 'filter_compute_peak_time.py']
       )
