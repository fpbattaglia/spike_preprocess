from distutils.core import setup, Extension


setup (name = 'sonIO',
       version = '0.1',
       author = 'Francesco P. Battaglia',
       description = 'Functions to read spike2 file ',
       packages = ['sonIO'],
       package_dir = {'sonIO': '.'}
       )
