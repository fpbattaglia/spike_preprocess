from setuptools import setup

setup(
    name='spike_preprocess',
    version='0.1',
    packages=['sonIO', 'neuralynx_io', 'spikeFeatures', 'spikePreprocess'],
    url='https://github.com/fpbattaglia/spike_preprocess',
    license='GPLv3',
    author='Francesco Battaglia ',
    author_email='fpbattaglia@gmail.com',
    description='"old" spike sorting pipeline, used e.g. for Arbab et al (2018), Cabral et al. (2014)'
)
