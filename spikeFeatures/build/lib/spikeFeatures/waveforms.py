"""Feature group compunting the energy, and the first three principal components"""
# TODO add output for mean waveform, adn support for variable number of PCs
import re
from numpy import *
import scipy
import scipy.linalg

from feature import FeatureGroup, listAvailableFeatures


class waveforms(FeatureGroup):
	"""Feature group compunting avg and std waveforms"""
	_nChans = 4 
	_featureStems = ('AvWaveform: ', 'StdWaveform: ')
	_availableFeatures = []
	_parameters = {}
	_paramHelp = {}
	for f in _featureStems:
		for c in range(1, _nChans+1):
			_availableFeatures.append(f + str(c))
	_availableFeatures = tuple(_availableFeatures)
	
	def __init__(self, info):
		super(waveforms, self).__init__(info)
		
		self._nPoints = info['nPoints']
		print('nPoints: ' +str(self._nPoints))
		self._nWvPoints = info['nWvPoints']
		self._nProbes = info['nProbes']		
		self._nProbePoints =  self._nWvPoints/self._nProbes
		
		self.meanWv = zeros((self._nProbes, self._nProbePoints), float64) # mean nwaveform
		self.mean2Wv = zeros((self._nProbes, self._nProbePoints), float64) # mean squared waveform
		
	
	def pass1(self, t, wv, istart, iend, ts = None, **dx):
		"""first pass compute energy and mean waveforms"""
		# first pass generate mean waveforms 
		
		wv2 = wv*wv
		self.meanWv = self.meanWv + sum(wv, 0)
		self.mean2Wv = self.mean2Wv + sum(wv2, 0)
	
	def endPass1(self):
		"""at the end of the first pass the mean waveform is normalized"""
		self.meanWv = self.meanWv / self._nPoints
		self.mean2Wv = self.mean2Wv / self._nPoints
		self.mean2Wv = sqrt(self.mean2Wv - self.meanWv*self.meanWv)
				

	def finalize(self):
		"""packs the results in the features dict"""
		avWaveformFeatures = sorted([x for x in self._availableFeatures  if x.startswith('AvWaveform') ])
		for f in avWaveformFeatures:
			q = (re.match(r'.*(\d+).*', f))
			i = int(q.groups()[0])-1
			self._features[f] = reshape(self.meanWv[i, :], (self._nProbePoints,))
			
		stdWaveformFeatures = sorted([ x for x in self._availableFeatures if x.startswith('StdWaveform') ])
		for f in stdWaveformFeatures:
			q = (re.match(r'.*(\d+).*', f))
			i = int(q.groups()[0])-1
			self._features[f] = reshape(self.mean2Wv[i, :], (self._nProbePoints,))
	
	
