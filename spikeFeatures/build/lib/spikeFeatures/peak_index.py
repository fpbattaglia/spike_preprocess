"""Feature group compunting the energy, and the first three principal components"""
# TODO add output for mean waveform, adn support for variable number of PCs
import re
from numpy import *
import scipy
import scipy.linalg

from feature import FeatureGroup, listAvailableFeatures


class PeakIndex(FeatureGroup):
	"""Feature group compunting the energy, and the first three principal components"""
	_nChans = 4 
	_featureStems = ('Peak Index: ',)
	_availableFeatures = []
	for f in _featureStems:
		for c in range(1, _nChans+1):
			_availableFeatures.append(f + str(c))
	_availableFeatures = tuple(_availableFeatures)
	
	def __init__(self, info):
		super(PeakIndex, self).__init__(info)
		
		self._nPoints = info['nPoints']
		self._nWvPoints = info['nWvPoints']
		self._nProbes = info['nProbes']		
		self._nProbePoints =  self._nWvPoints/self._nProbes
		
		self.peakIndex = zeros((self._nPoints, self._nProbes), float64)
	
	def pass1(self, t, wv, istart, iend, ts = None, **dx):
		"""first pass compute peakIndex, finishing the calculations"""
		# first pass generate mean waveforms 
		self.peakIndex[istart:iend, :] = argmax(wv, 2)
	
	def finalize(self):
		"""packs the results in the features dict"""
		for f in sorted(self._availableFeatures):
			q = (re.match(r'.*(\d+).*', f))
			i = int(q.groups()[0])-1
			self._features[f] = reshape(self.peakIndex[:, i], (self._nPoints, 1))		
	
	
