"""Feature group computing the 50%  peak spike width"""

import re
from numpy import *
import scipy

from feature import FeatureGroup, listAvailableFeatures


class SpikeWidth50(FeatureGroup):
	"""Feature group compunting the energy, and the first three principal components"""
	_nChans = 4 
	_featureStems = ('Spike Width50: ',)
	_availableFeatures = []
	for f in _featureStems:
		for c in range(1, _nChans+1):
			_availableFeatures.append(f + str(c))
	_availableFeatures = tuple(_availableFeatures)
	
	def __init__(self, info):
		super(SpikeWidth50, self).__init__(info)
		
		self._nPoints = info['nPoints']
		self._nWvPoints = info['nWvPoints']
		self._nProbes = info['nProbes']		
		self._nProbePoints =  self._nWvPoints/self._nProbes
		
		self.spikeWidth = zeros((self._nPoints, self._nProbes), float64)
	
	def pass1(self, t, wv, istart, iend, ts = None, **dx):
		"""first pass compute peakIndex, finishing the calculations"""
		# dimensions are (nSpike, nChan, nPoint)
		
		peak = wv.max(2)
		spikeWidth = zeros(peak.shape)
		for i in range(peak.shape[0]):
			for ch in range(peak.shape[1]):
				if peak[i, ch] <= 0:
					spikeWidth[i,ch] = 0
					continue
				pp = nonzero(wv[i,ch,:] > 0.5 * peak[i, ch])
				pp = pp[0]
				spikeWidth[i,ch] = pp[-1] - pp[0]
		self.spikeWidth[istart:iend, :] = spikeWidth

	def finalize(self):
		"""packs the results in the features dict"""
		for f in sorted(self._availableFeatures):
			q = (re.match(r'.*(\d+).*', f))
			i = int(q.groups()[0])-1
			self._features[f] = reshape(self.spikeWidth[:, i], (self._nPoints, 1))		
	
	
