"""Feature group compunting peaks"""
import re
from numpy import *
import scipy
import scipy.linalg

from feature import FeatureGroup, listAvailableFeatures


class Peak(FeatureGroup):
	"""Feature Group computing peak"""
	_nChans = 4 
	_featureStems = ('Peak: ',)
	_availableFeatures = []
	for f in _featureStems:
		for c in range(1, _nChans+1):
			_availableFeatures.append(f + str(c))
	_availableFeatures = tuple(_availableFeatures)
	


	def __init__(self, info):
		super(Peak, self).__init__(info)
		self._nPoints = info['nPoints']
		self._nWvPoints = info['nWvPoints']
		self._nProbes = info['nProbes']		
		self._nProbePoints =  self._nWvPoints/self._nProbes
		
		self.peak =  zeros((self._nPoints, self._nProbes), float64)
	
	def pass1(self,  t, wv, istart, iend, ts = None, **dx):
		"""first pass computing peak """
		self.peak[istart:iend, :] = wv.max(2)
	
	def endPass1(self):
		"""empty stub"""
		pass
	
	def finalize(self):
		"""packs the results in the features dict"""
		peakFeatures = sorted([x for x in self._availableFeatures if x.startswith('Peak') ])
		for f in peakFeatures:
			q = (re.match(r'.*(\d+).*', f))
			i = int(q.groups()[0])-1
			self._features[f] = reshape(self.peak[:, i], (self._nPoints, 1))
