"""Feature group compunting the energy, and the first three principal components"""
# TODO add output for mean waveform, adn support for variable number of PCs
import re
from numpy import *
import scipy
import scipy.linalg

from feature import FeatureGroup, listAvailableFeatures


class EnergyPCA(FeatureGroup):
	"""Feature group compunting the energy, and the first three principal components"""
	_nChans = 4 
	_featureStems = ('Energy: ', 'PC1: ', 'PC2: ', 'PC3: ')
	_npca = 3
	_availableFeatures = []
	_parameters = {'standardize': False}
	_paramHelp = {'standardize': "if True, waveforms are normalized by energy before PCA is performed"}
	_doPca = True
	for f in _featureStems:
		for c in range(1, _nChans+1):
			_availableFeatures.append(f + str(c))
	_availableFeatures = tuple(_availableFeatures)
	
	def __init__(self, info):
		super(EnergyPCA, self).__init__(info)
		
		self._nPoints = info['nPoints']
		self._nWvPoints = info['nWvPoints']
		self._nProbes = info['nProbes']		
		self._nProbePoints =  self._nWvPoints/self._nProbes
		
		self.meanWv = zeros((self._nProbes, self._nProbePoints), float64) # mean nwaveform
		self.mean2Wv = zeros((self._nProbes, self._nProbePoints), float64) # mean squared waveform
		self.meanWvStd = zeros((self._nProbes, self._nProbePoints), float64) # mean standardized waveform
		self.mean2WvStd = zeros((self._nProbes, self._nProbePoints), float64) # mean squared
		
		self.energy = zeros((self._nPoints, self._nProbes), float64)
		if self._doPca:
			self.cov = zeros((self._nProbes, self._nProbePoints, self._nProbePoints), float64)
			self.pca = zeros((self._nPoints, self._nProbes, self._npca), float64)
	
	def pass1(self, t, wv, istart, iend, ts = None, **dx):
		"""first pass compute energy and mean waveforms"""
		# first pass generate mean waveforms 
		
		wv2 = wv*wv
		self.meanWv = self.meanWv + sum(wv, 0)
		self.mean2Wv = self.mean2Wv + sum(wv2, 0)
		self.energy[istart:iend, :] = sqrt(sum(wv2, 2))
		if self.standardize:
			nrg = reshape(self.energy[istart:iend, :],
						  (len(t), self._nProbes, 1))
			nrg = repeat(nrg, self._nProbePoints, 2)
			wvstd = wv.astype(float64) / nrg
			self.meanWvStd += sum(wvstd, 0)
			self.mean2WvStd += sum(wvstd*wvstd, 0)
	
	def endPass1(self):
		"""at the end of the first pass the mean waveform is normalized"""
		self.meanWv = self.meanWv / self._nPoints
		self.mean2Wv = self.mean2Wv / self._nPoints
		self.mean2Wv = sqrt(self.mean2Wv - self.meanWv*self.meanWv)
		
		self.meanWvStd = self.meanWvStd / self._nPoints
		self.mean2WvStd = self.mean2WvStd / self._nPoints
		self.mean2Wvstd = sqrt(self.mean2WvStd - self.meanWvStd*self.meanWvStd)
		
	def pass2(self, t, wv, istart, iend, ts = None, **dx):
		"""second pass generates covariance matrix"""
		if not self._doPca:
			return
		if self.standardize:
 			nrg = reshape(self.energy[istart:iend, :],
						  (len(t), self._nProbes, 1))
			nrg = repeat(nrg, self._nProbePoints, 2)
			wv = wv.astype(float64) / nrg
			wv = wv - self.meanWvStd
		else:
			wv = wv - self.meanWv
		
		for i in range(self._nProbes):
			self.cov[i, :, :] += dot(transpose(wv[:, i, :]), wv[:, i, :])
		
	def endPass2(self):
		"""at the end of the second pass, the covariance matrix is normalized, and PCs are computed"""
		if not self._doPca:
			return
		self.cov = self.cov / self._nPoints
		self.pc = []
		
		# compute principal components
		for i in range(self._nProbes):
			(u, ev, ppc) = scipy.linalg.svd(self.cov[i,:,:])
			self.pc.append(transpose(ppc)[:, 0:self._npca])
	
	def pass3(self, t, wv, istart, iend, ts = None, **dx):
		"""third pass generates the principal components""" 
		if not self._doPca:
			return
		
		for i in range(self._nProbes):
			self.pca[istart:iend,i,:] = dot(wv[:,i,:], self.pc[i])
			
	def endpass3(self):
		"""does nothing"""
		pass
	
	def finalize(self):
		"""packs the results in the features dict"""
		energyFeatures = sorted([x for x in self._availableFeatures if x.startswith('Energy') ])
		for f in energyFeatures:
			q = (re.match(r'.*(\d+).*', f))
			i = int(q.groups()[0])-1
			self._features[f] = reshape(self.energy[:, i], (self._nPoints, 1))
		
		if self._doPca:
			pcFeatures = sorted([x for x in self._availableFeatures if x.startswith('PC')])
			for f in pcFeatures:
				q = (re.match(r'.*(\d+).*(\d+).*', f))
				npc = int(q.groups()[0])-1
				ch = int(q.groups()[1])-1
				self._features[f] = reshape(self.pca[:,ch,npc], (self._nPoints, 1))
	
	
class Energy(EnergyPCA):
	"""Feature group compunting the energy, and the first three principal components"""
	_featureStems = ('Energy: ')
	_availableFeatures = []
	_parameters = {'standardize': False}
	_paramHelp = {'standardize': "if True, waveforms are normalized by energy before PCA is performed"}
	_doPca = False
	for f in _featureStems:
		for c in range(1, EnergyPCA._nChans+1):
			_availableFeatures.append(f + str(c))
	_availableFeatures = tuple(_availableFeatures)

