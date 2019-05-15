"""classes computing several different time-related features"""

from numpy import *

from feature import FeatureGroup

class TimeInFile(FeatureGroup):
	"""Time of spike (in 1/10000 s units) as labeled in each file from the data acquisition system"""
	_availableFeatures = ('Time in File (1/10000 s)', )
	def __init__(self, info):
		super(TimeInFile, self).__init__(info)
		self.nPoints = info['nPoints']
		self.tt = zeros((self.nPoints, 1), float64)
	
	def pass1(self, t, wv, istart, iend, ts = None, **dx):
		"""reads in the time"""
		self.tt[istart:iend] = reshape(t, (iend-istart, 1))
	
	def finalize(self):
		"""packs the time in the features dict"""
		self._features[self._availableFeatures[0]] = self.tt


class TimeKlusters(FeatureGroup):
	"""Compute the time (made monotonic across files) in units suitable for Klusters"""
	_availableFeatures = ('Time (Klusters Units)', )
	
	def __init__(self, info):
		super(TimeKlusters, self).__init__(info)
		self.nPoints = info['nPoints']
		self.tt = zeros((self.nPoints, 1), float64)
		self.conv = 100 / float(info['sampleInterval'] * 1e6)
	
	def pass1(self, t, wv, istart, iend, ts=None, **dx):
		"""reads in and convert the time"""
		self.tt[istart:iend] = reshape(ts, (iend-istart, 1)) * self.conv

	def finalize(self):
		"""packs the time in the features dict"""
		self._features[self._availableFeatures[0]] = self.tt
		
class TimeShifted(FeatureGroup):
	"""Gives the time in 1/10000 s, monotonic across files, plus the id of the files containing the spike"""
	_availableFeatures = ('Time (1/10000 s)', 'File ID', )
	
	def __init__(self, info):
		super(TimeShifted, self).__init__(info)
		self.nPoints = info['nPoints']
		self.tt = zeros((self.nPoints, 1), float64)
		self.fid = zeros((self.nPoints, 1))
	
	def pass1(self, t, wv, istart, iend, ts=None,  **dx):
		"""reads in the times"""
		fid = dx['fid']
		self.tt[istart:iend] = reshape(ts, (iend-istart, 1))
		self.fid[istart:iend] = reshape(fid, (iend-istart, 1))
		
	def finalize(self):
		"""packs the time in the features dict"""
		self._features[self._availableFeatures[0]] = self.tt
		self._features[self._availableFeatures[1]] = self.fid
		
	
	

