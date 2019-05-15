"""module defining the parent class for all feature groups"""
# TODO add support for help text for features

import sys

class ParamDict(dict):
	"""helper class defining a dictionary with fixed keys, helping to find mispelled keys early, 
	and containing a further dictionary with help information"""
	class __HelpDict(dict):
		"""member class for the the help dictionary, with read-only keys"""
		def __init__(self, defaultDict):
			super(self.__class__, self).__init__(defaultDict)
			self.writable = True
			
		def __setitem__(self, name, value):
			if self.writable:
				if self.has_key(name):
					dict.__setitem__(self, name, value)
				else:
					raise KeyError("unrecognized parameter " + name)
			else:
				raise KeyError("Help dictionary is read only")
		
	def __init__(self, defaultDictionary, helpDict=None, cls=None):
		super(ParamDict, self).__init__(defaultDictionary)
		self.cls = cls
		self.__helpDict = self.__HelpDict(dict.fromkeys(self.keys()))
		if helpDict:
			v = [x in self.keys() for x in helpDict.keys()]
			if not all(v):
				raise KeyError("help dict contains unrecognized parameters")
			self.__helpDict.update(helpDict)
		self.__helpDict.writable = False
	
	def __setitem__(self, name, value):
		if self.has_key(name):
			dict.__setitem__(self, name, value)
		else:
			raise KeyError("unrecognized parameter " + name)
	
	@property
	def help(self):
		"""returns the help dictionary"""
		return self.__helpDict
	


class FeatureGroup(object):
	"""parent class to be subclassed by all the classes implementing featureGroups"""
	_availableFeatures = ( )
	_parameters = {}
	_paramHelp = {}
	
	@classmethod
	def paramDict(cls):
		"""returns the parameter dictionary"""
		return ParamDict(cls._parameters, cls._paramHelp, cls)
		
	@classmethod 
	def featuresAvailable(cls):
		"""returns the available features"""
		return cls._availableFeatures
	
	@classmethod 	
	def nPasses(cls):
		"""returns the required number of passes through the data, will work for subclasses as well, 
		using a little introspection (looks at which functions are overridden in the subclass)"""
		nPasses = 0
		
		if cls.pass1.im_func != FeatureGroup.pass1.im_func:
			nPasses = 1
			if cls.pass2.im_func != FeatureGroup.pass2.im_func:
				nPasses = 2
				if cls.pass3.im_func != FeatureGroup.pass3.im_func:
					nPasses = 3
		return nPasses
	
	def __init__(self, info, params=None):
		self._features = dict.fromkeys(self._availableFeatures)
		self.info = info
		# parse param dict
		
		# use default values if not provided
		if not params:
			params = self.paramDict()
		if not isinstance(params, ParamDict) or params.cls != self.__class__:
			raise ValueError("params must be a ParamDict for the appropriate class")
		
		for (n, v) in params.iteritems():
			self.__setattr__(n, v)
	def pass1(self, t, wv, istart, iend, ts = None, **dx):
		"""default implementation of pass1, raises a NotImplementedError"""
		raise NotImplementedError
	
	def pass2(self, t, wv, istart, iend, ts = None, **dx):
		"""default implementation of pass2, dos nothing"""
		pass
	
	def pass3(self, t, wv, istart, iend, ts = None, **dx):
		"""default implementation of pass2, dos nothing"""
		pass
	
	def endPass1(self):
		"""default implementation of endPass1, dos nothing"""
		pass
	
	def endPass2(self):
		"""default implementation of endPass2, dos nothing"""
		pass
	
	def endPass3(self):
		"""default implementation of endPass3, dos nothing"""
		pass
	
	def finalize(self):
		"""default implementation of finalize, dos nothing"""
		raise NotImplementedError
	
	@property
	def features(self):
		"""returns the feature dictionary"""
		if self._features is None:
			raise ValueError("Computation of features is not finished yet")
		return self._features
	
	@property
	def isMaskable(self):
		return False
	

def listAvailableFeatures(name):
	"""lists the available features, by looking at subclasses of FeatureGroup in the current workspace"""
	features = {}
	g = sys.modules[name].__dict__
	for (n, i) in g.iteritems():
		try:
			for c in i.__bases__:
				if c == FeatureGroup:
					features[n] = i.featuresAvailable()
		except: 
			pass
	return features
	

class MaskedFeatureGroup(object):
	"""wrapper class for a FeatureGroup, to compute features only on a subset of points, 
	indicated by non-zero values in mask"""
	def __init__(self, computeCls, info, params = None, mask = None):
		nPoints = info['nPoints']
		
		if mask and len(mask) != nPoints:
			raise ValueError("mask must have nPoints element")
		newInfo = info.copy()
		newInfo['nPoints'] = sum(mask != 0)
		self.compute = computeCls(info, params)
		self.mask = mask
	
	def doPass(self, nPass, t, wv, istart, iend, ts = None, **dx):
		"""does a 'filtered' pass  """
		if nPass > self.compute.nPasses():
			return 
		msk = self.mask[istart:iend] != 0 
		msk = reshape(msk, (len(msk), ))
		newT = reshape(t[msk], (sum(msk), 1))
		newWv = wv[msk,...]
		if ts:
			newTs = reshape(t[msk], (sum(msk), 1))
		else:	
			newTs = None
					
		if nPass == 1:
			self.compute.pass1(newT, newWv, self.istart[nPass], self.istart[nPass]+sum(msk), newTs, dx)
		elif nPass == 2:
			self.compute.pass3(newT, newWv, self.istart[nPass], self.istart[nPass]+sum(msk), newTs, dx)
		elif nPass == 3:
			self.compute.pass3(newT, newWv, self.istart[nPass], self.istart[nPass]+sum(msk), newTs, dx)
	def pass1(self, t, wv, istart, iend, ts = None, **dx):
		"""docstring for pass1"""
		self.doPass(self, 1, t, wv, istart, iend, ts = None, **dx)

	def pass2(self, t, wv, istart, iend, ts = None, **dx):
		"""docstring for pass1"""
		self.doPass(self, 2, t, wv, istart, iend, ts = None, **dx)
	
	def pass3(self, t, wv, istart, iend, ts = None, **dx):
		"""docstring for pass1"""
		self.doPass(self, 3, t, wv, istart, iend, ts = None, **dx)

	def __getattr__(self, name):
		"""docstring for __getattr__"""
		return 	getattr(self.compute, name)
	

# d = ParamDict({'a': 1, 'b': 2}, {'a': 'uno', 'b': 'due', 'c': 'tre'})
# 
# print d['a']
# d['b'] = 4
# print d['b']
# 
# print d.h['a']