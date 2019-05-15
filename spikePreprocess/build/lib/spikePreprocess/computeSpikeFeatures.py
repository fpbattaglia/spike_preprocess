#!/usr/local/bin/python
"""engine for computation of spike pre-processing features.  features are computed by passing through the 
spike files multiple times so that loading the waveform files into memory all at once is not necessary"""


## \file doSpikePreprocess.pt
## \brief methods to compute spike features from spike file, for use
## with Klustakwik clusters
##
## As a stand-alone program, this code is called as follows
## doSpikePreprocess root File1 [File2] [File3] ...
## root is used as the root for the  names of the output files
## the other input arguments form a list of spike files, assumed to have
## the same channel layout, and that are concatenated together, as a
## workaround for the time limitations in Spike2
## The outputs are
## root.fet.n: a feature file for each (tetrode) channel, in a format
## suitable for spike2
## root.spk.n: a file containing all the spike waveforms, in a format
## suitable for klusters
## root.par.n: a parameter file for klusters, in the old text format
## root.tim.n a fiel containing the times of each spike, to be used,
## for example, for conversion back in SON format



from numpy import *
import sonIO
import neuralynx_io
from spikeFeatures import *
import warnings 


###########################################################################

def pickFeatureGroups(features, modName):
	"""find the FeatureGroups, new algorithm"""
	d = listAvailableFeatures(modName)
	groups = d.keys()
#	print groups
	groups.sort(lambda x, y: eval ('len(' + x + '.featuresAvailable()) < len( ' + y + '.featuresAvailable()) '  ))
	fSet = set(features)
	groupsUsed = {}
	featuresFound = []
	for g in groups:
		fg = d[g]
#		print('fg: ' + str(fg))
		for f in features:
			if f in fg:
				featuresFound.extend(fg)
				groupsUsed[f] = g
		if fSet.issubset(set(featuresFound)): # meaning that we're done, don't look for other feature groups, 
			#probably more expensive
			break
#	print('groups' + str(groupsUsed.values()))
#	print('fSet ' + str(fSet))
#	print('featuresFound '+str(set(featuresFound)))
	if not fSet.issubset(set(featuresFound)):
		raise ValueError('Features '+str(fSet-set(featuresFound))[3:]+' not found')
	groups = list(set(groupsUsed.values()))
	gClass = []
	for g in groups:
		gClass.append(eval (g))
	return (gClass, groupsUsed) 



def writeTimFile(fname, features):
	"""write the 'tim' file that is used to resync the data for use in matlab"""
	tt = features['Time in File (1/10000 s)']
	tts = features['Time (1/10000 s)']
	ffid = features['File ID']
	
	(nPoints, dummy) = shape(tt)
	times = concatenate((tt, tts, ffid), 1)
	times = times.astype(int32)
	th = file(fname, 'w')
	fmtString = "%d\t%d\t%d\n"
	for i in range(nPoints):
		th.write(fmtString % tuple((times[i, :].tolist())))
	th.close()


def writeFetFile(fname, features, featuresFet):
	"""write the feature file as needed by KlustaKwik and Klusters"""
	fd = features[featuresFet[0]]
	for f in featuresFet[1:]:
		fd = concatenate((fd, features[f]), 1)
	
	fd = fd.astype(int64)
	(nPoints, nFeatures) = shape(fd) 
		
	fh = file(fname, 'w')
	fh.write("%d\n" % (nFeatures))
	
	fmtString = ""
	for i in range(nFeatures):
		fmtString += "%d\t"
	fmtString = fmtString[:-1] + '\n'
	
	for i in range(int(nPoints)):
		fh.write(fmtString % tuple((fd[i, :].tolist())))
	fh.close()

		
class FeatureComputer(object):
	"""FeatureComputer maintains a list of computed features
	parameters
	spikeFiles: a list of names of spike files 
	spkname: name of the spk fiel (for Klusters), if null, it won't be  written
	params: a dict of parameters to be passed to the feature group modules
	interleave: a time (in seconds) to be used as a "pause" between separate files when the are merged together
	"""
	def __init__(self, spikeFiles, spkname=None, params=None, interleave=120, shiftTimes = None, 
		startTime = None, stopTime = None):
		self.spikeFiles = spikeFiles
		self.spkname = spkname
		if not params:
			params = {}
		self.params = params
		self.interleave = interleave
		self.mask = None
		self.computedFeatures = {}
		self.clean = {}
		self.shiftTimes = shiftTimes
		self.startTime = startTime
		self.stopTime = stopTime
		
		
	def isComputedFeature(self, feature):
		"""returns True if f is a feature computed an in a clean state, false otherwise"""
		if self.computedFeatures.has_key(feature) and self.clean[feature]:
			return True
		else:
			return False
	
	def setMask(self, mask, sync=False):
		"""set a mask (vector of indices of "good spikes"), and applies it to all the featuers that are already 
		computed, if sync is True, then all features that need it will be recomputed """
		maxMask = max(mask)
		if len(self.computedFeatures) > 0:
			v = self.computedFeatures.values()[0]
			nPoints = len(v)
			if maxMask > nPoints - 1 or min(mask) < 0:
				raise ValueError('mask contains negative or out-of-range values')
				
		self.mask = mask
		(groups, groupsUsed) = pickFeatureGroups(self.computedFeatures.keys(), __name__)
		for f in self.computedFeatures.keys():
			v = self.computedFeatures[f]
			v = v[mask]
			self.computedFeatures[f] = v
			if not eval(groupsUsed[f]).isMaskable:
				self.clean[f] = False
		if sync:
			self.computeFeatures(self.computedFeatures.keys()) 
		
	def setMaskFromFeatures(self, featureMask ):
		"""set a mask for a number of features, 
		featureMask is a dict of items where values are pair of points corresponding to upper and lower limits"""
		if len(featureMask) > 0:
			featuresForMask = featureMask.keys()
			self.computeFeatures(featuresForMask)
			a = ones(shape(self.computedFeatures[featuresForMask[0]]), bool)
			for f in featuresForMask:
				if not featureMask[f][0] is None:
					a = a & (self.computedFeatures[f] > featureMask[f][0])
				if not featureMask[f][1] is None:
					a = a & (self.computedFeatures[f] < featureMask[f][1])
			# now make the mask
			mask = where(a)[0]
			self.setMask(mask, sync=True)
	
	def setMaskFromCluFile(self, cluFilename, clu):
		"""set a mask from a clu file for the spikes with an index clu"""
		clufile = fromfile(cluFilename, sep = '\n')
		clufile = clufile[1:]
		mask = (clufile == int(clu)).nonzero()[0]
		mask.tofile('mask.txt', '\n')
#		print shape(mask)
		self.setMask(mask, sync=True)
		
	def computeFeatures(self, featuresToCompute, realignPeaks = False):
		"""generates all the features in a 
	 """
		
		
		# the spike file handle	
		if self.spkname:
			try:
				spkh = file(self.spkname, 'w')	
			except TypeError:
				warnings.warn('Cannot write spike file ' + self.spkname)
				spkh = None
		else:
			spkh = None
		
		# the metadata over the (first) spike file
		info = self.spikeFiles[0].info()
		
		# set as many arguments as possible in the generator so that we don't have to think about it again
		print ("computeSpikeFeatures: shiftTimes: "+ str(self.shiftTimes))
		genData = lambda  yieldShiftedTime = False: \
		self.spikeFiles[0].genTTData(self.spikeFiles, nProbes = info['nProbes'], interleave = 120, \
			yieldShiftedTime = yieldShiftedTime, shiftTimes = self.shiftTimes)
		
		# augment metadata with total (across all spikes files) number of spikes 	
		if not self.mask is None:
			nPoints = len(self.mask)
		else:
			nPoints = info['nPoints']
			for s in self.spikeFiles[1:]:
				ii = s.info()
				nPoints += ii['nPoints']
		info['nPoints'] = nPoints
		featuresReady = []
		for f in featuresToCompute:
			if self.clean.has_key(f) and self.clean[f]:
				featuresReady.append(f)
		for f in featuresReady:
			featuresToCompute.remove(f)
		if len(featuresToCompute) == 0: 
			return 
#		print('featuresToCompute: '+ str(featuresToCompute))
		# create the feature computation groups
		(fGroupsClasses, dummy) = pickFeatureGroups(featuresToCompute, __name__)
		
		# search the feature groups containing the featuresToCompute
		fGroups = []
		for i in fGroupsClasses:
			try:
				fGroups.append(i(info, self.params[i.__name__]))
			except KeyError:
				fGroups.append(i(info))
		# and how many passes we need 
		nPasses = max([x.nPasses for x in fGroups]) 
		
		
		# now run the different passes
		if not self.mask is None:
			nMaskedSpikes = 0
		
		if nPasses > 0 or spkh:
			for (t, wv, istart, iend, ts, fid, wvRaw) in genData(yieldShiftedTime = True):
				if realignPeaks:
					(wv,wvRaw) = self.__realignPeaks(wv)
				if not self.mask is None:
					mskNow = self.mask[(self.mask >= istart) & (self.mask < iend)] - istart
					t = t[mskNow]
					if not ts is None:
						ts = ts[mskNow]
						fid = fid[mskNow]
					wv = wv[mskNow, :, :]
					istart = nMaskedSpikes
					nMaskedSpikes += len(mskNow)
					iend = nMaskedSpikes
				for fg in fGroups:
					fg.pass1(t, wv, istart, iend, ts=ts, fid=fid)
				if spkh:
					if not self.mask is None:
						wvRaw = wvRaw[mskNow]
					spkh.write(wvRaw.tostring())
			if spkh:
				spkh.close()
			for fg in fGroups:
				fg.endPass1()
		
		if not self.mask is None:
			nMaskedSpikes = 0	
		if nPasses > 1:
			for (t, wv, istart, iend, wvRaw) in genData():
				if realignPeaks:
					(wv,wvRaw) = self.__realignPeaks(wv)
				if not self.mask is None:
					mskNow = self.mask[(self.mask >= istart) & (self.mask < iend)] - istart
					t = t[mskNow]
					wv = wv[mskNow,:,:]
					istart = nMaskedSpikes
					nMaskedSpikes += len(mskNow)
					iend = nMaskedSpikes
				for fg in fGroups:
					fg.pass2(t, wv, istart, iend)
			for fg in fGroups:
				fg.endPass2()
		
		if not self.mask is None:
			nMaskedSpikes = 0	
		if nPasses > 2:
			if realignPeaks:
				(wv,wvRaw) = self.__realignPeaks(wv)
			for (t, wv, istart, iend, wvRaw) in genData():
				if not self.mask is None:
					mskNow = self.mask[(self.mask >= istart) & (self.mask < iend)] - istart
					t = t[mskNow]
					wv = wv[mskNow, :, :]
					istart = nMaskedSpikes
					nMaskedSpikes += len(mskNow)
					iend = nMaskedSpikes				
				for fg in fGroups:
					fg.pass3(t, wv, istart, iend)
			for fg in fGroups:
				fg.endPass3()
		
		for fg in fGroups:
			fg.finalize()
		
		# assemble the features dict
		features = {}
		for f in featuresToCompute:
			for fg in fGroups:
				if f in fg.features.keys():
					features[f] = fg.features[f]
		for f in features.keys():
			self.clean[f] = True
		self.computedFeatures.update(features)
	
	def __realignPeaks(self, wv):
		"""docstring for __realignPeaks"""
		from scipy import interpolate
		pi = argmax(wv, 2)
		pm = wv.max(2)
		pChan = argmax(pm,1)
		(ll,sz) = shape(pi)
		#pi = pi[:,pChan]
		pi = pi.flatten()[4*arange(0,ll)+pChan]
		#print pi
		#sys.exit(2)
		pi = pi-8
		sh = shape(wv)
		for i in range(sh[0]):
			for j in range(sh[1]):
				if pi[i] > 0:
					wv[i,j,0:(sh[2]-pi[i])] = wv[i,j,pi[i]:]
					# fill with splines
					x = arange(0, (sh[2]-pi[i]))
					tck = interpolate.splrep(x, wv[i,j,x], s=0, k=1 )
					x = arange((sh[2]-pi[i]), sh[2])
					wv[i,j,x] = interpolate.splev(x, tck, der=0)
				elif pi[i] < 0:
					wv[i,j,-pi[i]:] = wv[i,j,0:(sh[2]+pi[i])]
					# fill with splines
					x = arange(-pi[i], sh[2])
					tck = interpolate.splrep(x,wv[i,j,x], s=0)
					x = arange(0,-pi[i])
					wv[i,j,x] = interpolate.splev(x, tck, der=0)
		wvRaw = transpose(wv, [0, 2, 1])
		wvRaw = reshape(wvRaw, (sh[0], sh[1]*sh[2])).astype(int16)
		return (wv, wvRaw)
	
# write parameter file, in the XML format, as described by the klusters 1.6.2 docs
	
def makeXMLParFile(spikeFiles, parname, nFeatures, extFeatures = None):
	"""make the xml parameter file needed by KlustaKwik"""
	import re
	
	if extFeatures == None:
		extFeatures = {}
		
	parh = file(parname, "w")
	info = spikeFiles[0].info()
	info.update(extFeatures)
	from xml.dom.minidom import Document
	doc = Document()
	param = doc.createElement("parameters")
	doc.appendChild(param)
	
	acqSys = doc.createElement("acquisitionSystem")
	param.appendChild(acqSys)
	
	elem = doc.createElement("nBits")
	acqSys.appendChild(elem)		
	t = doc.createTextNode("16")
	elem.appendChild(t)
	
	elem = doc.createElement("nChannels")
	acqSys.appendChild(elem)
	t = doc.createTextNode("4")
	elem.appendChild(t)
	
	elem = doc.createElement("samplingRate")
	acqSys.appendChild(elem)
	sf = int(float(info['SamplingFrequency']))
	t = doc.createTextNode(str(sf))
	elem.appendChild(t)
	
	elem = doc.createElement("voltageRange")
	acqSys.appendChild(elem)
	adb = info['ADBitVolts']
	mtch = re.match('^(\S*)', adb)
	adb = float(mtch.group(1))
	vr = adb * (2 ** 15) *1000 
	t = doc.createTextNode(str(vr)) 
	elem.appendChild(t)
	
	elem = doc.createElement("amplification")
	acqSys.appendChild(elem)
	t = doc.createTextNode("1000") # amplification is a token value as the NLynx header 
	# is expressed in voltage at the probe values
	elem.appendChild(t)
	
	elem = doc.createElement("offset")
	acqSys.appendChild(elem)
	t = doc.createTextNode("0") # FIXME offset is taken to be zero at least for the moment 
	elem.appendChild(t)
	
	elemAnat = doc.createElement("anatomicalDescription")
	param.appendChild(elemAnat)
	
	elemChanGroups = doc.createElement("channelGroups")
	elemAnat.appendChild(elemChanGroups)
	
	elemGroup = doc.createElement("group")
	elemChanGroups.appendChild(elemGroup)
	
	elemChan = doc.createElement("channel")
	elemGroup.appendChild(elemChan)
	t = doc.createTextNode("0")
	elemChan.appendChild(t)
	
	elemChan = doc.createElement("channel")
	elemGroup.appendChild(elemChan)
	t = doc.createTextNode("1")
	elemChan.appendChild(t)
	
	elemChan = doc.createElement("channel")
	elemGroup.appendChild(elemChan)
	t = doc.createTextNode("2")
	elemChan.appendChild(t)
	
	elemChan = doc.createElement("channel")
	elemGroup.appendChild(elemChan)
	t = doc.createTextNode("3")
	elemChan.appendChild(t)
	
	elemAnat = doc.createElement("spikeDetection")
	param.appendChild(elemAnat)
	
	elemChanGroups = doc.createElement("channelGroups")
	elemAnat.appendChild(elemChanGroups)
	
	elemGroup = doc.createElement("group")
	elemChanGroups.appendChild(elemGroup)
	
	elemChannels = doc.createElement("channels")
	elemGroup.appendChild(elemChannels)
	
	elemChan = doc.createElement("channel")
	elemChannels.appendChild(elemChan)
	t = doc.createTextNode("0")
	elemChan.appendChild(t)
	
	elemChan = doc.createElement("channel")
	elemChannels.appendChild(elemChan)
	t = doc.createTextNode("1")
	elemChan.appendChild(t)
	
	elemChan = doc.createElement("channel")
	elemChannels.appendChild(elemChan)
	t = doc.createTextNode("2")
	elemChan.appendChild(t)
	
	elemChan = doc.createElement("channel")
	elemChannels.appendChild(elemChan)
	t = doc.createTextNode("3")
	elemChan.appendChild(t)
	
	elemNSamples = doc.createElement("nSamples")
	elemGroup.appendChild(elemNSamples)
	t = doc.createTextNode(str(info['WaveformLength'])) 
	elemNSamples.appendChild(t)
	
	elemPeakSampleIndex = doc.createElement("peakSampleIndex")
	elemGroup.appendChild(elemPeakSampleIndex)
	t = doc.createTextNode(str(info['AlignmentPt'])) 
	elemPeakSampleIndex.appendChild(t)
	
	elemNFeatures = doc.createElement("nFeatures")
	elemGroup.appendChild(elemNFeatures)
	t = doc.createTextNode(str(nFeatures)) 
	elemNFeatures.appendChild(t)
	
	
	doc.writexml(parh, '\n', '    ')
	parh.close()


