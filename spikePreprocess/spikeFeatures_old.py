#!/usr/local/bin/python

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
import scipy
import scipy.linalg

import sonIO
import neuralynxIO


###########################################################################

def spikeFeatures(spikeFiles, chan, fname, spkname, parname = None, timname = None, npca=3, nProbes=4, writeEnergy = False, standardize=False, interleave = 120):
	""" \brief generates spike features
	generates feature file for KlustaKwik
 computed features are, in the order
 energy1 pca1_1 pca2_1 ... pcaN_1 energy2 pca2_2 ... pcaN_2 ... ... time
 if writeEnergy is true, otherwise
 pca1_1 pca2_1 pcaN_1 ... pca1_M pca2_M ... pcaN_M
 (first index running across PC number second index across probes
 \param spikeFiles a list of spike files object (SON for CED or TTFiles to concatenate. The time is shifted
 for each file by the end time fo the previous file + interleave, so that a
 monotone time is obtained
 \param chan the channel number to read
 \param fname the filename for the output, in format that can be used
	as a klustakwik/klusters feature file
 \param nProbes the number of probes per electrode (defaults to 4, for
	tetrodes)
 \param writeEnergy if \a True the "energy" (or L2 norm) is included
	as a feature
 \param standardize if \a True the waveform is normalized to unitary
	L2 norm before PCA calculation
 \param interleave time in seconds to be left between the last spike of one
    fiel and the first spiek of the second one
 """
	
	# the spike file handle
	spkh = file(spkname, 'w')	
	info = spikeFiles[0].info()
	fileType = 'TTfile'
	genData = lambda interleave=0, yieldShiftedTime=False: spikeFiles[0].genTTData(spikeFiles, chan, nProbes=nProbes, interleave = interleave, 
		yieldShiftedTime=yieldShiftedTime)
		
	nPoints = info['nPoints']
	nWvPoints = info['nWvPoints']
	sampleInterval  = info['sampleInterval'] * 1e6
	
	for s in spikeFiles[1:]:
		ii = s.info()
		nPoints += ii.nPoints

	# pre-allocate memory for the features 
	
	nProbePoints =  nWvPoints/nProbes
	meanWv = zeros((nProbes, nProbePoints), float64) # mean nwaveform
	mean2Wv = zeros((nProbes, nProbePoints), float64) # mean squared waveform
	meanWvStd = zeros((nProbes, nProbePoints), float64) # mean standardized waveform
	mean2WvStd = zeros((nProbes, nProbePoints), float64) # mean squared
													 # standardized waveform
	energy = zeros((nPoints, nProbes), float64)
	cov = zeros((nProbes, nProbePoints, nProbePoints), float64)

	# first pass generate mean waveforms and writes spike waveforms in
	# a format suitable for klusters
	
	for (t, wv, istart, iend, wvRaw) in genData():
		np = len(t)
		wv2 = wv*wv
		meanWv = meanWv + sum(wv,0)
		mean2Wv = mean2Wv + sum(wv2,0)
		energy[istart:iend, :] = sqrt(sum(wv2,2))

		if spkh:
			spkh.write(wvRaw.tostring())
		
		if standardize:
			nrg = reshape(energy[istart:iend,:],
						  (len(t), nProbes, 1))
			nrg = repeat(nrg, nProbePoints, 2)
			wvstd = wv.astype(float64) / nrg
			meanWvStd += sum(wvstd,0)
			mean2WvStd += sum(wvstd*wvstd,0)

	
	spkh.close()
	meanWv = meanWv / nPoints
	mean2Wv = mean2Wv / nPoints
	mean2Wv = sqrt(mean2Wv - meanWv*meanWv)
	
	meanWvStd = meanWvStd / nPoints
	mean2WvStd = mean2WvStd / nPoints
	mean2Wvstd = sqrt(mean2WvStd - meanWvStd*meanWvStd)
	
	#second pass generates covariance matrix

	
	for (t, wv, istart, iend, wvRaw) in genData():
		if standardize:
			nrg = reshape(energy[istart:iend,:],
						  (len(t), nProbes, 1))
			nrg = repeat(nrg, nProbePoints, 2)
			wv = wv.astype(float64) / nrg
			wv = wv - meanWvStd
		else:
			wv = wv - meanWv
		
		for i in range(nProbes):
			cov[i,:,:] += dot(transpose(wv[:,i,:]), wv[:,i,:])

	
	cov = cov / nPoints
	
	sd = []
	pc = []
	
	# compute principal components
	for i in range(nProbes):
#		 sds= sqrt(diagonal(cov[i,:,:]))
#		 sd.append(reshape(sds, (len(sds),1)))
#		 cov[i,:,:] /= dot(sd[i], transpose(sd[i]))
		(u, ev, ppc) = scipy.linalg.svd(cov[i,:,:])
		pc.append(transpose(ppc)[:,0:npca])




	
	# third pass generates the principal components and times, both in the
	# actual and shifted version
	
	pca = zeros((nPoints, nProbes, npca), float64)
	tt = zeros((nPoints,1))
	tts = zeros((nPoints,1), float64)
	ffid = zeros((nPoints,1))
	
	for (t, wv, istart, iend, ts, fid, wvRaw) in genData(interleave = interleave,
											 yieldShiftedTime = True):
		#collect all times
		tt[istart:iend] = reshape(t, (iend-istart, 1))
		tts[istart:iend] = reshape(ts, (iend-istart, 1))
		ffid[istart:iend] = reshape(fid, (iend-istart, 1))
		
		for i in range(nProbes):
			pca[istart:iend,i,:] = dot(wv[:,i,:], pc[i])



	
	if writeEnergy:
		nFeatures = nProbes*(npca+1)
		fd = concatenate((reshape(energy, (nPoints, nProbes,1 )),
		pca) , 2)
		
		fd = reshape(fd, (nPoints, nFeatures))
	
	else:
		nFeatures = nProbes*npca
		fd = reshape(pca, (nPoints, nFeatures))

	
	ttsSave = tts # save the version in 1/10000s which will
	tts = tts * 100 / sampleInterval # tt arrives in 1/10000 s convert to microsec

	
	fd = concatenate((fd, tts), 1)
	fd = fd.astype(int32)
	nFeatures += 1
	
	fh = file(fname, 'w')
	
	fh.write("%d\n" % (nFeatures))
	
	fmt_string = ""
	for i in range(nFeatures):
		fmt_string += "%d\t"
	
	fmt_string = fmt_string[:-1] + '\n'

	
	for i in range(nPoints):
		fh.write(fmt_string % tuple((fd[i,:].tolist())))
	
	fh.close()
	
	
	if timname:
		
		times= concatenate((tt, ttsSave, ffid), 1)
		times = times.astype(int32)
		th = file(timname, 'w')
		fmt_string = "%d\t%d\t%d\n"
		for i in range(nPoints):
			th.write(fmt_string % tuple((times[i,:].tolist())))
		th.close()
		
	
	# write parameter file, in the XML format, as described by the klusters 1.6.2 docs

	
	if parname:
		if fileType == 'TTfile':
			makeXMLParFile(spikeFiles[0], parname, nFeatures)
		else:
			raise ValueError("XML param file is not implemented yet for SMR files")


def makeXMLParFile(sf, parname, nFeatures):
	parh = file(parname, "w")
	info = sf.info()
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
	t = doc.createTextNode(str(info['SamplingFrequency'])) 
	elem.appendChild(t)

	elem = doc.createElement("voltageRange")
	acqSys.appendChild(elem)
	adb = info['ADBitVolts']
	adb = float(adb[:(adb.find(' '))])
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
	t = doc.createTextNode("0") # TODO
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
	elemGroup.appendChild(elemChannels)
	t = doc.createTextNode("0")
	elemChan.appendChild(t)

	elemChan = doc.createElement("channel")
	elemGroup.appendChild(elemChannels)
	t = doc.createTextNode("1")
	elemChan.appendChild(t)

	elemChan = doc.createElement("channel")
	elemGroup.appendChild(elemChannels)
	t = doc.createTextNode("2")
	elemChan.appendChild(t)

	elemChan = doc.createElement("channel")
	elemGroup.appendChild(elemChannels)
	t = doc.createTextNode("3")
	elemChan.appendChild(t)

	elemNSamples = doc.createElement("nSamples")
	elemGroup.appendChild(elemNSamples)
	t = doc.createTextNode(str(info['WaveformLength'])) # TODO
	elemNSamples.appendChild(t)

	elemPeakSampleIndex = doc.createElement("peakSampleIndex")
	elemGroup.appendChild(elemPeakSampleIndex)
	t = doc.createTextNode(str(info['AlignmentPt'])) # TODO
	elemPeakSampleIndex.appendChild(t)

	elemNFeatures = doc.createElement("nFeatures")
	elemGroup.appendChild(elemNFeatures)
	t = doc.createTextNode(str(nFeatures)) # TODO
	elemNFeatures.appendChild(t)


	doc.writexml(parh, '\n', '    ')
	parh.close()

if __name__ == '__main__':
	import sys
	name = sys.argv[1]

	
	files = sys.argv[2:]
	spikeFiles = []
	if files[-3:] == 'smr':
		fileType = 'son'
	else:
		fileType = 'nlynx'
	
	
	if fileType == 'son':
		for fname in files:
			spikeFiles.append(sonIO.SONfile(fname))
	elif fileType == 'nlynx':
		for fname in files:
			spikeFiles.append(neuralynxIO.TTFile(fname))
	
	
	
	if fileType == 'son':
		adcMarkChans = spikeFiles[0].channelsOfKind('ADCMarkerChannel')
	elif fileType == 'nlynx':
		adcMarkChans = (0,)
	for chan in adcMarkChans:
		info = spikeFiles[0].info(chan)
		nChanStr  = info['title'][2:]
		fname = name + ".fet." + nChanStr #principal components file name needed
		# by klusters
		spkname = name + ".spk." + nChanStr # spike waveform file, needed by
		# klusters
		parname = name + ".par." + nChanStr
		timname = name + ".tim." + nChanStr
		if fileType == 'son':
			files = []
			for sp in spikesFiles:
				files.append(sonIO.TTfile(sp))
		else:
			files = spikeFiles
		spikeFeatures(files, chan, fname, spkname, parname, timname, npca=3)



	
	
