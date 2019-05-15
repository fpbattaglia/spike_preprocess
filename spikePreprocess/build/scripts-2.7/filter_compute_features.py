#!/home/tara/miniconda3/envs/spikesorting/bin/python

from numpy import *
import sonIO
import neuralynx_io
from spikeFeatures import *
import spikePreprocess.computeSpikeFeatures as computeSpikeFeatures

import re


FEATURES_PER_TRODE = ( 'PC1: ', 'PC2: ', 'PC3: ')
FEATURES_TIM = ('Time in File (1/10000 s)', 'Time (1/10000 s)', 'File ID' )
FEATURES_MASK = {'Peak Index: ': (4, 14)}
N_PROBES = 4

def runKKFeaturesFromFiles(name, files, nProbes, featuresPerProbe, featuresTim, featuresMask, 
						   featuresMat=None, output='KK', extFeatures = None, realignPeaks = False,
						   shiftTimes = None, startTime = None, stopTime = None):
	"""run (as a simple script) all the steps necessary to prepare for running KK"""
	fetFeatures = []
	maskFeatures = {}
	if featuresMask is not None:
		fM = featuresMask.keys()
	else:
		fM = []
	
	for c in range(1, nProbes+1):
		for f in featuresPerProbe:
			fetFeatures.append(f + str(c))
			
			
	for f in fM:
		if f[-2:] == ': ':
			for c in range(1, nProbes+1):
				maskFeatures[f + str(c)] = featuresMask[f]
		else:
			maskFeatures[f] = featuresMask[f]
	
	
	if output == 'KK':
		fetFeatures.append('Time (Klusters Units)')
	fetFeatures = tuple(fetFeatures)
#	print fetFeatures	
	featuresToCompute = []
	featuresToCompute.extend(fetFeatures)
	if output == 'KK':
		featuresToCompute.extend(featuresTim)
	
	

	#TODO add mean waveforms
	
	# create data file handler objects 
	spikeFiles = []
	fileType = None
	if files[0][-3:] == 'smr':
		fileType = 'son'
	elif files[0][-3:] == 'ntt':
		fileType = 'nlynx'
	elif files[0][-3:] == 'spk':
		fileType = 'KK'
	if not fileType:
		raise ValueError('unrecognized file name '+ files)	
	if fileType == 'son':
		for fname in files:
			spikeFiles.append(sonIO.SONfile(fname))
	elif fileType == 'nlynx':
		for fname in files:
			spikeFiles.append(neuralynx_io.TTFile(fname))
	elif fileType == 'KK':
		for fname in files:
			spikeFiles.append(klusters_io.TTfile(fname))
			
	if fileType == 'son':
		adcMarkChans = spikeFiles[0].channelsOfKind('ADCMarkerChannel')
	elif fileType == 'nlynx' or fileType == 'KK':
		adcMarkChans = (0,)
		
	for chan in adcMarkChans:
		spkname = None
		if output == 'KK':
			fname = name + ".fet.1" #principal components file name needed
			# by klusters
			spkname = name + ".spk.1" # spike waveform file, needed by
			# klusters
			parname = name + ".xml"
			timname = name + ".tim.1"
		elif output == 'mat':
			import pymatfile
			matFname = name + '_features.mat'
			matf = pymatfile.pymatfile(matFname, 'w')
		if fileType == 'son':
			ttFiles = []
			for sp in spikeFiles:
				ttFiles.append(sonIO.TTFile(sp))
		else:
			ttFiles = spikeFiles
		computer = computeSpikeFeatures.FeatureComputer(ttFiles, spkname, params=None, interleave=120, shiftTimes = shiftTimes,
												startTime = startTime, stopTime = stopTime)
		computer.setMaskFromFeatures(maskFeatures)
		computer.computeFeatures(featuresToCompute, realignPeaks)
		features = computer.computedFeatures
		if output == 'KK':
			computeSpikeFeatures.writeFetFile(fname, features, fetFeatures)
			computeSpikeFeatures.writeTimFile(timname, features)
			computeSpikeFeatures.makeXMLParFile(spikeFiles, parname, len(featuresPerProbe), extFeatures)
		elif output == 'mat':
			for key in fetFeatures:
				value = features[key]
#				print shape(value)
				matVarname = re.sub(': ', '_', key)
				matVarname = re.sub(' ', '', matVarname)
#				print matVarname
				matf.writeItem(value.astype(float64), matVarname)

			matf.close()



if __name__ == '__main__':
	import sys, os, getopt

	def usage():
		print("""
USAGE:
filter_compute_features.py [options] file1 file2 file3 ...

options:
-n prefix (prefix for the output files)
-f feature (feature to compute, can be called many times)
-o output (output type can be 'mat' for matlab and 'KK' for KlustaKwik)
--peakindexmask (apply peak index mask)
--withshifttimes files are given in the format file1 tshift1 file2 tshift2 file3 tshift3 (for sake of joining files)
--samplingfrequency=sr set the sampling frequency to sr 
    (needed for Cheetah 4 in which sampling rate is not reported in the header)
--starttime=time process starting from this time 
--stoptime=time process up to this time 
""")

	
	if len(sys.argv) == 1:
		os.chdir('testData')
		runKKFeaturesFromFiles(name = 'tt1', files = ('tt1_test.ntt', ), nProbes=N_PROBES, 
			featuresPerProbe=FEATURES_PER_TRODE, featuresTim=FEATURES_TIM, featuresMask=FEATURES_MASK)
		sys.exit(0)
		
	runName = sys.argv[1]
	runFiles = sys.argv[2:]
	try:
		(opts, args) = getopt.gnu_getopt(sys.argv[1:], "n:f:o:m:", ['peakindexmask', 'samplingfrequency=',\
		 'realignpeaks', 'withshifttimes', 'starttime=', 'stoptime='])
	except getopt.GetoptError, err:
		# print help information and exit:
		print str(err) 
		usage()
		sys.exit(2)
	
	runFeatures = []
	runOutput = 'KK'
	runMask = None
	extFeatures = {}
	realignPeaks = False
	haveShiftTimes = False
	startTime = None
	stopTime = None
	for (key, value) in opts:
		if key == '-n':
			runName = value
		elif key == '-f':
#			print ('#'+ value +'#')
			if value == 'KKstd':
				runFeatures.extend(FEATURES_PER_TRODE)
			else:
				runFeatures.append(value)
		elif key == '-o':
			if value in ('KK', 'mat'):
				runOutput = value
			else:
				print('unrecognized output type ' + value )
			 	usage()
				sys.exit(2)
		elif key == '--peakindexmask':
			runMask = FEATURES_MASK
		elif key == '--samplingfrequency':
			extFeatures['SamplingFrequency'] = value
			extFeatures['WaveformLength'] = 32
		elif key == '--realignpeaks':
			realignPeaks = True
		elif key == '--withshifttimes':
			haveShiftTimes = True
		elif key == '--starttime':
			startTime = value
		elif key == '--stoptime':
			stopTime = value
		else:
			print('unrecognized option ' + key)
			usage()
			sys.exit(2)
	
	if not runFeatures:
		runFeatures = FEATURES_PER_TRODE
	
	if not (startTime is None and stopTime is None):
		runMask['Time (1/10000 s)'] = (startTime, stopTime)
	if haveShiftTimes:
		runFiles = []
		shiftTimes = []
		while args:
			runFiles.append(args.pop(0))
			shiftTimes.append(int(args.pop(0)))
		shiftTimes.pop(0)
		shiftTimes.append(0)
	else:
		runFiles = args
		shiftTimes = None
	print ('shiftTimes = ' + str(shiftTimes))
	print ('runName = ' + runName)
	print ('runFiles = ')
	print (runFiles)
	print ('runOutput = ' + runOutput)
	print ('runMask = ')
	print (runMask)
	print ('runFeatures = ')
	print (runFeatures)
	runKKFeaturesFromFiles(name = runName, files = runFiles, nProbes=N_PROBES, 
		featuresPerProbe=runFeatures, featuresTim=FEATURES_TIM, featuresMask=runMask,
		output=runOutput, extFeatures=extFeatures, realignPeaks=realignPeaks, shiftTimes = shiftTimes,
		startTime=startTime, stopTime=stopTime)
	
