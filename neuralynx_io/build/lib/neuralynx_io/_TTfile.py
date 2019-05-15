#!/usr/bin/env python

## \file _TTfile.py
""" \brief This module contains data structures to read neuralynx
tetrode (.ntt extension), as generated from newer NT-based versions of cheetah
"""

__revision__ = "0.1"

from numpy import *

import re, os, os.path, time, sys

class TTFile(object):
	"""\brief This class implements input for the moment, and one day
	output from/to a Neuralynx tetrode (TT) file
	
	"""
	def __init__(self, name, mode='r', header=None):
		"""
		\brief opens a TT file and creates a TT fiel object
		\param name: the (path) name of the file to be opened
		\param mode: can be 'r' for read (default), and 'w' for write
		\param header: a header for a file to be written
		"""
		
		if mode != 'r' and mode != 'w':
			raise ValueError('invalid value for argument "mode"')
		
		self.mode = mode
		if mode == 'r':
			# read mode constructor
			# determine endianness of the machine
			self.isBigEndian = dtype(int32).str[0] == '>'
			self.filename = name
			self.fh = file(name, 'r')
			
			# determines if there is an old style header
			self.headerInfo = {}
			firstLine = self.fh.readline(80)
			if firstLine[:13] == "%%BEGINHEADER": # there's a header
				raise IOError("It seems like an old SUN file, which is not implemented yet")
##				 while firstLine[:11] != "%%ENDHEADER":
##					 self.fh.readline(80)
##				 self.offset = self.fh.tell()
			elif firstLine[:35] == "######## Neuralynx Data File Header" :
				self.__processHeader()
			else:
				self.offset = 0
			
			if self.headerInfo.has_key('RecordSize'):
				self.recordSize = int(self.headerInfo['RecordSize'])
			else:
				self.recordSize = 8 + 10 * 4 + 2 * 128
			self.__makeChannelInfo()
			
			self.fh.seek(self.offset)
		elif mode == 'w':
			raise NotImplementedError("Write mode not implemented yet.")
			# we can only write new style neuralynx file
			self.filename = name
			self.fh = file(name, 'w')
			if header:
				self.fh.write(header)
			else:
				self.fh.write("######## Neuralynx Data File Header\n")
				self.fh.write('## File Name ' + os.path.abspath(self.filename))
			self._lastTime = 0
	
	def firstTime(self):
		if self.mode == 'r':
			return self.__chanList[0].firstTime
		else:
			raise ValueError("it's a write mode file")
	
	def lastTime(self):
		if self.mode == 'r':
			return self.__chanList[0].lastTime
		else:
			raise ValueError("it's a write mode file")
	
	def writeSpikes(self, tStamps, wv):
		"""\brief writes a block of spikes
		   \param t: an int 64 numpy array of timestamps
		   \param wv: an int15 numpy array of waveforms 
		"""
		t = tStamps
		if self.mode != 'w':
			raise ValueError('Not a write mode file')
		if t[0] < self._lastTime:
			raise ValueError('block contains timestamps preceding the last timestamp currently in file')
		self._lastTime = t[-1]
		nSpikesToWrite = len(t)
		for i in range(nSpikesToWrite):
			self.fh.write(t[i].tostring())
			self.fh.write(wv[i, :].tostring())
	
	def close(self):
		"""\brief closes a file
		"""
		self.fh.close()
	
	def readSpikes(self, firstSpike, lastSpike, timeUnit = 0.0001):
		"""
		\brief reads a range of spikes
		\param firstSpike ordinal of first requested spike
		\param lastSpike ordinal of last requested spike (-1 for the last spike in the file)
		\param timeUnit desired time unit for the timestamps
		\return (t, wv), where t is a vector if int64 timestamps and wv is a (nChannels*Samples) x nSpikes
		array of waveform values
		"""
		if self.mode != 'r':
			raise ValueError('Not a read mode file')
		
		tConv = 1000000. * timeUnit
		nChannels = self.__chanList[0].nProbes
		nSamplesPerSpike = self.__chanList[0].waveformLength
		
		if lastSpike > self.__chanList[0].nPoints - 1 or lastSpike == -1:
			lastSpike  = self.__chanList[0].nPoints - 1

		if firstSpike > self.__chanList[0].nPoints - 1:
			firstSpike  = self.__chanList[0].nPoints - 1

		nSpikesToRead = lastSpike-firstSpike+1


		t = zeros(nSpikesToRead, int64)
		wv = zeros([nSpikesToRead, (nChannels*nSamplesPerSpike)], int16)

		self.fh.seek(self.offset+firstSpike*self.recordSize)

		for i in range(nSpikesToRead):
			t[i] = fromstring(self.fh.read(8), int64)
			dw = self.fh.read(40)
			wv[i,:] =fromstring(self.fh.read((nChannels*nSamplesPerSpike*2)), int16)

		t /= tConv
		return (t, wv)
	
	def genTTData(self, spikeFiles = [], nSpikesAtOnce=50000, nProbes=4, interleave = 0, yieldShiftedTime=False,
	 	shiftTimes=None):
		"""
		\brief generator yielding a block of spikes with timestamps and waveforms
		\param spikeFiles: a list of further spike files (including or not the one the TTfile object is based on), to 
			be concatenated
		\param nSpikesAtOnce (default 50000): the number of spikes to be read at each call of the generator
		\param nProbes (default 4): the number of probes (4 for a tetrode)
		\param interleave: a interval in time between to consecutive concatenated files (so that the concatenation is 
		easily recognizable)
		\param yieldShiftedTime: if True, the shifted time for concatenated file, together with a file Id will also be 
		yielded
		"""
		if shiftTimes:
			shiftTimes = list(shiftTimes) # work on a copy so that we don't destroy information for later use 
		print ("genTTdata init: shiftTimes: ", str(shiftTimes))
		timeShift = 0
		offsetTot = 0
		ns = 0
		chan = 0
		nWvPoints = self.__chanList[chan].nWvPoints
		nProbePoints = nWvPoints / nProbes
		spFiles = spikeFiles
		if shiftTimes:
			totShiftTimes = 0
		
		if not (len(spFiles) > 0 and self is spFiles[0]):
			spFiles.insert(0, self)
		for tt in spFiles:
			offset = 0
			nReads = int(ceil(float(tt.__chanList[chan].nPoints) / nSpikesAtOnce))
			for i in range(nReads):
				(t, wvRaw) = tt.readSpikes(offset, offset+nSpikesAtOnce-1)
				nPoints = shape(wvRaw)[0]			

				wv = reshape(wvRaw, [nPoints, nProbePoints, nProbes])
				wv = transpose(wv, [0, 2, 1])
				wv = wv.astype(float32)
				if yieldShiftedTime:
					ts = t+ timeShift
					fId = ns * ones(shape(ts))
					yield (t, wv, offsetTot, offsetTot+len(t), ts, fId, wvRaw)
				else:
					yield (t, wv, offsetTot, offsetTot+len(t), wvRaw)
				offset += len(t)
				offsetTot += len(t)
			if yieldShiftedTime:
				print ("genTTdata: yst: shiftTimes: "+ str(shiftTimes))
				if shiftTimes:
					timeShift += shiftTimes.pop(0)
				elif len(spFiles) > 1:
					raise ValueError("no ShiftTimes")
			ns += 1
	
	def info(self, chan=0):
		d = self.headerInfo
		d.update(self.__chanList[0].__dict__)
		d['nProbes'] = 4 # hack that will only work for tetrode configurations 
		if not d.has_key('WaveformLength'):
			d['WaveformLength'] = 32 # default value needed for the chetah 4 files 
		return d 

	
	def __processTimeLinesHeader(self, headerLine, openClose):
		mtch = re.match(r'## Time ' + openClose + r'.*\(m/d/y\): (.*)  \(h:m:s\.ms\) (.*)', headerLine)
		if not mtch:
			mtch = re.match(r'## Time ' + openClose + r'.*\(m/d/y\): (.*)  At Time: (.*)', headerLine)
		if not mtch:
			raise IOError('Unexpected line in file header')
		dateTimeStr = ' '.join((mtch.group(1), mtch.group(2)))
		mtch = re.match('(.*)\..*', dateTimeStr)
		if mtch:
			dateTimeStr = mtch.group(1)
		date = time.strptime(dateTimeStr, '%m/%d/%Y %H:%M:%S')
		return date
	
	def __processHeader(self):
		"""
		internal function, process the file header, based on the
		header format as in Cheetah 5.0.2
		"""
		self.offset = 16384
		
		hinfo = {}
		# the original filename
		headerLine = self.fh.readline(80).rstrip()
		print(headerLine)
		mtch = re.match(r'## File Name[:\s]+(.*)$', headerLine)
		hinfo['OrigFileName'] = mtch.group(1)
		if not mtch:
			raise IOError('Unexpected line in file header')
		mtch = re.match(r'.*\\(.*)\\.*', hinfo['OrigFileName'])
		hinfo['Session'] = mtch.group(1)
		
		# time opened
		headerLine = self.fh.readline(80).rstrip()
		date = self.__processTimeLinesHeader(headerLine, 'Opened')
		hinfo['dateTimeOpened'] = date
		
		# time closed
		headerLine = self.fh.readline(80).rstrip()
		try:
			date = self.__processTimeLinesHeader(headerLine, 'Closed')
		except IOError:
			date = ''
		hinfo['dateTimeClosed'] = date
		
		# further lines
		headerLine = self.fh.readline(160).rstrip().lstrip(' 	')
		mtch = re.match('\s*(.*)', headerLine)
		headerLine = mtch.group(1)
		while (len(headerLine) == 0) or (headerLine[0] == '-')  :
			print headerLine
			if len(headerLine) == 0:
				pass
			elif headerLine[:8] == '-Feature':
				mtch = re.match(r'-Feature\s+(.*)\s+(.*)', headerLine)
				hinfo['Feature_' + mtch.group(1)] = mtch.group(2)
			elif headerLine[0] == '-':
				mtch = re.match(r'\s*-(\S*)\s+(.*)', headerLine)
				hinfo[mtch.group(1)] = mtch.group(2)
			headerLine = self.fh.readline(160).rstrip().lstrip()
			#mtch = re.match('\s*(.*)', headerLine)
			#headerLine = mtch.group(1)
		print 'Last Line #' + headerLine + '#'
		print hinfo
		self.headerInfo = hinfo
	
	def __makeChannelInfo(self):

		__chanList = []
		scInfo = self.__ShortChannelInfo()

		fileLength = os.stat(self.filename).st_size
		scInfo.nPoints = (fileLength-self.offset) / self.recordSize
		scInfo.kind = 'TT'
		scInfo.nBlocks = 1 # no meaning for blocks in Neuralynx files

		if self.headerInfo.has_key('NumADChannels'):
			scInfo.nProbes = int(self.headerInfo['NumADChannels'])
		else:
			scInfo.nProbes = 4
		if self.headerInfo.has_key('WaveformLength'):
			scInfo.waveformLength = int(self.headerInfo['WaveformLength'])
		else:
			scInfo.waveformLength = 32
		scInfo.nWvPoints = scInfo.nProbes * scInfo.waveformLength

		if self.headerInfo.has_key('AcqEntName'):
			scInfo.title = self.headerInfo['AcqEntName']

		else:
			rt = os.path.basename(self.filename)
			(rt, ext) = os.path.splitext(rt)
			scInfo.title = rt

		if self.headerInfo.has_key('SamplingFrequency'):
			scInfo.sampleInterval = 1. / float(self.headerInfo['SamplingFrequency'])
		else:
			scInfo.sampleInterval = 1./32000


		self.fh.seek(self.offset)
		scInfo.firstTime = fromstring(self.fh.read(8), int64)

		self.fh.seek(-self.recordSize, 2)
		scInfo.lastTime = fromstring(self.fh.read(8), int64)
		self.lastTimes = scInfo.lastTime

		self.fh.seek(self.offset)

		__chanList.append(scInfo)
		self.__chanList = __chanList
	
	class __ShortChannelInfo:
		"""class to summarize channel information 
		"""
		def __init__(self):
			pass
	
	


if __name__ == '__main__':
	
	s = TTFile(sys.argv[1])
	(t, wv) = s.readSpikes(0, 1000)