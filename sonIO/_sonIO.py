#!/usr/bin/env python


## \file sonIO.py
""" \brief This modules contains data structures to read spike2 SON files
   generating several types of output
"""
## mostly read stuff

__revision__ = "0.1"


from numpy import *


import struct, re, sys

class TTFile(object):
	"""\brief This class implements the interface for one tetrode"""
	def __init__(self, name, chan):
		"""name can be a SONfile object, or a file name, in the latter case the corresponding SONfile 
		object is created."""
		
		if isinstance(name, SONfile):
			self.sonFile = name
		else:
			self.sonFile = SONfile(name)
		self.chan = chan
	
	def firstTime(self):
		return self.sonFile.firstTime()
	
	def lastTime(self):
		return self.sonFile.lastTime()
	
	def readSpikesByBlocks(self, firstBlock = 0, lastBlock = -1, timeUnit = 0.0001):
		"""\brief reada chunk of spikes by block
		
		This function return blocks of data (that is, blocks in the
		CED file) from the TT file. This is done in order to
		save memory. data is raw in the sense that it is kept in the
		original data format (int32 for time, int16 for data) data are
		returned as a bidimensional matrix, and the waveform from
		multiple probes are not sorted out. The only processing
		performed is conversion of times to desired time unit
		\param firstBlock (\c int)  the first block in the file to be
		read (defaults to 0 = first block)
		\param lastBlock (\c int) one past the last block to be read
		(defaults to -1 = last block)
		\param timeUnit (\c float) the desired output time unit (in
		seconds)
		\return (t, marks, wv) where t is a (nPointsx1) array, marks
		is a (nPoints x 4) array, wv is a (nPointsxnWvPoints)"""
		
		(t, marks, wv) = self.sonFile.getADCMarkerChannelRawBlocks(self.chan, firstBlock, lastBlock, timeUnit)
		return (t, wv)
	
	def genTTData(self, spikeFiles=[], nBlocksAtOnce=50, nProbes=4, interleave = 0, yieldShiftedTime=False,
		yieldMarkers=False):
		""" generator yielding SON data in blocks to avoid excessive memory
		comsumption, if a non-null file handle is supplied as spikeHandle,
		will write the spike waveform data in a format suitable for klusters in that handles
		"""

		son = self.sonFile
		chan = self.chan
		offset = 0
		timeShift = 0
		ns = 0
		if not (len(spikeFiles) > 0 and self is spikeFiles[0]):
			spikeFiles.insert(0, son)
		for son in spikeFiles:
			
			info = son.info(chan)
			nBlocks = info['nBlocks']
			nWvPoints = info['nWvPoints']
			nProbePoints = nWvPoints / nProbes


			for i in range(0,nBlocks, nBlocksAtOnce):
				(t, marks, wvRaw) = son.getADCMarkerChannelRawBlocks(chan,
																  firstBlock=i,
																  lastBlock =
																  i+nBlocksAtOnce)
				nPoints = shape(wvRaw)[0]
				wv = reshape(wvRaw, [nPoints, nProbePoints, nProbes])
				wv = transpose(wv, [0,2,1])
				wv = wv.astype(float32)
				if yieldShiftedTime:
					ts = t+ timeShift
					fId = ns * ones(shape(ts))
					yield (t, wv, offset, offset+len(t), ts, fId, wvRaw)
				elif yieldMarkers:
					yield (t, wv, offset, offset+len(t), wvRaw, marks)
				else:
					yield (t, wv, offset, offset+len(t), wvRaw)
				offset += len(t)

			timeShift += (max(son.lastTimes) +interleave) * 10000
			ns += 1
	
	def info(self):
		d = self.sonFile.info(self.chan)
		d['nProbes'] = 4
		return d
	


class SONfile(object):
	"""\brief This class implentes input for the moment, and noe day
	output from/to a SON file
	This class is an implementation of a reader for CED son files in
	python, providing functionalities for
	- wavemark files (spike data)
	- ADC channels (EEG)
	- Marker channels
	""" 
	
	class __SONchannel:
		"""data structure for one channel in the SON file  """
		def __init__(self):
			(self.delSize,
			self.nextDelBlock,
			self.firstblock,
			self.lastblock,
			self.blocks,
			self.nExtra,
			self.preTrig,
			self.free0,
			self.phySz,
			self.maxData) = [0 for i in range(10)]
			self.filename = ''
			self.comment = []
			(self.maxChanTime,
			 self.lChanDvd,
			 self.phyChan) = [0, 0, 0]
			(self.idealRate,
			 self.kind,
			 self.pad) = [0, 0, 0]
			self.title = ''
			(self.scale,
			 self.offset,
			 self.units,
			 self.divide,
			 self.interleave,
			 self.min,
			 self.max,
			 self.initLow,
			 self.nextLow) = [0 for i in range(9)] 
		
	
	class __BlockHeader:
		"""data structure for the block header in the SON file. The block
		    header will be read and kept in memory upon file opening
		"""
		def __init__(self, nblocks):
			"""allocate space for the block header data
			"""
			self.predBlock = zeros(nblocks, int32) # offset of
			# previous block in file
			self.succBlock = zeros(nblocks, int32) # offset of next
			# block for channel in file
			self.startTime = zeros(nblocks, int32) # time in time
			#	 units for first timestamp in block
			self.endTime = zeros(nblocks, int32) # time in time
			# units for end timestamp in 
			self.items = zeros(nblocks, int16) # number of items
			# in the block
			self.nBlocks = nblocks
			self.empty = True
	
	def __init__(self, name):
		"""\brief opens a SON file and creates a SON object		
		takes a filename as an argument, try to open for read.
		raises IOError on failure
		\param name: the name of the file to open
		"""	
		
		self.channelKind = {}
		self.channelKind[1] = 'ADCChannel'
		self.channelKind[2] = 'EventChannel'
		self.channelKind[3] = 'EventChannel'
		self.channelKind[4] = 'EventChannel'
		self.channelKind[5] = 'MarkerChannel'
		self.channelKind[6] = 'ADCMarkerChannel'
		self.channelKind[7] = 'RealMarkerChannel'
		self.channelKind[8] = 'TextMarkerChannel'
		self.channelKind[9] = 'RealWaveChannel'		
		self.__filename = name
		self.__fh = file(name, 'rb')
		self.__processHeader()
		self.__channelInfo = []
			
			
		for i in range(self.__channels):
			self.__channelInfo.append(self.__readChannelInfo(i))
			
			
		self.____BlockHeader = []
		for i in range(self.__channels):
			self.____BlockHeader.append(self.__readBlockHeader(i))
		
		
		self.__chanList = self.__makeChannelList()
		
		self.lastTimes = []
		for i in range(self.__channels):
			self.lastTimes.append(self.lastTime(i))
		
		self.lastTimes = array(self.lastTimes)
	
	def ticksToSecondsFactor(self):
		""" \brief return the factor that has to be multipled by the time in
		 file to get seconds
		"""
		
		if self.__systemID >= 6:
			conv = self.__usPerTime * self.__dTimeBase
		else:
			conv = self.__usPerTime * 1e-6

		return float(conv)
	
	def firstTime(self, chan):
		return 0
		
		
	def lastTime(self, chan):
		"""\brief returns the last time in seconds for channel chan
		"""
		conv = self.ticksToSecondsFactor()
		if len(self.____BlockHeader[chan].endTime):
			t = conv * self.____BlockHeader[chan].endTime[-1]
		else:
			t = 0
		return t
	
	def sampleInterval(self, chan):
		""" \brief returns the claimed sample interval
		 (inverse of sampling
		rate) in microseconds
		"""

		info = self.__channelInfo[chan]
#		if not (info.kind in (1, 5, 6, 7, 9)):
#			raise IOError('Invalid Channel Type')
		if self.__systemID < 6:
			if hasattr(info, 'divide'):
				_sampleInterval = info.divide * \
								  self.__usPerTime * self.__timePerADC * 1e-6
			else:
				_sampleInterval = None
		else:
			_sampleInterval = info.lChanDvd*self.__usPerTime*self.__dTimeBase

		return _sampleInterval
	
	def getADCMarkerChannelRawBlocks(self, chan, firstBlock = 0, lastBlock = -1, timeUnit = 0.0001):
		"""\brief yields block of raw data from a ADCMarker (wavemark, or
		spikes) channel

		This function return blocks of data (that is, blocks in the
		CED file) from a ADCMarker channel. This is done in order to
		save memory. data is raw in the sense that it is kept in the
		original data format (int32 for time, int16 for data) data are
		returned as a bidimensional matrix, and the waveform from
		multiple probes are not sorted out. The only processing
		performed is conversion of times to desired time unit
		\param chan (\c int) the number of the channel to read
		\param firstBlock (\c int)  the first block in the file to be
		read (defaults to 0 = first block)
		\param lastBlock (\c int) one past the last block to be read
		(defaults to -1 = last block)
		\param timeUnit (\c float) the desired output time unit (in
		seconds)
		\return (t, marks, wv) where t is a (nPointsx1) array, marks
		is a (nPoints x 4) array, wv is a (nPointsxnWvPoints) 
		"""
		fh = self.__fh
		header = self.____BlockHeader[chan]
		info = self.__channelInfo[chan]


		nBlocks = info.blocks
		nWvPoints = info.nExtra/2
		nMarkers = 4
		headerSize = 20
		conv = float(self.ticksToSecondsFactor())
		scale = float(info.scale) / 6553.6
		offset = info.offset

		if lastBlock == -1 or lastBlock > nBlocks:
			lastBlock = nBlocks

		nMarks = sum(header.items[firstBlock:lastBlock].astype(int32))
		t = zeros(nMarks, float64)
		tout = zeros(nMarks, int32)
		wv = zeros([nMarks, nWvPoints], int16)
		marks = zeros([nMarks, nMarkers], int8)

		c = 0

		if firstBlock == 0:
			off = info.firstblock+headerSize
		else:
			off = header.succBlock[firstBlock-1]+headerSize

		for b in range(firstBlock,lastBlock):
			fh.seek(off)
			off = header.succBlock[b]+headerSize
			for i in range(header.items[b]):
				(tc,) = struct.unpack('=i', fh.read(4))
				t[c] = tc*conv
				marks[c, :] = fromstring(fh.read(4), int8)
				wv[c, :] = fromstring(fh.read(nWvPoints*2), int16)
				c += 1


		tout = t/timeUnit
		tout = tout.astype(int32)


		return (tout, marks, wv)
	
	def getADCMarkerChannel(self, chan, nProbes=4):
		"""\brief get an ADC Marker Channel output

		Returns all the spikes relative to a certain channel, as
		float64 (double) arrays, times are output in seconds, waveform
		data are returned in Volts
		\param chan the number of desired channel
		\param nProbes the number of probes (defaults to 4, as for a
		tetrode)
		\return - t a (nPoints x 1) column double array of times, in seconds
		- marks (nPoints x 4) int8 array of markers, for example
		associated to spike sorting
		- wv a (nPoints x nProbes x NProbePoints) array containing the
		waveforms, waveforms relative to the different probes are
		disentangled
		.
		"""
		fh = self.__fh
		header = self.____BlockHeader[chan]
		info = self.__channelInfo[chan]
		nMarks = sum(header.items.astype(int32))
		nBlocks = info.blocks
		nWvPoints = info.nExtra/2
		nProbePoints = nWvPoints/nProbes
		nMarkers = 4
		headerSize = 20
		conv = float(self.ticksToSecondsFactor())
		scale = float(info.scale) / 6553.6
		offset = info.offset

		t = zeros(nMarks, float64)
		wv = zeros([nMarks, nProbes, nProbePoints], float64)
		marks = zeros([nMarks, nMarkers], int8)

		c = 0
		off = info.firstblock+headerSize
		for b in range(nBlocks):
			fh.seek(off)
			off = header.succBlock[b]+headerSize
			for i in range(header.items[b]):
				(tc,) = struct.unpack('=i', fh.read(4))

				t[c] = tc*conv
				marks[c, :] = fromstring(fh.read(4), int8)
				wv[c, :, :] = transpose(reshape(
					fromstring(fh.read(nWvPoints*2), int16),
				[nProbePoints, nProbes]))
				c += 1

		wv = wv * scale + offset
		return (t, marks, wv)
	
	def getADCChannel(self, chan):
		"""\brief get an ADC  Channel output

		returns all data relative to an ADC (EEG) channel
		\param chan the number of desired channel
		\return - t a (nPoints x 1) double array of times, in seconds
		- wv a (nPoints x 1) double array of AD readouts, in Volts
		.
		"""
		fh = self.__fh
		header = self.____BlockHeader[chan]
		info = self.__channelInfo[chan]

		nMarkers = sum(header.items.astype(int32))
		nBlocks = info.blocks
		headerSize = 20
		conv = float(self.ticksToSecondsFactor())
		sampleInterval = self.sampleInterval(chan)
		scale = float(info.scale) / 6553.6
		offset = info.offset
		t = zeros(nMarkers, float64)
		wv = zeros(nMarkers, float64)

		c = 0
		off = info.firstblock+headerSize
		for b in range(nBlocks):
			fh.seek(off)
			off = header.succBlock[b]+headerSize
			firstTime = conv * header.startTime[b]
			lastTime = conv * header.endTime[b]
			items = header.items[b]
			realSampleInterval =  float(lastTime-firstTime) /  (items-1)


			t[c:(c+items)] = arange(firstTime-realSampleInterval, 
					   lastTime-1e-10,
					   realSampleInterval,
					   float64) + realSampleInterval
			wv[c:(c+items)] = fromstring(fh.read(2*items), int16)
			c += items

		wv = wv*scale +offset

		return (t, wv)
	
	def getMarkerChannel(self, chan):
		"""\brief get a Marker channel data

		returns all data in a marker channel
		\param chan (\c int) the number of desired channel
		\return - t a (nPoints x 1 ) double array of times, in seconds
		- marks (nPoints x 4) int8 array of markers
		.

		"""
		fh = self.__fh
		header = self.____BlockHeader[chan]
		info = self.__channelInfo[chan]
		nMarks = sum(header.items.astype(int32))
		nBlocks = info.blocks
		nMarkers = 4
		headerSize = 20
		conv = float(self.ticksToSecondsFactor())

		t = zeros(nMarks, float64)
		marks = zeros([nMarks, nMarkers],  int8)
		c = 0
		off = info.firstblock+headerSize
		for b in range(nBlocks):
			fh.seek(off)
			off = header.succBlock[b]+headerSize

			for i in range(header.items[b]):
				t[c] = conv * fromstring(fh.read(4), int32)
				hh = fromstring(fh.read(4), int8)
				marks[c, :] = hh
				c += 1

		return (t, marks)
	
	def getEventChannel(self, chan):
		"""\brief get a Event channel data

		returns all data in a marker channel
		\param chan (\c int) the number of desired channel
		\return - t a (nPoints x 1 ) double array of times, in seconds
		.

		"""
		fh = self.__fh
		header = self.____BlockHeader[chan]
		info = self.__channelInfo[chan]
		nMarks = sum(header.items.astype(int32))
		nBlocks = info.blocks
		headerSize = 20
		conv = float(self.ticksToSecondsFactor())

		t = zeros(nMarks, float64)
		c = 0
		off = info.firstblock+headerSize
		for b in range(nBlocks):
			fh.seek(off)
			off = header.succBlock[b]+headerSize

			for i in range(header.items[b]):
				t[c] = conv * fromstring(fh.read(4), int32)
				c += 1

		return t
	
	def isValidChannel(self, chan):
		"""
		\brief whether a channel is valid/contains data
		
		\param chan the number of channel to query
		\return \a True if channel as valid data in it, \a False otherwise 
		"""
		
		info = self.__channelInfo[chan]
		header = self.____BlockHeader[chan]

		if info.kind == 0:
			return False

		if header.empty:
			return False

		return True
	
	def channelsOfKind(self, kind):
		""" \brief returns list of channels in file of a certain kind

		\param kind: one of the channelKind identifiers:
		- ADCChannel
		- EventChannel
		- MarkerChannel
		- ADCMarkerChannel
		- RealMarkerChannel
		- TextMarkerChannel
		- RealWaveChannel
		.
		see also the SON fiel format specification from CED
		\return a list of channel numbers (all the ones of that type
		contained in the file) 
		"""
		return filter(lambda x: self.__chanList[x].kind == kind,
					  self.__chanList.keys())
	
	def showHeaderInfo(self):
		"""\brief prints out all the info contained in the file header
		"""
		print ("""system ID %d
copyright = %s
creator = %s
		""" % (self.__systemID, self.__copyright, self.__creator))

		print ("""usPerTime = %d
timePerADC = %d
filestate = %d
firstdata = %d
channels = %d
chansize = %d
extraData = %d
bufferSize = %d
osFormat = %d
maxFTime = %d
dTimeBase = %g
timeDateDetail = %d
timeDateYear = %d
		""" % (self.__usPerTime, # short
		 self.__timePerADC, # short
		 self.__filestate, # short
		 self.__firstdata, # int
		 self.__channels, # short
		 self.__chansize, # short
		 self.__extraData, #short
		 self.__bufferSize, # short
		 self.__osFormat, # short
		 self.__maxFTime, # int
		 self.__dTimeBase, #double
		 self.__timeDateDetail, # uchar
		 self.__timeDateYear # short
		 ))

		print "file comments:\n"
		print self.__fileComments	
	
	def info(self, chan):
		d = {}
		d['kind'] = self.__channelInfo[chan].kind
		d['peakSampleIndex'] = 8
		d.update(self.__chanList[chan].__dict__)
		d['SamplingFrequency'] = 1000000. / self.sampleInterval(chan)
		scale = 2**15 * float(self.__channelInfo[chan].scale) / 6553.6
		d['voltageRange'] = scale
		return d
	
	def __makeChannelList(self):
		"""
		build the information summaries that we output, really for use
		in the constructor only 
		""" 
		__chanList = {}
		conv = self.ticksToSecondsFactor()
		for i in range(self.__channels):
			if self.isValidChannel(i):
				scInfo = self.__ShortChannelInfo()
				info  = self.__channelInfo[i]
				header = self.____BlockHeader[i]
				scInfo.title = info.title
				scInfo.firstTime = 0
				scInfo.endTime = -1
				
				if len(header.startTime) > 0:
					scInfo.firstTime = conv * header.startTime[0]
					scInfo.endTime = conv * header.endTime[-1]
				scInfo.nPoints = sum(header.items.astype(int32))
				scInfo.kind = self.channelKind[info.kind]
				scInfo.nBlocks = info.blocks
				if scInfo.kind == 'ADCMarkerChannel':
					scInfo.nWvPoints = info.nExtra/2
				if scInfo.kind in  ('ADCMarkerChannel', 'ADCChannel'):
					scInfo.sampleInterval = self.sampleInterval(i)
				
				__chanList[i] = scInfo
		return __chanList
	
	def __processHeader(self):
		""" internal function: reads the SON header """
		fh = self.__fh
		fh.seek(0)
		
		(self.__systemID,) = struct.unpack('h', fh.read(2))
		self.__copyright = fh.read(10)
		self.__creator = fh.read(8)
		#print "struct size = %d\n" % struct.calcsize('=hhhihhhhhidBh')
		(self.__usPerTime, # short
		 self.__timePerADC, # short
		 self.__filestate, # short
		 self.__firstdata, # int
		 self.__channels, # short
		 self.__chansize, # short
		 self.__extraData, #short
		 self.__bufferSize, # short
		 self.__osFormat, # short
		 self.__maxFTime, # int
		 self.__dTimeBase, #double
		 self.__timeDateDetail, # uchar
		 self.__timeDateYear # short
		 ) = struct.unpack('=hhhihhhhhidBh', fh.read(35))

		self.__headerPad = fh.read(52)
		self.__channels += 1
		point = fh.tell()
		fc = []


		for i in range(5):
			(bytes,) = struct.unpack('B', fh.read(1))
			fc.append(fh.read(bytes))
			point = point + 80
			fh.seek(point)
			
		self.__fileComments = tuple(fc)
		
	
	def __readChannelInfo(self, chan):
		"""
		reads the info about a channel
		"""
		
		if chan > self.__channels:
			raise ValueError("channel %d does not exist in this file" %
						  (chan,))
		# Offset due to file header and preceding channel headers
		base = 512+(140*(chan-1));   

		fh = self.__fh
		fh.seek(base)

		ch = self.__SONchannel()
		ch.filename = self.__filename

		(ch.delSize,
		ch.nextDelBlock,
		ch.firstblock,
		ch.lastblock,
		ch.blocks,
		ch.nExtra,
		ch.preTrig,
		ch.free0,
		ch.phySz,
		ch.maxData) = struct.unpack('=hiiihhhhhh', fh.read(26))

		(bytes,) = struct.unpack('=B', fh.read(1))
		point = fh.tell()
		ch.comment = fh.read(bytes)
		fh.seek(point+71)

		(ch.maxChanTime,
		 ch.lChanDvd,
		 ch.phyChan, bytes) = struct.unpack('=iihB',
		 fh.read(struct.calcsize('=iihB')))

		point = fh.tell()
		ch.title = fh.read(bytes)
		fh.seek(point+9)


		(ch.idealRate,
		 ch.kind,
		 ch.pad) = struct.unpack('=fBb',
		 fh.read(struct.calcsize('=fBb')))



		if ch.kind in (1, 6):
			(ch.scale,
			 ch.offset,
			 bytes) = struct.unpack('=ffB',
									fh.read(struct.calcsize('=ffB')))
			point = fh.tell()
			ch.units = fh.read(bytes)
			fh.seek(point+5)

			if self.__systemID < 6:
				(ch.divide,) = struct.unpack('=h',
				fh.read(struct.calcsize('=h')))
			else:
				(ch.interleave,) = struct.unpack('=h',
				fh.read(struct.calcsize('=h')))
		elif ch.kind in (7,9):
			(ch.min,
			 ch.max,
			 bytes) = struct.unpack('=ffB',
									fh.read(struct.calcsize('=ffB')))
			point = fh.tell()
			ch.units = fh.read(bytes)
			fh.seek(point+5)

			if self.__systemID < 6:
				(ch.divide,) = struct.unpack('=h',
				fh.read(struct.calcsize('=h')))
			else:
				(ch.interleave,) = struct.unpack('=h',
				fh.read(struct.calcsize('=h')))
		
		elif ch.kind in (4,):
			(ch.initLow, ch.nextLow) = struct.unpack('=BB',
			fh.read(struct.calcsize('=BB')))


		return ch
	
	def __readBlockHeader(self, chan):
		"""internal function reads the header of a block """
		chanInfo = self.__channelInfo[chan]
		nBlocks = chanInfo.blocks
		header = self.__BlockHeader(nBlocks)
		fh = self.__fh

		# get to the first block, if block non existing, return an
		# empty header
		nextBlock = chanInfo.firstblock
		# if they are marker channels we'll consider them non empty
		# even if they don't have anything in them
		if nextBlock == -1 and chanInfo.kind != 5:
			return header

		# start loading blocks
		for i in range(nBlocks):
			if nextBlock == -1:
				raise IOError("""Blocks ended before intended for
channel """ + str(chan))
			
			
			fh.seek(nextBlock)
			(header.predBlock[i],
			 header.succBlock[i],
			 header.startTime[i],
			 header.endTime[i],
			 chanNumber,
			 header.items[i]) = struct.unpack('=iiiihh',
			 fh.read(struct.calcsize('=iiiihh')))
			if chanNumber != chan:
				raise IOError('block header with wrong channel number found \
				for channel ' + str(chan))
			nextBlock = header.succBlock[i]
			
			
		header.empty = False
		return header
	
	class __ShortChannelInfo:
		"""class to summarize channel information """
		def __init__(self):
			pass
		
	







#####################################


	
if __name__ == '__main__':


	
	s = SONfile(sys.argv[1])
	#s.showHeaderInfo()
	
