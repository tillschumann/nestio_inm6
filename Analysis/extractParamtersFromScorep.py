import re
import numpy as np
import matplotlib.pyplot as plt
import time


p_sd_update = re.compile("([A-Z]*)[ \t]*([0-9]*)[ \t]*([0-9]*).*spike_detector::update.*")
p_sd_spike = re.compile("([A-Z]*)[ \t]*([0-9]*)[ \t]*([0-9]*).*record_event_internal_flush.*")
p_update = re.compile("([A-Z]*)[ \t]*([0-9]*)[ \t]*([0-9]*).*threaded_update_openmp_UPDATE.*")
p_deliver = re.compile("([A-Z]*)[ \t]*([0-9]*)[ \t]*([0-9]*).*threaded_update_openmp_DELIVER.*")

p_multi_update = re.compile("([A-Z]*)[ \t]*([0-9]*)[ \t]*([0-9]*).*Multimeter::update.*")
p_multi_handles = re.compile("([A-Z]*)[ \t]*([0-9]*)[ \t]*([0-9]*).*Multimeter::handle.*")
p_multi_printvalues = re.compile("([A-Z]*)[ \t]*([0-9]*)[ \t]*([0-9]*).*Multimeter::print_value_.*")
p_multi_printvalue = re.compile("([A-Z]*)[ \t]*([0-9]*)[ \t]*([0-9]*).*RecordingDevice::print_value.*")

fobj_in = open("/scratch/schumann/traces_small.txt")


Tbegin_sd_block = {}
Tend_sd_block = {}
spikes_sd_block = {}
multi_values = {}
multi_inverval = {}
multi_handles = {}

last_spike = {}
last_multi = {}
deliver_start = {}

numberOfMultimeters = {}
numberOfSpikedetectors = {}

###
##
## separate nodes and threads: 
##	node = int(LOC[-32:],2)
##	thread = int(LOC[:-32],2)
deliver_duration = []
deadtime_sd_update = []
deadtime_multi_update = []
numberOfSpikes = []
numberOfValues = []
numberOfHandles = []
numberOfMultiCallsPerInterval = []
numberOfSDPerThread = []
numberOfMultisPerThread = []



#deliver_duration = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#    deliver_duration[i][j] = []
#deadtime_sd_update = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#		deadtime_sd_update[i][j] = []
#deadtime_multi_update = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#		deadtime_multi_update[i][j] = []
#numberOfSpikes = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#    numberOfSpikes[i][j] = []
#numberOfValues = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#    numberOfValues[i][j] = []
#numberOfHandles = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#    numberOfHandles[i][j] = []
#numberOfMultiCallsPerInterval = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#    numberOfMultiCallsPerInterval[i][j] = []
#numberOfSDPerThread = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#    numberOfSDPerThread[i][j] = []
#numberOfMultisPerThread = np.empty((32,8), dtype=list)
#for i in range(32):
#  for j in range(8):
#    numberOfMultisPerThread[i][j] = []


#entries = sum(1 for line in fobj_in)
entries = 438271456.0
i=0.0
c=0

start_time = time.time()

for line in fobj_in:
  if c>100000:
    cur_time = time.time()
    print "%f percent\t %i minutes to go" %(i/entries*100, (cur_time-start_time)/i*(entries-i)/60)
    c=0
  else:
    c+=1
  i+=1 	
  
  #
  #
  #
  if ("threaded_update_openmp_DELIVER" in line):
    m_deliver = p_deliver.match(line)
    if m_deliver:
      KEY = m_deliver.group(1)
      LOC = int(m_deliver.group(2))
      TIME = long(m_deliver.group(3))
      if KEY == "ENTER":
	deliver_start[LOC] = TIME
      else:
	node = int(bin(LOC+pow(2,38))[-32:],2)
	thread = int(bin(LOC+pow(2,38))[3:-32],2)
	#print "%i %i" %(node,thread)
	#deliver_duration[node][thread].append(TIME-deliver_start[LOC])
	deliver_duration.append(TIME-deliver_start[LOC])
  #
  #
  #
  if ("threaded_update_openmp_UPDATE" in line):
    m_update = p_update.match(line)
    if m_update:
      KEY = m_update.group(1)
      LOC = int(m_update.group(2))
      TIME = long(m_update.group(3))
      if KEY == "ENTER":
	last_spike[LOC]=-1
	last_multi[LOC]=-1
	numberOfSpikedetectors[LOC]=0
	numberOfMultimeters[LOC]=0
      else:
	node = int(bin(LOC+pow(2,38))[-32:],2)
	thread = int(bin(LOC+pow(2,38))[3:-32],2)
	#numberOfSDPerThread[node][thread].append(numberOfSpikedetectors[LOC])
	#numberOfMultisPerThread[node][thread].append(numberOfMultimeters[LOC])
	numberOfSDPerThread.append(numberOfSpikedetectors[LOC])
	numberOfMultisPerThread.append(numberOfMultimeters[LOC])
      
    #
    #
    #
  if ("spike_detector::update" in line):
    m_sd_update = p_sd_update.match(line)
    if m_sd_update:
      KEY = m_sd_update.group(1)
      LOC = int(m_sd_update.group(2))
      TIME = long(m_sd_update.group(3))
      node = int(bin(LOC+pow(2,38))[-32:],2)
      thread = int(bin(LOC+pow(2,38))[3:-32],2)
      if KEY == "ENTER":
	numberOfSpikedetectors[LOC]=numberOfSpikedetectors[LOC]+1
	#Tbegin_sd_block[LOC] = TIME
	spikes_sd_block[LOC] = 0
	if last_spike[LOC]>0:
	  #deadtime_sd_update[node][thread].append(TIME-last_spike[LOC])
	  deadtime_sd_update.append(TIME-last_spike[LOC])
	last_spike[LOC] = TIME
      else:
	#Tend_sd_block[LOC] = TIME
	#numberOfSpikes[node][thread].append(spikes_sd_block[LOC])
	numberOfSpikes.append(spikes_sd_block[LOC])
	#print "location=%s:\t numberOfSpikes=%i"%(LOC,spikes_sd_block[LOC])
	
	
    #
    #
    #
  if ("record_event_internal_flush" in line):
    m_sd_spike = p_sd_spike.match(line)
    if m_sd_spike:
      KEY = m_sd_spike.group(1)
      LOC = int(m_sd_spike.group(2))
      TIME = long(m_sd_spike.group(3))
      if KEY == "ENTER":
	spikes_sd_block[LOC]=spikes_sd_block[LOC]+1
	  	
    #
    #
    #
  if ("Multimeter::update" in line):
    m_multi_update = p_multi_update.match(line)
    if m_multi_update:
      KEY = m_multi_update.group(1)
      LOC = int(m_multi_update.group(2))
      TIME = long(m_multi_update.group(3))
      node = int(bin(LOC+pow(2,38))[-32:],2)
      thread = int(bin(LOC+pow(2,38))[3:-32],2)
      if KEY == "ENTER":
	numberOfMultimeters[LOC]=numberOfMultimeters[LOC]+1
	#Tbegin_sd_block[LOC] = TIME
	multi_handles[LOC] = 0
	if last_multi[LOC]>0:
	  #deadtime_multi_update[node][thread].append(TIME-last_multi[LOC])
	  deadtime_multi_update.append(TIME-last_multi[LOC])
	last_multi[LOC] = TIME
      else:
	#Tend_sd_block[LOC] = TIME
	#numberOfHandles[node][thread].append(multi_handles[LOC])
	numberOfHandles.append(multi_handles[LOC])
	#print "location=%s:\t multi_handles=%i"%(LOC,multi_handles[LOC])

	    	
    #
    #
    #
  if ("Multimeter::handle" in line):
    m_multi_handles = p_multi_handles.match(line)
    if m_multi_handles:
      KEY = m_multi_handles.group(1)
      LOC = int(m_multi_handles.group(2))
      TIME = long(m_multi_handles.group(3))
      if KEY == "ENTER":
	multi_handles[LOC] = multi_handles[LOC] +1
	multi_inverval[LOC] = 0
      else:
        node = int(bin(LOC+pow(2,38))[-32:],2)
	thread = int(bin(LOC+pow(2,38))[3:-32],2)
	#numberOfMultiCallsPerInterval[node][thread].append(multi_inverval[LOC])
	numberOfMultiCallsPerInterval.append(multi_inverval[LOC])
		  
		  
	      
	  
	
    #
    #
    #
  if ("Multimeter::print_value_" in line):
    m_multi_printvalues = p_multi_printvalues.match(line)
    if m_multi_printvalues:
      KEY = m_multi_printvalues.group(1)
      LOC = int(m_multi_printvalues.group(2))
      TIME = long(m_multi_printvalues.group(3))
      if KEY == "ENTER":
	multi_inverval[LOC] = multi_inverval[LOC] +1
	multi_values[LOC] = 0
      else:
        node = int(bin(LOC+pow(2,38))[-32:],2)
	thread = int(bin(LOC+pow(2,38))[3:-32],2)
	#numberOfValues[node][thread].append(multi_values[LOC])
	numberOfValues.append(multi_values[LOC])

		
		  
	
  #
  #
  #
  if ("RecordingDevice::print_value" in line):
    m_multi_printvalue = p_multi_printvalue.match(line)
    if m_multi_printvalue:
      KEY = m_multi_printvalue.group(1)
      LOC = int(m_multi_printvalue.group(2))
      TIME = long(m_multi_printvalue.group(3))
      if KEY == "ENTER":
	multi_values[LOC] = multi_values[LOC] +1
      
      
fobj_in.close()

numberOfSpikes = np.array(numberOfSpikes)
deadtime_sd_update = np.array(deadtime_sd_update)
numberOfSDPerThread = np.array(numberOfSDPerThread)
deliver_duration = np.array(deliver_duration)

numberOfMultisPerThread = np.array(numberOfMultisPerThread)
numberOfValues = np.array(numberOfValues)
numberOfHandles = np.array(numberOfHandles)
numberOfMultiCallsPerInterval = np.array(numberOfMultiCallsPerInterval)
deadtime_multi_update = np.array(deadtime_multi_update)


#print "========== RESULT =========="
#print "spike: mean=%f variance=%f" % (numberOfSpikes[0][0].mean(),numberOfSpikes.var())
#print "spikedetector dead time: mean=%f variance=%f" % (deadtime_sd_update.mean(),deadtime_sd_update.var())
#print "numberOfSpikedetectorsPerThread: mean=%f variance=%f" % (numberOfSDPerThread.mean(),numberOfSDPerThread.var())
#print "deliver_duration: mean=%f variance=%f" % (deliver_duration.mean(),deliver_duration.var())
#print "numberOfValues: mean=%f variance=%f" % (numberOfValues.mean(),numberOfValues.var())
#print "numberOfMultiCallsPerInterval: mean=%f variance=%f" % (numberOfMultiCallsPerInterval.mean(),numberOfMultiCallsPerInterval.var())
#print "deadtime_multi_update: mean=%f variance=%f" % (deadtime_multi_update.mean(),deadtime_multi_update.var())
#print "numberOfMultisPerThread: mean=%f variance=%f" % (numberOfMultisPerThread.mean(),numberOfMultisPerThread.var())

#print "======= CHECK ========="
#print "spike: max=%f min=%f" % (numberOfSpikes.max(),numberOfSpikes.min())
#print "spikedetector dead time: max=%f min=%f" % (deadtime_sd_update.max(),deadtime_sd_update.min())
#print "numberOfSpikedetectorsPerThread: max=%f min=%f" % (numberOfSDPerThread.max(),numberOfSDPerThread.min())
#print "deliver_duration: max=%f min=%f" % (deliver_duration.max(),deliver_duration.min())
#print "numberOfValues: max=%f min=%f" % (numberOfValues.max(),numberOfValues.min())
#print "numberOfMultiCallsPerInterval: max=%f min=%f" % (numberOfMultiCallsPerInterval.max(),numberOfMultiCallsPerInterval.min())
#print "deadtime_multi_update: max=%f min=%f" % (deadtime_multi_update.max(),deadtime_multi_update.min())
#print "numberOfMultisPerThread: max=%f min=%f" % (numberOfMultisPerThread.max(),numberOfMultisPerThread.min())

np.savez('/scratch/schumann/parameters_small.npz', numberOfSpikes,deadtime_sd_update,numberOfSDPerThread,deliver_duration,numberOfValues,numberOfMultiCallsPerInterval,deadtime_multi_update,numberOfHandles,numberOfMultisPerThread)
