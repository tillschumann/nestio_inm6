import re
import numpy as np
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) <2):
  sys.exit();


npzfile = np.load(sys.argv[1])
numberOfSpikes= npzfile['arr_0']
deadtime_sd_update= npzfile['arr_1']
numberOfSDPerThread= npzfile['arr_2']
deliver_duration= npzfile['arr_3']

#numberOfValues= npzfile['arr_4']
#numberOfMultiCallsPerInterval= npzfile['arr_5']
#deadtime_multi_update= npzfile['arr_6']
#numberOfHandles= npzfile['arr_7']
#numberOfMultisPerThread= npzfile['arr_8']

print "========== RESULT =========="
print "spike: mean=%f variance=%f" % (numberOfSpikes.mean(),numberOfSpikes.var())
print "spikedetector dead time: mean=%f variance=%f" % (deadtime_sd_update.mean(),deadtime_sd_update.var())
print "numberOfSpikedetectorsPerThread: mean=%f variance=%f" % (numberOfSDPerThread.mean(),numberOfSDPerThread.var())
print "deliver_duration: mean=%f variance=%f" % (deliver_duration.mean(),deliver_duration.var())

print "========== BOUNDARIES ======"
print "spike: min=%f max=%f" % (numberOfSpikes.min(),numberOfSpikes.max())
print "spikedetector dead time: min=%f max=%f" % (deadtime_sd_update.min(),deadtime_sd_update.max())
print "numberOfSpikedetectorsPerThread: min=%f max=%f" % (numberOfSDPerThread.min(),numberOfSDPerThread.max())
print "deliver_duration: min=%f max=%f" % (deliver_duration.min(),deliver_duration.max())

f, (ax1,ax2,ax3,ax4)=plt.subplots(4,1)


ax1.plot(numberOfSpikes)
ax1.set_title("spikes per update call")

ax2.plot(deadtime_sd_update)
ax2.set_title("deattime between spikedetector update calls")

ax3.plot(numberOfSDPerThread)
ax3.set_title("number of spikedetectors per thread")

ax4.plot(deliver_duration)
ax4.set_title("duration of deliver call")

f.subplots_adjust(hspace=0.5)
plt.show()

#plt.hist(deadtime_sd_update, bins=1000)
#plt.title("deadtime between spikedetector update calls")
#plt.show()
