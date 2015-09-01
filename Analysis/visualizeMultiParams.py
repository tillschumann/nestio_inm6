import re
import numpy as np
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) <2):
  sys.exit();


npzfile = np.load(sys.argv[1])

#numberOfSpikes= npzfile['arr_0']
#deadtime_sd_update= npzfile['arr_1']
#numberOfSDPerThread= npzfile['arr_2']
#deliver_duration= npzfile['arr_3']

numberOfValues= npzfile['arr_4']
numberOfMultiCallsPerInterval= npzfile['arr_5']
deadtime_multi_update= npzfile['arr_6']
numberOfHandles= npzfile['arr_7']
numberOfMultisPerThread= npzfile['arr_8']

print numberOfValues
print numberOfMultiCallsPerInterval
print deadtime_multi_update
print numberOfHandles
print numberOfMultisPerThread

print "========== RESULT =========="
print "numberOfValues: mean=%f variance=%f" % (numberOfValues.mean(),numberOfValues.var())
print "multi dead time: mean=%f variance=%f" % (deadtime_multi_update.mean(),deadtime_multi_update.var())
print "numberOfMultiCallsPerInterval: mean=%f variance=%f" % (numberOfMultiCallsPerInterval.mean(),numberOfMultiCallsPerInterval.var())
print "numberOfHandles: mean=%f variance=%f" % (numberOfHandles.mean(),numberOfHandles.var())
print "numberOfMultisPerThread: mean=%f variance=%f" % (numberOfMultisPerThread.mean(),numberOfMultisPerThread.var())

print "========== BOUNDARIES ======"
print "numberOfValues: min=%f max=%f" % (numberOfValues.min(),numberOfValues.max())
print "multi dead time: min=%f max=%f" % (deadtime_multi_update.min(),deadtime_multi_update.max())
print "numberOfMultiCallsPerInterval: min=%f max=%f" % (numberOfMultiCallsPerInterval.min(),numberOfMultiCallsPerInterval.max())
print "numberOfHandles: min=%f max=%f" % (numberOfHandles.min(),numberOfHandles.max())
print "numberOfMultisPerThread: min=%f max=%f" % (numberOfMultisPerThread.min(),numberOfMultisPerThread.max())

f, (ax1,ax2,ax3,ax4,ax5)=plt.subplots(5,1)


ax1.plot(numberOfValues)
ax1.set_title("number of written values")

ax2.plot(deadtime_multi_update)
ax2.set_title("deattime between multimeters update calls")

ax3.plot(numberOfMultiCallsPerInterval)
ax3.set_title("number of record calls per interval")

ax4.plot(numberOfHandles)
ax4.set_title("number of handler calls")

ax5.plot(numberOfMultisPerThread)
ax5.set_title("number of multis per thread")

f.subplots_adjust(hspace=0.5)
plt.show()
