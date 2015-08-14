import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import Tkinter, tkFileDialog
from matplotlib.lines import Line2D


configurations = {}

class Browser:
    def __init__(self):
	self.text = None
	self.selected = None
	self.text = ax.text(0.02, 0.98, 'selected: none',
				transform=ax.transAxes, va='top', fontsize=8)
    def onload(self, artist):
	if self.selected is not None:
	  self.selected.set_linewidth(0.7)
	artist.set_linewidth(1.2)
	self.selected = artist
	self.gid = artist.get_gid()
	self.update()
	
    def onpick(self, event):
	if isinstance(event.artist, Line2D):
	  self.selected.set_linewidth(0.7)
	  event.artist.set_linewidth(1.2)
	  self.selected = event.artist
	  self.gid = event.artist.get_gid()
	  self.update()
	
    def update(self):
        if self.gid is None: return

        dataind = self.gid
        
        print configurations[dataind]

        config = configurations[dataind]
	self.text.set_text("""
driver: %s
bufferSize: %s
numerOfThreads: %s
numberOfProcesses: %s
numberOfSpikeDetectorsPerThread: %s
numberOfMultimetersPerThread: %s
spikesPerDector: %s
samplingIntervalsOfMeter: %s
numberOfValuesWrittenByMeter: %s
deadTimeSpikeDetector: %s
deadTimeMultimeters: %s
deadTimeDeliver: %s"""
	  %(config[0],config[1],config[2],config[3],config[4],config[5],config[6],config[7],config[8],config[9],config[10],config[11]))

        fig.canvas.draw()

class Index:
    def __init__(self, browser):
	self.ind = 0
	self.browser = browser
    def next(self, event):
        root = Tkinter.Tk()
	root.withdraw()
	dir_path = tkFileDialog.askdirectory(parent=root,initialdir="/home/schumann/clusters/hambach.inm.kfa-juelich.de/test/nestio",title='Please select a directory')
	print dir_path


	filename_write = '%s/benchfile_write.csv'%dir_path
	filename_sleep = '%s/benchfile_sleep.csv'%dir_path
	filename_deliver = '%s/benchfile_deliver.csv'%dir_path
	filename_sync = '%s/benchfile_sync.csv'%dir_path

	reader_write=csv.reader(open(filename_write,"rb"),delimiter=';')
	x=list(reader_write)
	result_write=np.array(x).astype('float')

	reader_sleep=csv.reader(open(filename_sleep,"rb"),delimiter=';')
	x=list(reader_sleep)
	result_sleep=np.array(x).astype('float')

	reader_deliver=csv.reader(open(filename_deliver,"rb"),delimiter=';')
	x=list(reader_deliver)
	result_deliver=np.array(x).astype('float')

	reader_sync=csv.reader(open(filename_sync,"rb"),delimiter=';')
	x=list(reader_sync)
	result_sync=np.array(x).astype('float')

	Tw = result_write.T
	Tsl = result_sleep.T
	Tsy = result_sync.T
	Td = result_deliver.T
	
	gid = "dataset_%i"%self.ind
        artist_lst = ax.plot(np.mean(Tw-Tsl+Tsy, axis=1), picker=True, gid=gid) # have two versions of read.py only diff is /1000
	#artist_lst = ax.plot(np.mean(Tw-Tsl+Tsy, axis=1)/1000, picker=True, gid=gid) # validate order by comparing sum below to real execution time
        
        print "sum timings: ", np.sum(np.sum(Tw+Tsl+Tsy+Td, axis=1))
        
        reader_configuration = csv.reader(open('%s/config.csv'%dir_path,"rb"),delimiter=';')
        for row in reader_configuration:
	  configurations[gid] = row
	self.browser.onload(artist_lst[0])
        
        self.ind += 1 
        
        plt.draw()
        
        
fig, (ax) = plt.subplots(1,1)
plt.subplots_adjust(bottom=0.1)

ax.set_xlabel("iteration")
ax.set_ylabel("mean (threads and nodes) of timings (write+sync-sleep) in milliseconds")

browser = Browser()
fig.canvas.mpl_connect('pick_event', browser.onpick)

callback = Index(browser)
#axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.02, 0.07, 0.05])
bnext = Button(axnext, 'Load')
bnext.on_clicked(callback.next)
#bprev = Button(axprev, 'Previous')
#bprev.on_clicked(callback.prev)

plt.show()
