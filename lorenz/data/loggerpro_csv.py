import numpy

class loggerpro:
    def __init__(self, channels=1):
        self.channels = channels + 1 #channels + time
    
    def load(self, fname):
        runs = 0
        data = None #[run][channel][iteration]
        with open(fname, "r") as fin:
            for i, line in enumerate(fin):
                if i == 0:
                    line = line.split(",")
                    runs = len(line)//self.channels
                    data = [[[] for channel in range(self.channels)] for run in range(runs)]
                else:
                    line = line.split(",")
                    for run in range(runs):
                        for channel in range(self.channels):
                            data[run][channel].append(float(line[channel + run*self.channels]))
            data = numpy.array(data)
        return data
