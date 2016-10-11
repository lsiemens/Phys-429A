import numpy

class tektronix:
    def __init__(self):
        self.header_len = 17
        self.channels = 0
        self.points = 0
        self.yoffset = None
        self.scale = None
        self.xunits = ""
        self.yunits = ""
        
        self.channel_index_offset = 6

    def _load_header(self, fname):
        with open(fname, "r") as fin:
            for i, line in enumerate(fin):
                line = line.split(",")
                if i == 0:
                    self.channels = int(len(line)/6)
                    self.points = int(line[1]) - self.header_len
                    self.yoffset = [0.0 for _ in range(self.channels)]
                    self.scale = [[0.0 for _ in range(2)] for _ in range(self.channels)]
                if i == 7:
                    self.yunits = line[1]
                if i == 8:
                    for channel in range(self.channels):
                        self.scale[channel][1] = float(line[1])
                        line = line[self.channel_index_offset:]
                if i == 9:
                    for channel in range(self.channels):
                        self.yoffset[channel] = float(line[4])
                        line = line[self.channel_index_offset:]
                if i == 10:
                    self.xunits = line[1]
                if i == 11:
                    for channel in range(self.channels):
                        self.scale[channel][0] = float(line[1])
                        line = line[self.channel_index_offset:]

    def load_file(self, fname):
        self._load_header(fname)
        data = numpy.zeros(shape=(self.channels, 2, self.points))
        with open(fname, "r") as fin:
            for i, line in enumerate(fin):
                if i < self.header_len:
                    pass
                else:
                    line = line.split(",")
                    for channel in range(self.channels):
                        data[channel,0,i-self.header_len] = float(line[3])
                        data[channel,1,i-self.header_len] = float(line[4])
                        line = line[self.channel_index_offset:]
        return data

#from matplotlib import pyplot
#Tek = tektronix()
#data = Tek.load_file("../data/s1y5x5.1.csv")
#data = Tek.load_file("../data/s3y8.9x8.9.1.csv")

#data[2,1] = -0.5*(data[2,1] - 0.65)#set3

#for i in range(Tek.channels):
#    pyplot.plot(data[i,0], data[i,1])
#pyplot.show()

#xy = data[0,1]*data[1,1]
#xyh = xy*(1-0.0075*(data[0,1]**2 + data[1,1]**2))
#xyh = xy*(1+0.0015*(data[0,1]**2 + data[1,1]**2))
#pyplot.plot(data[0,0],-0.07*xy)
#pyplot.plot(data[0,0],-0.07*xyh)
#pyplot.plot(data[0,0],data[2,1])
#pyplot.show()

#pyplot.plot(-0.07*xy,data[2,1])
#pyplot.plot(-0.07*xyh,data[2,1])
#pyplot.show()
