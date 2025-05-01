class Window():
    def __init__(self,chrm,start,stop):
        self.chrm = chrm
        self.start = start
        self.stop = stop

    def overlap(self,other):
        if self.chrm == other.chrm:
            if self.start >= other.start and self.start <= other.stop:
                # print("1")
                # print(self.start,other.start)
                # print(self.start,other.stop)
                return True
                
            elif self.start<other.start and  self.stop > other.start:
                #print("2")
                return True
            else:
                return False
        else:
            return False

    def __str__(self):
        return f"{self.chrm}:{self.start}-{self.stop}"

    def __repr__(self):
        return self.__str__()

def parse_file(f):
    
    wnds = []
    with open(f) as inp:
        header = True
        for line in inp:
            if header:
                header = False
                continue
            else:
                data = line.strip().split(":")
                chrm = data[0]
                start = int(data[1].split("-")[0])
                stop = int(data[1].split("-")[1])
                wnds.append(Window(chrm,start,stop))   

    return wnds


def parse_cpgs(f):
    
    wnds = []
    with open(f) as inp:
        header = True
        for line in inp:
            if header:
                header = False
                continue
            else:
                data = line.strip().split("\t")
            
                chrm = data[0]
                start = int(data[1])
                stop = int(data[1])
                wnds.append(Window(chrm,start,stop))   

    return wnds

setA = parse_cpgs("dmers_positions.bed") # positions of the DMRs or set of positions A that need to be checked

setB = parse_file("pcadat_positions.txt") ## positions of the DMRs or set of positions B that need to be checked


for wp in setA:
    for wd in setB:
        if wd.overlap(wp):
            print(wp,wd)
            