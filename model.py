class Sample:
    def __init__(self,sequence,insertPos,intronPos):
        self.gDNA1 = OrientedSample(sequence,insertPos,intronPos,True,True)
        self.gDNA2 = OrientedSample(sequence,insertPos,intronPos,True,False)
        self.cDNA1 = OrientedSample(sequence,insertPos,intronPos,False,True)
        self.cDNA2 = OrientedSample(sequence,insertPos,intronPos,False,False)

class OrientedSample:
    def __init__(self,sequence,insertPos,intronPos,hasIntron,isForward):
        insertArea = sequence[insertPos[0]-1:insertPos[1]]
        print(insertArea)
        if not hasIntron:
            insertArea = sequence[insertPos[0]-1:intronPos[0]-1] + sequence[intronPos[1]:insertPos[1]]
        if not isForward:
            insertArea = insertArea[::-1]
        self.sequence = sequence[:insertPos[0]-1] + insertArea + sequence[insertPos[1]:]
        self.hasIntron = hasIntron
        self.isForward = isForward
