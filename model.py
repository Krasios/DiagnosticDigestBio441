class Sample:

    def __init__(self,sequence,insertPos,intronPos,circular):
        self.circular = circular
        self.gDNA1 = OrientedSample(sequence,insertPos,intronPos,True,True)
        self.gDNA2 = OrientedSample(sequence,insertPos,intronPos,True,False)
        self.cDNA1 = OrientedSample(sequence,insertPos,intronPos,False,True)
        self.cDNA2 = OrientedSample(sequence,insertPos,intronPos,False,False)

class OrientedSample:

    def __init__(self,sequence,insertPos,intronPos,hasIntron,isForward):
        insertArea = sequence[insertPos[0]-1:insertPos[1]]
        comp = {'A':'T','C':'G','T':'A','G':'C'}
        if not hasIntron:
            beginEx = insertPos[0]-1
            insertArea = ""
            for intron in intronPos:
                insertArea += sequence[beginEx:intron[0]-1]
                beginEx = intron[1]
            insertArea += sequence[beginEx:insertPos[1]]
        if not isForward:
            insertArea = ''.join([comp[x] for x in insertArea[::-1]])
        introns = intronPos
        if hasIntron and not isForward:
            introns = [(insertPos[1]+insertPos[0]-x[1],insertPos[1]+insertPos[0]-x[0]) for x in intronPos]
        self.introns = introns
        self.sequence = sequence[:insertPos[0]-1] + insertArea + sequence[insertPos[1]:]
        self.hasIntron = hasIntron
        self.isForward = isForward

class RestrictionEnzyme:

    def __init__(self, sequence, cutPos, buffers):
        self.getPatterns(sequence,cutPos)
        self.buffers = buffers

    def getPatterns(self,sequence,cutPos):
        comp = {'A':'T','C':'G','T':'A','G':'C'}
        self.topPattern = self.getRegexPattern(sequence)
        self.topCut = len(sequence) + cutPos[0]
        self.bottomPattern = self.getRegexPattern(''.join([comp[b] if b in comp.keys() else b for b in sequence[::-1]]))
        self.bottomCut = 0 - cutPos[1]

    def getRegexPattern(self,sequence):
        slc = {
        'A':'A',
        'T':'T',
        'G':'G',
        'C':'C',
        'B':'[^A]',
        'V':'[^T]',
        'H':'[^G]',
        'D':'[^C]',
        'W':'[AT]',
        'R':'[AG]',
        'M':'[AC]',
        'K':'[TG]',
        'Y':'[TC]',
        'S':'[GC]',
        'N':'[ATGC]'}
        return '(?=('+''.join([slc[b] for b in sequence])+'))'