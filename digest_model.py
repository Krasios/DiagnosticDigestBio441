import re

class SingleDigest:

    def __init__(self,sample,enzyme):
        self.gDNA1Cuts = self.processSample(sample.gDNA1,enzyme)
        self.gDNA2Cuts = self.processSample(sample.gDNA2,enzyme)
        self.cDNA1Cuts = self.processSample(sample.cDNA1,enzyme)
        self.cDNA2Cuts = self.processSample(sample.cDNA2,enzyme)
        self.gDNA1Frags = self.processCuts(sample.circular,sample.gDNA1.sequence,self.gDNA1Cuts)
        self.gDNA2Frags = self.processCuts(sample.circular,sample.gDNA2.sequence,self.gDNA2Cuts)
        self.cDNA1Frags = self.processCuts(sample.circular,sample.cDNA1.sequence,self.cDNA1Cuts)
        self.cDNA2Frags = self.processCuts(sample.circular,sample.cDNA2.sequence,self.cDNA2Cuts)

    def processSample(self,sample,enzyme):
        topMatches = re.finditer(enzyme.topPattern,sample.sequence)
        topCuts = [(i.start(1)+enzyme.topCut) for i in topMatches]
        bottomMatches = re.finditer(enzyme.bottomPattern,sample.sequence)
        bottomCuts = [(i.start(1)+enzyme.bottomCut) for i in bottomMatches]
        topCuts.extend([x for x in bottomCuts if x not in topCuts])
        return sorted(topCuts)

    def processCuts(self,circular,sequence,cuts):
        frag = []
        for i in range(1,len(cuts)):
            frag.append(cuts[i]-cuts[i-1])
        if circular:
            frag.append(len(sequence)-cuts[-1]+cuts[0])
        else:
            frag.append(cuts[0])
            frag.append(len(sequence)-cuts[-1])
        return sorted(frag)
