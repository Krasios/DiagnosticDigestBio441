from model import *
from digest_model import *
from fragment_ranking import *
from itertools import combinations

def runSuggester(seqFile,enzymeFile,insertPositions,intronPositions,circular,ladder):
    sFile = open(seqFile, 'r')
    ladders = [[10000,8000,6000,5000,4000,3000,2500,2000,1500,1000,750,500,253],[1500,1000,900,800,700,600,500,400,300,200,100]]
    ladder = ladders[ladder]
    seq = sFile.read()
    intronPositions = [(int(x[0]),int(x[1]))for x in intronPositions]
    sample = Sample(str(seq),insertPositions,intronPositions,circular)
    lines = [line.rstrip('\n').split(',') for line in open(enzymeFile,'r')]
    enzymes = {}
    for line in lines:
        enzymes[line[0]] = RestrictionEnzyme(str(line[1]),[int(line[2]),int(line[3])],'')
    digests = {}
    for enzymeName,enzyme in enzymes.items():
        digests[enzymeName] = SingleDigest(sample,enzyme)
    filteredDigests = {k:v for k,v in digests.items() if v.hasCut}
    enzymePairs = combinations(filteredDigests.keys(),2)
    for pair in enzymePairs:
        filteredDigests[pair[0]+" && "+pair[1]] = DoubleDigest(filteredDigests[pair[0]],filteredDigests[pair[1]],sample)  
    scoreId = {k:scoreIdentity(v,ladder) for k,v in filteredDigests.items()}
    scoreOri = {k:scoreOrientation(v,ladder) for k,v in filteredDigests.items()}
    topId = []
    if intronPositions:
        topId = sorted([(v,k) for k,v in scoreId.items()])[::-1][:10]
    topOri = sorted([(v,k) for k,v in scoreOri.items()])[::-1][:10]
    topIdString = '\n'.join([x[1] for x in topId]) if id else "No Introns Specified in Sample"
    topOriString = '\n'.join([x[1] for x in topOri])
    return {'identity': topId, 'orientation': topOri}