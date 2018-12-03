def scoreOrientation(digest, ladder):
    return (canDistinguish(sorted(set(digest.gDNA1Frags+digest.cDNA1Frags)),
        sorted(set(digest.gDNA2Frags+digest.cDNA2Frags)),ladder) +
    canDistinguish(sorted(set(digest.gDNA2Frags+digest.cDNA2Frags)),
        sorted(set(digest.gDNA1Frags+digest.cDNA1Frags)),ladder))/2

def scoreIdentity(digest, ladder):
    return (canDistinguish(sorted(set(digest.gDNA1Frags+digest.gDNA2Frags)), 
        sorted(set(digest.cDNA1Frags+digest.cDNA2Frags)),ladder) + 
    canDistinguish(sorted(set(digest.cDNA1Frags+digest.cDNA2Frags)),
        sorted(set(digest.gDNA1Frags+digest.gDNA2Frags)),ladder))/2

def closestTopLadder(frag,ladder):
    for i in range(len(ladder)):
        if ladder[i] < frag:
            return ladder[i-1]
    return ladder[0]

def closestBottomLadder(frag,ladder):
    for i in range(len(ladder)):
        if ladder[i] <= frag:
            return ladder[i]
    return ladder[-1]

def twoFragCompare(frag,other,ladder):
    topFragLadder = closestTopLadder(frag,ladder)
    bottomFragLadder = closestBottomLadder(frag,ladder)
    otherTopFragLadder = closestTopLadder(other,ladder)
    otherBottomFragLadder = closestBottomLadder(other,ladder)
    if topFragLadder == otherTopFragLadder and bottomFragLadder == otherBottomFragLadder:
        return (abs(other-frag)/(topFragLadder-bottomFragLadder+1))*(ladder.index(bottomFragLadder))
    else:
        topmost =  max(topFragLadder,otherTopFragLadder)
        bottommost = min(bottomFragLadder,otherBottomFragLadder)
        bandsBetween = ladder.index(bottommost)-ladder.index(topmost)
        return (abs(other-frag)/(topmost-bottommost+1))*(ladder.index(bottommost))+(ladder.index(bottommost))*bandsBetween

def canDistinguish(frags,otherFrags,ladder):
    score = 0
    for frag in frags:
        
        if frag in otherFrags or frag < ladder[-1]:
            continue
        closestTop = frag
        for x in otherFrags:
            if x > frag:
                closestTop = x
                break
        closestBottom = frag
        for x in otherFrags[::-1]:
            if x < frag:
                closestBottom = x
                break
        nearestFrag = closestBottom if ((closestBottom != frag) and (abs(closestTop-frag)>abs(closestBottom-frag))) else closestTop
        if otherFrags == []:
            nearestFrag = ladder[-1] if (abs(frag-ladder[-1])>abs(frag-ladder[0])) else ladder[0]
            score += twoFragCompare(frag,nearestFrag,ladder) * 2
        else:
            score += twoFragCompare(frag,nearestFrag,ladder)
    score /= (len(frags) if frags!=[] else 1)
    return score
    