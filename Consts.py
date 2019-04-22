import re

def replaceConstants(strd):
    #print strd
    #strd=strd+'+P'
    tep=re.findall(r'[\d\.]+(?:[eE][+-]?\d+)?',strd)    
    try:
        tepp=[float(x) for x in tep]
        #tep=map(float, )
        l=len(tep)
       # print l
        if  l > 0:
           sub="P";
           stringForReplacing=strd.replace(tep[0],sub,1)
           tempo=strd
           for i in range(0,l):
               # print i
                sub="P";
               # stringForReplacing tep[i]
                stringForReplacing=strd.replace(tep[i],sub,1)
                strd=stringForReplacing
                tempo=str(stringForReplacing)
                #print str(mghj)
           return (str(tempo),tepp)
    except:
        return (str(strd),'')
    print(strd)
    return (str(strd),'')


def replaceConstants2(strd):

    final_list = []
    for elem in strd.split():
        try:
            final_list.append(float(elem))
        except ValueError:
            pass
    bla=strd.split()
    print str(bla)
    l=len(final_list)
    print str(final_list)
   # print l
    if  l > 0:
       sub="P";
       tempString=strd.replace(final_list[0],sub,1)
       StringToReturn=tempString
       for i in range(0,l):
           # print i
            sub="P";
           # print tep[i]
            mghj=strd.replace(final_list[i],sub,1)
            strd=mghj
            StringToReturn=str(mghj)
       return str(StringToReturn)
    return str(strd)
