import FitnessFunctionsForCma
import cma


def Testarea(X,Y,func,nbParam,ff):
    #T=symbols('T')
    #print "a intrat in cma"
    funcu=func.replace("^","**")
    funca=funcu.replace("P","prez")
    rez=[0]
    if nbParam >1:
        if ff==1 :
            rez=cma.fmin(FitnessFunctionsForCma.ME, nbParam*[0], 0.5,args=(X,Y,funcu))
        else :
            rez=cma.fmin(FitnessFunctionsForCma.ME, nbParam*[0], 0.5,args=(X,Y,funcu))
        prez=rez[0].tolist()
        rezult=prez
    #print "a iesit din cma"
        return (prez,rez[1],funca)
    else :
        return (rez,1002,funca)




