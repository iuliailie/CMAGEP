from __future__ import division


import sys
import math
import re as REg
from sympy import *
#import numpy as np

def Evs(optimizedParam,function,terminale):
    print "a intrat in evs"
    lis=list(terminale)
    for i in range(0,len(terminale)):
            lis[i]=symbols(lis[i])
            #print lis
            #print ar

  #funSct='['+funce+' for M in X]'

    prez=optimizedParam

    print "param list "+ str(prez)
    #fun='P[0]*X**P[1]+P[2]*cos(P[3]*Y)'
    fun=function.replace("^","**")
    for i in range(0,len(terminale)):
        fun=fun.replace(str(lis[i]),"lis["+str(i)+"]")

    for i in range(0,len(prez)):
       fun=fun.replace("**P["+str(i)+"]","**int(P["+str(i)+"])",)

    for i in range(0,len(prez)):
        fun=fun.replace('P['+str(i)+']',"round(prez["+str(i)+"],3)")
    print fun
    #print "\n"
    fun=fun.replace(" ","")
    #print fun
    ka=''
    try:
        ka=N(eval(fun),2)
        con="conjugate"
        if con in str(ka) :
            ka=simplify(eval(fun))
    except:
        ka=nsimplify(eval(fun),rational=True)


    eq=str(ka).replace("**","^")

    fileObject = open('./Results/RezultateOptimizari.5.txt','a')
    fileObject.write(eq)
    fileObject.write('\n')
    fileObject.close();



#    for i in range(0,len(lis)):
#        eq=eq.replace(str(lis[i]),"1*"+str(lis[i]))  #adaug 1 la var care nu au nimic in fata
#
#    print eq

#    conj="conjugate"
#    if conj in r :
#        print "a iesit din evs conj"
#        return None

   # Consts=REg.findall(r"[-+]?\d*\.\d+|\d+",eq.replace(" ",""))
   # print Consts
    #Consts=[float(x) for x in Consts]
    #print Consts

    print "a iesit din evs"
    return (eq,optimizedParam)

