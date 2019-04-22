

from __future__ import division

from sympy import *

import os
import sys
import math
import timeout

import subprocess
import sys
import random as ran

#
#sys.stdout = open('Results/outPy/outputPythonEU'+str(ran.randint(0,1000))+'.txt', 'w')#redirecting output to file
#sys.stderr = open('Results/outPy/erorLogEU'+str(ran.randint(0,1000))+'.txt', 'w')#redirecting output to file
'''py_function.py - Python source designed to '''
'''demonstrate the use of python embedding'''




def Teste(stra):
    print "a intrat in simpy"
    #print stra+

    T=symbols('T')
    stro=stra.replace("^","**")
    print "stro="+str(stro)
    try:
#        bla=eval(stro)
#        print "bla="+str(bla)
#        res=simplify(bla)
#        print "res="+str(res)
        strasd='python ./timeout.py 10 \'python ./interm.py -i "'+stro+'"\''#setting timeout for response of simplify
       # print strasd
        res = subprocess.check_output(strasd, shell=True)
        #res=os.popen(strasd).read()
       # print "res="+str(res)
        stri=str(res[:-1]).replace("**","^")
        stri=stri.replace(" ","")
        print "res="+str(stri)
        con="conjugate"
        if con in stri :
            print "a iesit din simpy conj"
            return '!'+stra
              #  print stri

        print "a iesit din simpy"
        print "\n"
        return stri
    except:
        print "Unexpected error: ", sys.exc_info()[0]

        print "a iesit din simpy"
        print "\n"
        return '!'+stra










