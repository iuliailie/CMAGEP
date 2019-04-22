
from __future__ import division
import sys, getopt
import random as ran
sys.path.append('./')
sys.path.append('./sympy')
#sys.stderr = open('Results/outPy/erorLog'+str(ran.randint(0,1000))+'.txt', 'w')#redirecting output to file
from sympy import *
import math

def main(argv):

  #print lis
  #print ar
   stro=''
   try:
      opts, args = getopt.getopt(argv,"hi:")
   except getopt.GetoptError:
      sys.exit(2)
   for opt, arg in opts:
       stro = arg

   lis=list('ABDGIJKMNOQRUVWXYZbdjkuvwyz');#list of symbols, order counts!!
   for i in range(0,len(lis)):
              lis[i]=symbols(lis[i])
   funce=stro
   for i in range(0,len(lis)):
        funce=funce.replace(str(lis[i]),"lis["+str(i)+"]")
   #del conjugate
   #print(funce)
   bla=''
   try:
       bla=N(eval(funce),2)#function gets simplified
      # bla=bla.evalf(10)
       con="conjugate"
       if con in str(bla) :
          # print '\n aici\n'
           bla=simplify(eval(funce))
   except:
        bla=nsimplify(eval(funce),rational=True)



   print str(bla)

if __name__ == "__main__":
   main(sys.argv[1:])
