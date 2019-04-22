import cma
import FitnessFunctionsForCma
import random


def CallFunctie(X,Y,func,nbParam,ff,terminalStr,OptimizedBeforeConsts):
   # print('ff e'+str(ff))
    #T=symbols('T')
    #print "a intrat in cma"
#    print X
#    print"\n"
#    print Y
#    print terminalStr
#    print "\n"
#    print func
#    func='P[0]*W+P[1]*T+P[2]*P[3]**(P[4]*T)*P[5]**(P[6]*Y-P[7])'
  #  func='P[0]*t-P[1]*W+P[2]*log(P[3]*W)+P[4]*z/t+P[5]*t/(P[6]*y-P[7]*x +P[8])-P[9]'
    #func='P[0]*P[1]**(P[2]*T+P[3])+P[4]*W';
    #func='P[0]*log(P[1]*W)+(P[2]*A/((t/Z)+P[3]))+(P[4]*x/log(((P[5]*A*t)+P[6]*y)))'
    #func='P[0]*W+P[1]*log(P[2]*X-P[3]*t-P[4]*W+P[5]*Z)+(P[6]*G(P[7]*G*(P[8]*X-P[9]*Z))/P[10]*Z-P[11]'
#log(W)+(A/((t/g)+9))+(x/log(((A*t)+y)))
    #nbParam=12;
   # print OptimizedBeforeConsts
    #print('orf=',str(orf))
    print str(OptimizedBeforeConsts)
    funcu=func.replace("^","**")
    funca=funcu.replace("P","prez")
    print(str(nbParam))
    print "functia din cma= "+str(funca)
    rez=[[0],]
    rezTemp=[[0],]
    fit=1e50;
    conj="conjugate"
    if conj in funca :
        print "a iesit din evs conj"
        return (rez,2e50,funca,0)
    if nbParam >1:
        if ff==1 :
            print "RMSE \n"
            for i in range(0,1):
                if len(OptimizedBeforeConsts)>0:
                    PARAMS=OptimizedBeforeConsts
                else:
                    PARAMS=nbParam*[random.uniform(0, 10)]
                rezTemp=cma.fmin(FitnessFunctionsForCma.RMSE, PARAMS, 0.05,args=(X,Y,funcu,nbParam,terminalStr))
                if rezTemp[1] >1e50 :
                    break
                if rezTemp[1]<fit :
                    print rezTemp
                    fit=rezTemp[1]
                    rez=rezTemp

        elif ff==4:

            for i in range(0,1):
              #  OptimizedBeforeConsts=[ -0.15,0.23,0.678,3,-0.026,0.15]
                PARAMS=[]
                if len(OptimizedBeforeConsts)==nbParam:
                    print str(OptimizedBeforeConsts)
                    PARAMS=OptimizedBeforeConsts
                    if len(PARAMS)==1:
                        PARAMS.append(random.uniform(0, 1))

                  #  print PARAMS
                else:
                    for j in range(0,nbParam):
                        PARAMS.append(random.uniform(0, 1))
                        #print str(PARAMS[j])


                #print str(PARAMS)

                rezTemp=cma.fmin(FitnessFunctionsForCma.ME, PARAMS, 0.5,args=(X,Y,funcu,nbParam,terminalStr))
                #print rezTemp[1]
                if rezTemp[1]>1e50:
                    break
                if rezTemp[1]<fit :
                    fit=rezTemp[1]
                    rez=rezTemp
#                fileObject = open('./Results/RezultateOptimizareSmoth60.txt','a')
#                fileObject.write(','.join(map(str,rezTemp[0])))
#                fileObject.write('\n')
#                fileObject.close();


            prez=[0]
            try:
                #print rez
                prez=rez[0].tolist()
                me=FitnessFunctionsForCma.ME
                call=me.called
    #            print '\n test '
                print str(call)
                print '\n'+str(fit)

            except:
                return (PARAMS,2e50,funca,0)
        #print "a iesit din cma"
            return (prez,rez[1],funca,call)

        elif ff==5:
            print('a intrat la aic')

            for i in range(0,1):
              #  OptimizedBeforeConsts=[ -0.15,0.23,0.678,3,-0.026,0.15]
                PARAMS=[]
                if len(OptimizedBeforeConsts)>0:
                    print str(OptimizedBeforeConsts)
                    PARAMS=OptimizedBeforeConsts
                  #  print PARAMS
                else:
                    for j in range(0,nbParam):
                        PARAMS.append(random.uniform(0, 1))
                        #print str(PARAMS[j])


                #print str(PARAMS)

                rezTemp=cma.fmin(FitnessFunctionsForCma.AIC, PARAMS, 1,args=(X,Y,funcu,nbParam,terminalStr))
                #print rezTemp[1]
                if rezTemp[1]>1e50:
                    break
                if rezTemp[1]<fit :
                   # print rezTemp
                    fit=rezTemp[1]
                    rez=rezTemp


            prez=[0]
            try:
                #print rez
                prez=rez[0].tolist()
                me=FitnessFunctionsForCma.AIC
                call=me.called
    #            print '\n test '
                print str(call)
                print '\n'+str(fit)

            except:
                return (PARAMS,2e50,funca,0)
        #print "a iesit din cma"
            return (prez,rez[1],funca,call)
        elif ff==7:
            print('a intrat la aicx')

            for i in range(0,1):
              #  OptimizedBeforeConsts=[ -0.15,0.23,0.678,3,-0.026,0.15]
                PARAMS=[]
                if len(OptimizedBeforeConsts)>0:
                    print str(OptimizedBeforeConsts)
                    PARAMS=OptimizedBeforeConsts
                  #  print PARAMS
                else:
                    for j in range(0,nbParam):
                        PARAMS.append(random.uniform(0, 1))
                        #print str(PARAMS[j])


                #print str(PARAMS)

                rezTemp=cma.fmin(FitnessFunctionsForCma.AICx, PARAMS, 1,args=(X,Y,funcu,nbParam,terminalStr))
                #print rezTemp[1]
                if rezTemp[1]>1e50:
                    break
                if rezTemp[1]<fit :
                   # print rezTemp
                    fit=rezTemp[1]
                    rez=rezTemp


            prez=[0]
            try:
                #print rez
                prez=rez[0].tolist()
                me=FitnessFunctionsForCma.AICx
                call=me.called
    #            print '\n test '
                print str(call)
                #print '\n'+str(fit)

            except:
                return (PARAMS,2e50,funca,0)
        #print "a iesit din cma"
            return (prez,rez[1],funca,call)
        elif ff==6:
            print('a intrat la complexityMeasures')

            for i in range(0,1):
              #  OptimizedBeforeConsts=[ -0.15,0.23,0.678,3,-0.026,0.15]
                PARAMS=[]
                if len(OptimizedBeforeConsts)>0:
                    print str(OptimizedBeforeConsts)
                    PARAMS=OptimizedBeforeConsts
                  #  print PARAMS
                else:
                    for j in range(0,nbParam):
                        PARAMS.append(random.uniform(0, 1))
                        #print str(PARAMS[j])


                #print str(PARAMS)

                rezTemp=cma.fmin(FitnessFunctionsForCma.ComplexityMeasures, PARAMS, 1,args=(X,Y,funcu,nbParam,terminalStr))
                #cma.plot();
                #cma.show
                #print rezTemp[1]
                if rezTemp[1]>1e50:
                    break
                if rezTemp[1]<fit :
                   # print rezTemp
                    fit=rezTemp[1]
                    rez=rezTemp


            prez=[0]
            try:
                #print rez
                prez=rez[0].tolist()
                me=FitnessFunctionsForCma.ComplexityMeasures
                call=me.called
                print(call)
    #            print '\n test '
                print str(call)
                print '\n'+str(fit)

            except:
                return (PARAMS,2e50,funca,0)
        #print "a iesit din cma"
            return (prez,rez[1],funca,call)   
        elif ff==9:
            print('a intrat la bFF')

            for i in range(0,1):
              #  OptimizedBeforeConsts=[ -0.15,0.23,0.678,3,-0.026,0.15]
                PARAMS=[]
                if len(OptimizedBeforeConsts)>0:
                    print str(OptimizedBeforeConsts)
                    PARAMS=OptimizedBeforeConsts
                  #  print PARAMS
                else:
                    for j in range(0,nbParam):
                        PARAMS.append(random.uniform(0, 1))
                        #print str(PARAMS[j])


                #print str(Y)

                rezTemp=cma.fmin(FitnessFunctionsForCma.BFF, PARAMS, 1,args=(X,Y,funcu,nbParam,terminalStr))
                #cma.plot();
                #cma.show
                #print rezTemp[1]
                if rezTemp[1]>1e50:
                    break
                if rezTemp[1]<fit :
                   # print rezTemp
                    fit=rezTemp[1]
                    rez=rezTemp


            prez=[0]
            try:
                #print rez
                prez=rez[0].tolist()
                me=FitnessFunctionsForCma.BFF
                call=me.called
    #             print(call)
    # #            print '\n test '
    #             print str(call)
    #             print '\n'+str(fit)

            except:
                return (PARAMS,2e50,funca,0)
        #print "a iesit din cma"
            return (prez,rez[1],funca,call)   
        elif ff==10:
            print('a intrat la bFFnoSE')

            for i in range(0,1):
              #  OptimizedBeforeConsts=[ -0.15,0.23,0.678,3,-0.026,0.15]
                PARAMS=[]
                if len(OptimizedBeforeConsts)>0:
                    print str(OptimizedBeforeConsts)
                    PARAMS=OptimizedBeforeConsts
                  #  print PARAMS
                else:
                    for j in range(0,nbParam):
                        PARAMS.append(random.uniform(0, 1))
                        #print str(PARAMS[j])


                #print str(PARAMS)

                rezTemp=cma.fmin(FitnessFunctionsForCma.BFFNoSE, PARAMS, 1,args=(X,Y,funcu,nbParam,terminalStr))
                #cma.plot();
                #cma.show
                #print rezTemp[1]
                if rezTemp[1]>1e50:
                    break
                if rezTemp[1]<fit :
                   # print rezTemp
                    fit=rezTemp[1]
                    rez=rezTemp


            prez=[0]
            try:
                #print rez
                prez=rez[0].tolist()
                me=FitnessFunctionsForCma.BFFNoSE
                call=me.called
    #             print(call)
    # #            print '\n test '
    #             print str(call)
    #             print '\n'+str(fit)

            except:
                return (PARAMS,2e50,funca,0)
        #print "a iesit din cma"
            return (prez,rez[1],funca,call)   
    else :
        return (rez,2e50,funca,0)



