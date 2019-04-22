#include "Global.h"
#include "GEP.h"


using namespace std;

//here we set the data size and number of explanatory variables
const int number_of_fitness_cases=5000;
const int number_of_parameters=18;
const string allPossibleCharsForTermials="ABDGIJKMNOQRUVWXYZbdjkuvwyz";

double valueA[number_of_fitness_cases][number_of_parameters];
double Target[number_of_fitness_cases];

string writingFile;
string inputReadFile;
string targetReadFie;
string writingFileStats;
string currentResultFile;
double ValidationInput[20];
double ValidationTarget[20];





void oneGEPRun(int k, clock_t currentTime, int finali);
int editDistance(string *s, string *t);
int myMin(int a, int b, int c);
void simplifyTestFunc();
void optimizedGep();
void readInputFromFile(int  dataSize, string file, bool input);

long seedgen(void)
{
    long s, seed, pid;
    pid = getpid();
    s = time(NULL);    /* get CPU seconds since 01/01/1970 */
    seed = abs(((s*181)*((pid-83)*359))%104729);
    cout<<"seed "<<seed<<endl;
    return seed;
}


#include <string>
#include <sstream>
#include <vector>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}



int main()
{
    Py_Initialize();//initialize python module
    PyRun_SimpleString(
        "import sys\n"
        "sys.path.append('./')\n"
        "sys.path.append('./sympy')\n"      
    );

	//initate the random seed
    cout<<"a intrat in main";
    string currentDataSize="dataSizeDeInlocuit";
    stringstream ssDataSize;
    ssDataSize<<currentDataSize.substr(0,currentDataSize.length()-1);
    ssDataSize>>dataSize;
    dataSize=dataSize-5;
    //cout<<currentDataSize<<endl;
    cout<<"data size "<<dataSize<<endl;

    nParam=number_of_parameters;
    clock_t t1,t2;
    stringstream ss,ssStats;
    unsigned long int sec= time(NULL);     
    string currentDate="dataDeInlocuit";
    currentResultFile="currentNode";
	currentResultFile=currentResultFile.substr(0,currentResultFile.length()-3);
    ss<<"./Results/"<<currentDate<<"/Result"<<currentResultFile<<".csv";
    ssStats<<"./Results/"<<currentDate<<"/Stats"<<currentResultFile<<".txt";
    std::ofstream out("./Results/"+currentDate+"/Out"+currentResultFile+".txt");
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
      std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    const char *c_name = ss.str().c_str();
    string asda=ss.str();
    writingFile=c_name;
    const char *c_nameStats = ssStats.str().c_str();
    string asZda=ssStats.str();
    writingFileStats=c_nameStats;
    fstream outfile;
    outfile.open(writingFile,ios::out | ios::trunc);
    outfile<<"fitness,max fitness possible,structure,simplified structure,optimized strucure,tree size,optimizedFitness,type of fitness function,inputFile,number of generations,fitness evaluations,runtime"<<"\n";
    outfile.close();
    int k=0;   
    int simplifyTest=0;

    if(!simplifyTest)
    {
        float diff;

       k=0 ;

        while(k>=0)
        {
            t1=clock();
            oneGEPRun(9,t1,k);
            diff=(double)(clock() - t1)/CLOCKS_PER_SEC;
            k--;
            outfile.open(writingFile,fstream::in | fstream::out | fstream::app);
            outfile<<(int)diff<<"\n";           
            outfile.close();
            out.flush();
        }

    }


    Py_Finalize();

     std::cout.rdbuf(coutbuf); //reset to standard output again

    return 0;
}


void oneGEPRun(int k, clock_t currentTime, int finali)
{
	string currentFolder="folderDeInlocuit";
	string numarDeInl="numarDeInlocuit";
    cout<<numarDeInl<<endl;
    cout<<numarDeInl.at(0)<<endl;
    cout<<numarDeInl.at(0) - '0'<<endl;
  inputReadFile=currentFolder+"/XLagDataSet"+to_string(10*(numarDeInl.at(0) - '0')+finali)+".txt";
    targetReadFie=currentFolder+"/YLagDataSet"+to_string(10*(numarDeInl.at(0) - '0')+finali)+".txt";
     	/*inputReadFile=currentFolder+"/X.txt";
    targetReadFie=currentFolder+"/Y.txt"; */
    cout<<inputReadFile;
    readInputFromFile(dataSize, inputReadFile, true);  //true for input file false for target file
    readInputFromFile(dataSize,  targetReadFie, false);
    //inititate the function alphabet and variables
    //set the parameters for the chromosomes
    int populationSize=112;
    functionSet="+-*/^SLE";//SELsC +-*/^SLsCE;  +-*/^SLE
   terminals=allPossibleCharsForTermials.substr(0,number_of_parameters);
   terminalsToUse="ABGOQUVXY";//fara P H F T S L E C g
   constants="123456789";
    double *ptr=&valueA[0][0];
    actualData=&ptr;//assign the observations
    targetData=Target;
    int  headLength=5;
    int nrGene=3;
    char linkF='+';
    int maxNumberOfCycles=0;
    double maxTimeToRun=900;
    double timeToStartCMA=0;
    bool useCMA=true;
    double currentTimeRun;
    

    //generating the first population of chrmosomes randomly
    GEP * gep=new GEP();
    vector<Chromosome *> *ChromList= gep->generatePopulation(headLength,nrGene,linkF,populationSize);
    map<char,FunctionStructure> *localMap;
    Global *g=new Global;
    genesHeadLength=headLength;
    numberOfGenesInAChromosome=nrGene;
    g->mapCharactersWithFunctions();
    localMap=g->functionDictionary;
    //set the treshhold accepted for difference between the predicted and actual data
    gep->fitnessType=(FitnessType)(k);    
    gep->precisionForFitness=0.01;

//set parameteres for the mutation rates for the algorithm
    //double mutationValue=(double)2/((2*headLength+1)*nrGene);
    gep->mutationRate=0.5;//equivalent to 2 point mutation per chromosome
    gep->inversionRate=0.05;
    gep->onePRRate=0.4;
    gep->twoPRRate=0.3;
    gep->geneRRate=0.1;
    gep->IS_transpRate=0.05;
    gep->RIS_transpRate=0.1;
    gep->geneTranspRate=0.1;
    //first training cycle where we  evaluate  the expresion trees and the fitnesses in the chromosomes of the first population
    /*saving the settings of the runs*/
    string currentDate="dataDeInlocuit";
    fstream settingsFile;
    string settingsWritingFile;
    string fisierSettings="./Results/"+currentDate+"/Settings.txt";
    //string fisierSettings="./Results/Settings.txt";
    stringstream ss;
    ss<<fisierSettings;
    const char *c_n = ss.str().c_str();
    string asda=ss.str();
    settingsWritingFile=c_n;
    settingsFile.open(settingsWritingFile,fstream::in | fstream::out | fstream::trunc);
    settingsFile<<"samplesize "<<dataSize<<"\n";
    settingsFile<<"populationSize "<<populationSize<<"\n";
    settingsFile<<"functionSet "<<functionSet<<"\n";
    settingsFile<<"terminals "<<terminals<<"\n";
    settingsFile<<"constants "<<constants<<"\n";
    settingsFile<<"headLength "<<headLength<<"\n";
    settingsFile<<"nrGene "<<nrGene<<"\n";
    settingsFile<<"linkF "<<linkF<<"\n";
    settingsFile<<"maxNumberOfCycles "<<maxNumberOfCycles<<"\n";
    settingsFile<<"maxTimeToRun "<<maxTimeToRun<<"\n";
    settingsFile<<"timeToStartCMA "<<timeToStartCMA<<"\n";
    settingsFile<<"gep->fitnessType "<< gep->fitnessType<<"\n";
    settingsFile<<"gep->mutationRate "<<gep->mutationRate<<"\n";
    settingsFile<<"gep->geneTranspRate "<<gep->geneTranspRate<<"\n";
    settingsFile<<"gep->inversionRate "<<gep->inversionRate<<"\n";
    settingsFile<<"gep->onePRRate "<<gep->onePRRate<<"\n";
    settingsFile<<"gep->twoPRRate "<<gep->twoPRRate<<"\n";   
    settingsFile<<"gep->geneRRate "<< gep->geneRRate<<"\n";
    settingsFile<<"gep->IS_transpRate "<<gep->IS_transpRate<<"\n";
    settingsFile<<"gep->RIS_transpRate"<<gep->RIS_transpRate<<"\n";
    settingsFile.close();
    /*ends here*/
    gep->doTrainingCycle(ChromList,Target,valueA,localMap,false,true);
    vector<Chromosome *> *SelectedChromosomes;
    //if(currentTime>=timeToStartCMA)useCMA=true;
    SelectedChromosomes=gep->SelectChromosomesOnOptimizedFitness(ChromList,useCMA);//initial attribution of fitnesses, no optimization on this step useCMA
    delete ChromList;
    double currentMaxFitness=1e50;
    double MaxUnimproovedGenerationsAllowed=50;//we set how many generations can have the previous fitness before stop
    int NumberOfUnimproovedGenerations=0;
    //!!---------------begin GEP
    int j;

    for(j=0; j<maxNumberOfCycles; j++)
    {
//        cout<<"cycle "<<j<<endl;
        clock_t initTimeCycle=clock();

        if(SelectedChromosomes->at(0)->optimizedFitness<=SelectedChromosomes->at(0)->maxFitnessPossible)
        {
            //cout<<"best fitness found at run "<<j<<endl;
            clock_t endTimeCycle=clock();
            gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
            cout<<"--------------------"<<j<<endl;
            break;
        }

        for(unsigned int i=1; i<SelectedChromosomes->size(); i++)
        {
            if(getARandom(0.0,1.0)<gep->mutationRate)
              { SelectedChromosomes->at(i)->Mutate();//
                continue;
              }

            if(getARandom(0.0,1.0)<gep->inversionRate)
               {  SelectedChromosomes->at(i)->doInversion();
                continue;
              }

            int other=rand()%(SelectedChromosomes->size()-1)+1;

            if(getARandom(0.0,1.0)<gep->geneRRate)
              {   SelectedChromosomes->at(i)->GeneSwap(SelectedChromosomes->at(other));
                continue;
              }

            if(getARandom(0.0,1.0)<gep->IS_transpRate)
              {   SelectedChromosomes->at(i)->doISTransposition();
                continue;
              }

            if(getARandom(0.0,1.0)<gep->RIS_transpRate)
             {    SelectedChromosomes->at(i)->doRootTransposition();
                continue;
              }

            other=rand()%(SelectedChromosomes->size()-1)+1;

            if(getARandom(0.0,1.0)<gep->onePRRate)
             {    SelectedChromosomes->at(i)->doOnePointRecombination(SelectedChromosomes->at(other));
                continue;
              }

            other=rand()%(SelectedChromosomes->size()-1)+1;

            if(getARandom(0.0,1.0)<gep->onePRRate)
           {      SelectedChromosomes->at(i)->doTwoPointRecombination(SelectedChromosomes->at(other));
                continue;
              }
        }

        gep->doTrainingCycle(SelectedChromosomes,Target,valueA,localMap,false,false);
       
        //select best chromes
        vector<Chromosome *> * SelectedChromosomesOld=SelectedChromosomes;
        currentTimeRun=(double)((clock() - currentTime)/CLOCKS_PER_SEC);



        if(currentTimeRun>=timeToStartCMA)
            useCMA=true;

        SelectedChromosomes=gep->SelectChromosomesOnOptimizedFitness(SelectedChromosomesOld,useCMA);

        if(SelectedChromosomesOld!=SelectedChromosomes)
        {
            for(int l=0; l<SelectedChromosomesOld->size(); l++)
            {
                delete SelectedChromosomesOld->at(l);
            }

            delete SelectedChromosomesOld;
            //free(SelectedChromosomesOld);
        }//we delete the previous generation from memory

        currentTimeRun=(double)((clock() - currentTime)/CLOCKS_PER_SEC);

        if(maxTimeToRun<=currentTimeRun)
        {
            clock_t endTimeCycle=clock();
            gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
            break;
        }

        if(currentMaxFitness>SelectedChromosomes->at(0)->optimizedFitness)//we check how many generations we have without improvement
        {
            currentMaxFitness=SelectedChromosomes->at(0)->optimizedFitness;
            NumberOfUnimproovedGenerations=0;
        }
        else
        {
            NumberOfUnimproovedGenerations++;

            if(NumberOfUnimproovedGenerations>MaxUnimproovedGenerationsAllowed)
            {
                //cout<<"best fitness found at run "<<j<<endl;
                cout<<"--------------------"<<j<<endl;
                clock_t endTimeCycle=clock();
                gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
                break;
            }
        }

        if(SelectedChromosomes->at(0)->optimizedFitness<=SelectedChromosomes->at(0)->maxFitnessPossible)
        {
            cout<<"best fitness found at run "<<j<<endl;
            cout<<"--------------------"<<j<<endl;
            clock_t endTimeCycle=clock();
            gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
            break;
        }

        if(j==maxNumberOfCycles)
        {
            cout<<"best fitness found at run "<<j<<endl;
            cout<<"--------------------"<<j<<endl;
        }

        //cout<<j<<endl;
        clock_t endTimeCycle=clock();
        gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;

        currentTimeRun=(double)((clock() - currentTime)/CLOCKS_PER_SEC);

 //save best chromes of previous generation in hist file
        stringstream ssHist;
        ssHist<<"./Results/"+currentDate+"/History"+to_string(finali)+".csv";
        const char *c_nameHist = ssHist.str().c_str();
        fstream outfileHist ;
        outfileHist .open(c_nameHist,fstream::in | fstream::out | fstream::app);
        for(int indexH=0;indexH<10;indexH++)
        {outfileHist<<SelectedChromosomes->at(indexH)->optimizedMathExpression<<" , "<<SelectedChromosomes->at(indexH)->optimizedFitness;     
        }
        outfileHist<<"\n";

        
       
       
    }

	//for optimizing the last generation
	
		vector<Chromosome *> * SelectedChromosomesOld=SelectedChromosomes;  
        SelectedChromosomes=gep->SelectChromosomesOnOptimizedFitness(SelectedChromosomesOld,true);
        if(SelectedChromosomesOld!=SelectedChromosomes)
        {
            for(int l=0; l<SelectedChromosomesOld->size(); l++)
            {
                delete SelectedChromosomesOld->at(l);
            }

            delete SelectedChromosomesOld;
            
        }//we delete the previous generation from memory
     


    //!!-----------end GEP
//    cout<<"main 264";
    fstream outfile;
    outfile.open(writingFile,fstream::in | fstream::out | fstream::app);
    gep->doTrainingCycle(SelectedChromosomes,Target,valueA,localMap,true,false);
    outfile<<SelectedChromosomes->at(0)->fitness<<","<<SelectedChromosomes->at(0)->maxFitnessPossible<<",";
    outfile<<SelectedChromosomes->at(0)->mathExpression;
    outfile<<","<< SelectedChromosomes->at(0)->simplifiedMathExpression;
    outfile<<","<<SelectedChromosomes->at(0)->optimizedMathExpression;
    if(SelectedChromosomes->at(0)->optimizedParameters->size()>0)
        outfile<<","<<SelectedChromosomes->at(0)->optimizedParameters->size();//save tree size
    else
        outfile<<","<<SelectedChromosomes->at(0)->OrfSum;

    outfile<<","<<SelectedChromosomes->at(0)->optimizedFitness;
    outfile<<","<<gep->fitnessType;
	outfile<<","<<inputReadFile;
    outfile<<","<<(j+1);//saving number of generations in the results file
    outfile<<","<<gep->fitFuncCalls+NBCalls<<",";//saving number of fitness eval in the results file
    outfile.close();
    fstream outfileStats;
    outfileStats.open(writingFileStats,fstream::in | fstream::out | fstream::app);

    for(int indexstats=0; indexstats<=gep->statistics->size(); indexstats++)
    {
        outfileStats<<currentResultFile<<" "<<finali<<" "<<indexstats+1<<" ";

        if(indexstats<gep->statistics->size())
            for(int istats=0; istats<gep->statistics->at(indexstats).size(); istats++)
                outfileStats<<gep->statistics->at(indexstats).at(istats)<<" ";
        else
        {
            outfileStats<<SelectedChromosomes->at(0)->optimizedFitness<<" ";

            for(int istats=1; istats<gep->statistics->at(indexstats-1).size(); istats++)
                outfileStats<<gep->statistics->at(indexstats-1).at(istats)<<" ";
        }

        outfileStats<<endl;
    }

    outfileStats<<"\n"<<endl;
    outfileStats.close();

    /*//for testing purposes
        cout<<"-----------------"<<endl;
        double *pointer=SelectedChromosomes->at(0)->Values;
        for(int i=0;i<number_of_fitness_cases;i++)
        {
        cout<<endl<<*(pointer+i);
        }
        cout<<"-----------------"<<endl;
        */
    for(int i=0; i<SelectedChromosomes->size(); i++)
    {
        SelectedChromosomes-> at(i)->optimizedParameters->clear();
        // delete  SelectedChromosomes-> at(i)->optimizedParameters;
    }

    SelectedChromosomes->clear();
    gep->statistics->clear();
    delete gep, g;
    delete SelectedChromosomes;
    cout<<"GEP run done"<<endl;
}


void readInputFromFile(int  dataSize, string file, bool input)
{
    srand(seedgen());
    vector<string> STRING;
    ifstream infile;
    infile.open(file);
    vector<string> readStrings;
    string s;
    int t=0;

    while(infile.good()) // To get you all the lines.
    {
        STRING.push_back("");
        getline(infile,STRING.at(t)); // Saves the line in STRING.
        t++;
    }

    infile.close();
    //int k= getARandom((int)0,187);
    int k=0;

    for(int i=k; i<k+dataSize; i++)
    {
        istringstream f(STRING.at(i));

        if(input)
        {
            while(getline(f, s, '\r'))
            {
                readStrings=split(s.c_str(),'\t');
            }

            //cout<<readStrings.size()<<endl;

            for(int l=0; l<number_of_parameters; l++)
            {
                //cout<<readStrings.at(l)<<endl;
                valueA[i-k][l]=atof(readStrings.at(l).c_str());
            }
        }
        else
        {
            while(getline(f, s, '\r'))
            {
                readStrings.push_back(s);
            }
            Target[(i)]=atof(readStrings.at(i).c_str());

            /*if(atof(readStrings.at(i).c_str())>0)*/
            /*Target[(i)]=log (atof(readStrings.at(i).c_str()));
            else 
            Target[(i)]=1e-05;*/
        }
    }
}




int editDistance(string *s, string *t)
{
    int m,n;
    m=s->length();
    n=t->length();
    int d[m+1][n+1];

    for(int i=0; i<=m; i++)
        d[i][0]=i;

    for(int j=0; j<=n; j++)
        d[0][j]=j;

    for(int j=1; j<=n; j++)
        for(int i=1; i<=m; i++)
        {
            if(s->at(i-1)==t->at(j-1))
                d[i][j]=d[i-1][j-1];
            else
                d[i][j]=myMin(d[i-1][j] + 1,  // a deletion
                              d[i][j-1] + 1,  // an insertion
                              d[i-1][j-1] + 1 // a substitution
                             );
        }

//    cout<<endl;
//
//    for(int i=0; i<=m; i++)
//    {
//        for(int j=0; j<=n; j++)
//            cout<< d[i][j]<<" ";
//
//        cout<<endl;
//    }

    cout<<endl<<"distanta este "<<d[m][n];
    return d[m][n];
}


int myMin(int a, int b, int c)
{
    int temp=1e5;

    if(temp>a)
        temp=a;

    if(temp>b)
        temp=b;

    if(temp>c)
        temp=c;

//    cout<<endl<<a<<" "<<b<<" "<<c<<" "<<endl;
//    cout<<temp<<endl;
    return temp;
}






