#ifndef GEP_H
#define GEP_H
#include "Global.h"
#include "Chromosome.h"
#include "FitnessFunction.h"


#include <string>
#include <iostream>
#include <vector>
#include <utility>
using namespace std;
//extern const int nbParam;//number of parameters must be set here


class GEP
{

public:
    GEP();
    double precisionForFitness;
    double mutationRate;
    double inversionRate;
    double onePRRate; //one point recombination rate
    double twoPRRate; //two point rec. rate
    double geneRRate;//gene recomb rate
    double IS_transpRate;//IS transposition rate
    double RIS_transpRate;//IS transposition rate
    double geneTranspRate;
    vector<std::array<double, 7> > * statistics;
    double fitFuncCalls;
    FitnessType fitnessType;


    vector<Chromosome *>  *generatePopulation(int HeadLength,int geneNumber, char linkingFunction ,  int populationSize);

    int randSelector(double totalFitness,vector<Chromosome *> *ChromList);//random wheighted selection
    int randSelectorTournament(vector<Chromosome *> *ChromList);
    vector<Chromosome *> *SelectChromosomes(vector<Chromosome *> *ChromList);//selecting the new generation of chromosomes before mutations
    vector<Chromosome *> *SelectChromosomesOnOptimizedFitness(vector<Chromosome *> *ChromList, bool useOptimization);
    int randSelectorOpt(double totalFitness,vector<Chromosome *> *ChromList);
    int randSelectorTournamentOpt(vector<Chromosome *> *ChromList);

    template <size_t nbParam> int doTrainingCycle(vector<Chromosome *> *ChromList,double Target[],
            double valueA[][nbParam],
            map<char,FunctionStructure> *localMap,
            bool last, bool first)
    //has to be implemented here because it is a function template and these aren't compiled until called
    {
        double  *act;
        act=Target;
        double *p;

        for(unsigned int k=0; k<ChromList->size(); k++)
        {

            for(unsigned j=0; j<ChromList->at(k)->numberOfGenes; j++)
            {
                Gene * gene=(ChromList->at(k)->Genes->at(j));
                gene->GenerateNodeMap(localMap);//creating the corresponding expression tree for every gene in the current chromosome

            }
        }

        if(first)
        {

            for(int j=0; j<dataSize; j++)
            {
                p= valueA[j];
                ChromList->at(0)->computeChromosomeValue(p,localMap, j);
            }

            ChromList->at(0)->mathExpression=ChromList->at(0)->getChromosomeMathExpression();

            FitnessFunction * ff=new FitnessFunction();
            ff->type=this->fitnessType;
            ff->predicted=(ChromList->at(0)->Values);
            ff->actual=act;
            ff->sampleSize=dataSize;
            ff->maxProgramSize=ChromList->at(0)->numberOfGenes*(ChromList->at(0)->Genes->at(0)->headLength+ChromList->at(0)->Genes->at(0)->tailLength);
            ff->minProgramSize=ChromList->at(0)->numberOfGenes;
            ff->programSize=ChromList->at(0)->OrfSum;
            ff->precision=this->precisionForFitness;
            ChromList->at(0)->fitness=ff->ComputeFitness();//computing the fitness for each chromsome
            this->fitFuncCalls++;
            ChromList->at(0)->maxFitnessPossible=1e-06;
            delete ff;
            ChromList->at(0)->optimizedFitness=ChromList->at(0)->fitness;

        }

        if(!last)// !last if it's not the generation which contains the best fitted chromosome, for which we already know the fitness
        {
            for(unsigned int k=1; k<ChromList->size(); k++)
            {

                for(int j=0; j<dataSize; j++)
                {
                    p= valueA[j];
                   // cout<<valueA[j]<<endl;
                    ChromList->at(k)->computeChromosomeValue(p,localMap, j);
                }

                ChromList->at(k)->mathExpression=ChromList->at(k)->getChromosomeMathExpression();

                FitnessFunction * ff=new FitnessFunction();
                ff->type=this->fitnessType;
                ff->predicted=(ChromList->at(k)->Values);

//                for (int index=0;index<dataSize;index++)
//                {
//
//                    cout<<*(ChromList->at(k)->Values+index)<<endl;
//
//                }
                ff->actual=act;
                ff->sampleSize=dataSize;
                ff->maxProgramSize=ChromList->at(k)->numberOfGenes*(ChromList->at(k)->Genes->at(0)->headLength+ChromList->at(k)->Genes->at(0)->tailLength);
                ff->minProgramSize=ChromList->at(k)->numberOfGenes;
                ff->programSize=ChromList->at(k)->OrfSum;

                ff->precision=this->precisionForFitness;
                ChromList->at(k)->fitness=ff->ComputeFitness();//computing the fitness for each chromsome
                this->fitFuncCalls++;
                ChromList->at(k)->maxFitnessPossible=-1e06;
                delete ff;
                ChromList->at(k)->optimizedFitness=ChromList->at(k)->fitness;

            }

            //  cout<<"------------------------------"<<endl;
        }

        if(last)
        {

        }

        return 0;
    }

private:
};

#endif // GEP_H
