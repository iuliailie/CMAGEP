#include "Chromosome.h"

#include <string>
#include <iostream>
#include <vector>
using namespace std;
//bool operator< (Chromosome &cC1, Chromosome &cC2)
//{
//    return cC1.optimizedFitness <=cC2.optimizedFitness;
//}

Chromosome::Chromosome(int geneNb,char _linkFunction)
{
    numberOfGenes=geneNb;
    linkFunction=_linkFunction;
    Genes=new std::vector<Gene *>();
    Values=new double[dataSize];
    optimizedParameters=new std::vector<double> ();
}

Chromosome::Chromosome()
{
    Genes=new std::vector<Gene *> ();
    Values= new double[dataSize];
    optimizedParameters=new std::vector<double> ();
}

Chromosome *Chromosome::copyChromosome()
{
    Chromosome *crom=new Chromosome();
    crom->numberOfGenes=this->numberOfGenes;
    crom->linkFunction=this->linkFunction;
    crom->fitness=this->fitness;
    crom->OrfSum=this->OrfSum;
    crom->optimizedFitness=this->optimizedFitness;
    crom->optimizedMathExpression=this->optimizedMathExpression;
    crom->mathExpression=this->mathExpression;
    crom->simplifiedMathExpression=this->simplifiedMathExpression;
    crom->maxFitnessPossible=this->maxFitnessPossible;
    crom->optimizedParameters=this->optimizedParameters;

    for(int i=0; i<crom->numberOfGenes; i++)
    {
        Gene *copiedGene=0;
        copiedGene=this->Genes->at(i)->copyGene();
        crom->Genes->push_back(copiedGene);
    }

    return crom;
}

void Chromosome::computeChromosomeValue(double *p,map<char,FunctionStructure> *localMap,int j)

{
    // cout<<"for value="<<p[0];
    // cout<<endl;
    vector<double> vector1, vector2, vector3;//initialize working vectors
    this->OrfSum=0;

//
//    for(int i=0;i<3;i++)
//    {
////        for (int tip=0;tip<3;tip++)
////        { i+=tip;
//            cout<<*(p+i)<<" ";
//
////        }
//
//    cout<<endl;
//    }
//    cout<<endl<<"done"<<endl;

    for(int CromIterator=0; CromIterator<this->numberOfGenes; CromIterator++)
    {
        Gene *gene;
        gene=(this->Genes->at(CromIterator));
        Node *node=&(gene->NodeMap->at(0));
        //cout<<gene->geneString<<endl;
        double val=computeTreeValue(node,p,localMap);

        if(node->param!=1)
        {
            val=node->param*val;
            cout<<node->param<<"----"<<endl;
        }

        this->OrfSum=this->OrfSum+gene->ORF;

        if(vector1.size()==0)//if  we are at the first gene, the chromosome value is equal to the gene value
            vector1.insert(vector1.begin(),val);
        else//else we apply the linking function to the values from the former genes and the current one
        {
            vector2.insert(vector2.begin(),val);
            vector3.insert(vector3.begin(),0);
            transform(vector1.begin(),vector1.end(),vector2.begin(),vector3.begin(),localMap->at(this->linkFunction).pFunction);
            vector1.at(0)=vector3.at(0);
        }
    }

//    if(this->Values->size()==nr)
//        this->Values->clear();
//cout<<endl<<vector1.at(0);
    this->Values[j]=vector1.at(0);
    vector1.clear();
    vector2.clear();
    vector3.clear();
}


void Chromosome::Mutate()
{
    string functii=functionSet;
    string terminale=terminalsToUse;
    // select a random gene
    int k=rand()%this->numberOfGenes;
    // select a random position in the gene
    string tempHead=this->Genes->at(k)->_head;
    string tempTail=this->Genes->at(k)->_tail;
    // in the head we can change to a function or a terminal
    int index=rand()%(tempHead.size()+tempTail.size());

    if(index < this->Genes->at(k)->headLength)
    {
        // note that the first position in the head must be a function
        char temp;

        if(index == 0)
        {
            // mutate to a (different) function because the head can only remain a function
            temp=functii.at(rand()%functii.size());
            tempHead.replace(index,1,1,temp);//replace for length 1 for 1 time starting from index
        }
        else
        {
            // mutate to a (different) function or a terminal
            temp=(functii+terminale).at(rand()%((functii+terminale).size()));
            tempHead.replace(index,1,1,temp);
        }
    }
    else
    {
        // in the tail we can only mutate to a terminal
        char temp;
        temp=terminale.at(rand()%terminale.size());
        tempTail.replace(index-this->Genes->at(k)->headLength,1,1,temp);
    }

    Gene *old=this->Genes->at(k);
    this->Genes->at(k)=new Gene(tempHead,tempTail);
    delete old;
    // return true;.
}





void Chromosome::doInversion()
{
    int k=rand()%this->numberOfGenes;//randomly choose for which gene will the inversion be performed
    string head=this->Genes->at(k)->_head;
    string tail=this->Genes->at(k)->_tail;
    int i1=rand()%this->Genes->at(k)->headLength;
    int i2=rand()%this->Genes->at(k)->headLength+i1;
    string substr=head.substr(i1,i2-i1);
    std::reverse(substr.begin(),substr.end());
    head.replace(i1,i2-i1,substr);
    Gene *old=this->Genes->at(k);
    this->Genes->at(k)=new Gene(head,tail);
    delete old;
}




void Chromosome::doISTransposition()
{
    // select a random gene
    int k=rand()%this->numberOfGenes;
    Gene * randomGene=this->Genes->at(k);
    string head=randomGene->_head;
    string tail=randomGene->_tail;
    string localGeneString=randomGene->geneString;
    int beginPoint,endPoint,targetPoint;
    beginPoint=getARandom(1,localGeneString.length()-1);
    targetPoint=getARandom(1,randomGene->headLength-1);
    //make sure that the length for the IS does not surpass the length of the head
    endPoint=getARandom(beginPoint,beginPoint+head.length()-targetPoint-1);
    string IS=localGeneString.substr(beginPoint,endPoint);
    head.insert(targetPoint,IS);
    head.erase(head.length()-1-IS.length(),IS.length());
    Gene *old=this->Genes->at(k);
    randomGene=0;
    this->Genes->at(k)=new Gene(head,tail);
    delete old,randomGene;
}

void Chromosome::doRootTransposition()
{
    string functii=functionSet;
    // select a random gene
    int k=rand()%this->numberOfGenes;
    Gene * randomGene=this->Genes->at(k);
    string head=randomGene->_head;
    string tail=randomGene->_tail;
    string localGeneString=randomGene->geneString;
    int beginPoint,endPoint,targetPoint;
    beginPoint=getARandom(1,head.length()-1);
    int i;

    for(i=beginPoint; i<head.length(); i++)
    {
        if(isFunction(head.at(i)))
            break;
        else
        {
            if(i==head.length()-1)
                return;
        }
    }

    beginPoint=i;
    targetPoint=0;//transpose to the root
    //make sure that the length for the IS does not surpass the length of the head
    endPoint=getARandom(beginPoint,beginPoint+3);
    string RIS=localGeneString.substr(beginPoint,endPoint+1);

    if(head.length()>=RIS.length())
    {
        head.erase(0,RIS.length());
        head.insert(0,RIS);
    }

//    else{
//        head=RIS.substr(0,head.length());
//        int diff=RIS.length()-head.length();
//        tail.erase(0,diff);
//        tail.insert(0,RIS.substr(head.length(),RIS.length()-1));
//    }
//  if(head.length()+tail.length()==localGeneString.length())
//this->Genes->at(k)=new Gene(head,tail);
//else {
//
//i=10;
//}
}




void Chromosome::GeneSwap(Chromosome * Chrom2)
{
    // Choose a random gene from both chromosomes, and swap them
    int g1= rand() % this->numberOfGenes;
    int g2 = rand() % Chrom2->numberOfGenes;
    Gene *gene1=new Gene(this->Genes->at(g1)->_head,this->Genes->at(g1)->_tail);
    Gene *gene2=new Gene(Chrom2->Genes->at(g2)->_head,Chrom2->Genes->at(g2)->_tail);
    Gene *oldGene1,*oldGene2;
    oldGene1=this->Genes->at(g1);
    oldGene2=Chrom2->Genes->at(g2);
    Chrom2->Genes->at(g2)=gene1;
    this->Genes->at(g1)=gene2;
    gene1=0;
    gene2=0;
    delete gene1,gene2;
    delete oldGene1,oldGene2;
}

void Chromosome::GeneSwap(Chromosome * Chrom2, int geneNb)
{
    // Choose a random gene from both chromosomes, and swap them
    int g1= geneNb;
    int g2 =geneNb;
    Gene *gene1=new Gene(this->Genes->at(g1)->_head,this->Genes->at(g1)->_tail);
    Gene *gene2=new Gene(Chrom2->Genes->at(g2)->_head,Chrom2->Genes->at(g2)->_tail);
    delete this->Genes->at(g1),Chrom2->Genes->at(g2);
    Chrom2->Genes->at(g2)=gene1;
    this->Genes->at(g1)=gene2;
    gene1=0;
    gene2=0;//clean the pointers
    delete gene1,gene2;
}

void Chromosome::doOnePointRecombination(Chromosome * otherChrome)
{
    //cout<<"before--------------"<<endl;
//
//            cout<<"c1"<<endl;//!!pentru testare selectie/evolutie cromozomi
//        for(int k=0;k<this->numberOfGenes;k++)
//       {
//           cout<<this->Genes->at(k)->geneString<<endl;
//
//       }
//       cout<<endl;
//                cout<<"c2"<<endl;//!!pentru testare selectie/evolutie cromozomi
//        for(int k=0;k<otherChrome->numberOfGenes;k++)
//       {
//           cout<<otherChrome->Genes->at(k)->geneString<<endl;
//
//       }
//       cout<<endl;
    int g=getARandom(0,this->numberOfGenes-1);//choosing the gene
    int rP=getARandom(0,this->Genes->at(0)->geneString.length()-1);//choose crossover point

    for(int i=g+1; i<this->numberOfGenes; i++) //we swap the genes which come after the one in which the crossover is done
    {
        if(i!=g)
            this->GeneSwap(otherChrome,i);
    }

    string geneS1=this->Genes->at(g)->geneString;
    string geneS2=otherChrome->Genes->at(g)->geneString;
    string subStr1=geneS1.substr(rP);
    string subStr2=geneS2.substr(rP);
    geneS1.replace(rP,string::npos,subStr2);
    geneS2.replace(rP,string::npos,subStr1);
    int headLength=this->Genes->at(g)->headLength;
    delete this->Genes->at(g),otherChrome->Genes->at(g);
    this->Genes->at(g)=new Gene(geneS1.substr(0,headLength),geneS1.substr(headLength,string::npos));
    otherChrome->Genes->at(g)=new Gene(geneS2.substr(0,headLength),geneS2.substr(headLength,string::npos));
//cout<<"after--------------"<<g<<"------"<<rP<<"------"<<endl;
//          cout<<"c1"<<endl;//!!pentru testare selectie/evolutie cromozomi
//        for(int k=0;k<this->numberOfGenes;k++)
//       {
//           cout<<this->Genes->at(k)->geneString<<endl;
//
//       }
//       cout<<endl;
//                cout<<"c2"<<endl;//!!pentru testare selectie/evolutie cromozomi
//        for(int k=0;k<otherChrome->numberOfGenes;k++)
//       {
//           cout<<otherChrome->Genes->at(k)->geneString<<endl;
//
//       }
//       cout<<endl;
}

void Chromosome::doTwoPointRecombination(Chromosome * other)
{
    string  fullString1="";
    string fullString2="";

    for(int i=0; i<this->numberOfGenes; i++) //we generate a string which will contain all the genes' string for ease of double crossover
    {
        fullString1+=this->Genes->at(i)->geneString;
        fullString2+=other->Genes->at(i)->geneString;
    }

    int p1=getARandom(0,fullString1.length()-2);//we choose the two crossover points, the first can be at most the char before the last
    int p2=getARandom(p1+1,fullString1.length()-1);
    int hL=this->Genes->at(0)->headLength;//get headLength for future use
    string substr1=fullString1.substr(p1,p2-p1);//get the substrings which will be exchanged betweent the two chromes
    string substr2=fullString2.substr(p1,p2-p1);
    fullString1.replace(p1,p2-p1,substr2);
    fullString2.replace(p1,p2-p1,substr1);

    for(int i=0; i<this->numberOfGenes; i++)
    {
        delete this->Genes->at(i), other->Genes->at(i);
        string sH1,sH2,sT1,sT2;
        sH1=fullString1.substr((2*hL+1)*i,hL+(2*hL+1)*i-(2*hL+1)*i);
        sT1=fullString1.substr(hL+(2*hL+1)*i,2*hL+1+(2*hL+1)*i-(hL+(2*hL+1)*i));
        this->Genes->at(i)=new Gene(sH1,sT1);
        sH2=fullString2.substr((2*hL+1)*i,hL+(2*hL+1)*i-(2*hL+1)*i);
        sT2=fullString2.substr(hL+(2*hL+1)*i,2*hL+1+(2*hL+1)*i-(hL+(2*hL+1)*i));
        other->Genes->at(i)=new Gene(sH2,sT2);
    }
}

string Chromosome::getChromosomeMathExpression()
{
    string ChromExpression="";
//    cout<<"Chrom 360"<<endl;

    for(int k=0; k<this->numberOfGenes; k++)
    {
        string  geneExpression=printExpresie(&this->Genes->at(k)->NodeMap->at(0));

        if(k<this->numberOfGenes-1)
            ChromExpression+=geneExpression+this->linkFunction;
        else
            ChromExpression+=geneExpression;
    }

//    cout<<"Chrom 372"<<endl;
    return ChromExpression;
}





