#include "Gene.h"
#include "Global.h"
#include <string>
#include <iostream>
#include <vector>

using namespace std;
Gene::Gene(string head, string tail)
{
    _head=head;
    _tail=tail;
    headLength=_head.length();
    tailLength=_tail.length();
    geneString=_head+_tail;
    this->NodeMap=new  map<int,Node>();
}


Gene::Gene()
{
}


Gene  *Gene::copyGene()//we just need a new address because the rest of properties will be computed in the new cycle
{
    Gene *copiedGene=new Gene();
    copiedGene->_head=this->_head;
    copiedGene->_tail=this->_tail;
    copiedGene->headLength=copiedGene->_head.length();
    copiedGene->tailLength=copiedGene->_tail.length();
    copiedGene->geneString=copiedGene->_head+copiedGene->_tail;
    copiedGene->NodeMap=new  map<int,Node>();
    return copiedGene;
}




void Gene::getORF_Size(map<char,FunctionStructure> *localMap)
{
    string geneString=this->geneString;
    ORF=0;
    int orf=1;//we start with an orf size of one, assuming we only have a terminal to begin with

    for(int i=0; i<geneString.size(); i++)
    {
        if(isFunction(geneString.at(i))) //if the current character is a function we increase the orf size with the arity of the function
            //deducting 1 for the function itself whcih can be considered as a parameter for anther function
            orf+=localMap->at(geneString.at(i)).Arity-1;
        else
            orf--;//if it's a terminal we deduct the number of needed parameters with one

        if(orf==0)//when we don't need anymore characters to fill the need of parameters for all the function we return the current position in the gene adding one because the count starts at 0
        {
            ORF =i+1;
            break;
        }
    }

    if(ORF<=0)
        ORF=-1;//in case of a problem we return -1
}
string Gene::toString()
{
    return _head+_tail;
}

void Gene::GenerateNodeMap(map<char,FunctionStructure> *localMap)//we create the binary tree here for the current gene
{
    // this->NodeMap=new  map<int,Node>();
    this->getORF_Size(localMap);
    int marimeOrf=this->ORF;
    string geneString=this->toString();
    int *Mark=new int[marimeOrf];

    for(int i=0; i<marimeOrf; i++)
        Mark[i]=0;

    //  Node * nodejshd=&this->NodeMap[0];
    for(int i=marimeOrf-1; i>=0; i--)
    {
        if(this->NodeMap->size()<marimeOrf)
        {
            int j=i-1;

            if(!(isFunction(geneString[i]))&&Mark[i]==0&&marimeOrf!=1)
            {
                while(j>0)
                {
                    if(isFunction(geneString[j])&&(Mark[j]==0))
                        break;

                    j--;
                }

                Node *node1, *node2, *node;

                if(localMap->at(geneString[j]).Arity==1)
                {
                    node1=new Node(geneString[i]);
                    node2=new Node(0);
                    this->NodeMap->insert(std::pair<int,Node>(i,*node1));
                    Mark[i]++;
                }
                else
                {
                    node1=new Node(geneString[i]);
                    this->NodeMap->insert(std::pair<int,Node>(i,*node1));
                    Mark[i]++;
                    map<int,Node>::iterator it;

                    if((it=this->NodeMap->find((i-1)))==this->NodeMap->end())
                    {
                        node2=new Node(geneString[i-1]);
                        this->NodeMap->insert(std::pair<int,Node>(i-1,*node2));
                    }
                    else
                        node2=&this->NodeMap->at(i-1);

                    Mark[i-1]++;
                }

                node=new Node(geneString[j],node2,node1);

                if(localMap->at(geneString[j]).Arity==1)
                {
                    node->type=UnOperator;
                }

                Mark[j]++;
                this->NodeMap->insert(std::pair<int,Node>(j,*node));
                delete node,node1,node2;
                // i--;
            }
            else
                if(isFunction(geneString[i]) && Mark[i]==1)
                {
                    while(j>0)
                    {
                        if(isFunction(geneString[j])&&(Mark[j]==0))
                            break;

                        j--;
                    }

                    Node *node1, *node2, *node;

                    if(localMap->at(geneString[j]).Arity==1)
                    {
                        node1=&this->NodeMap->at(i);
                        node2=new Node(0);
                        Mark[i]++;
                    }
                    else
                        if(marimeOrf>1)
                        {
                            node1=&this->NodeMap->at(i);
                            Mark[i]++;
                            map<int,Node>::iterator it;

                            if((it=this->NodeMap->find((i-1)))==this->NodeMap->end())
                            {
                                node2=new Node(geneString[i-1]);
                                this->NodeMap->insert(std::pair<int,Node>(i-1,*node2));
                            }
                            else
                                node2=&this->NodeMap->at(i-1);

                            Mark[i-1]++;
                        }

                    node=new Node(geneString[j],node2,node1);

                    if(localMap->at(geneString[j]).Arity==1)
                    {
                        node->type=UnOperator;
                    }

                    this->NodeMap->insert(std::pair<int,Node>(j,*node));
                    Mark[j]++;
                    delete node,node1,node2;
                    // i--;
                }
                else
                    if(marimeOrf==1)
                    {
                        Node *node1;
                        node1=new Node(geneString[0]);
                        this->NodeMap->insert(std::pair<int,Node>(0,*node1));
                        Mark[0]++;
                        delete node1;
                    }
        }
    }

    delete [] Mark;
}



