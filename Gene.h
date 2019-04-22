#ifndef GENE_H
#define GENE_H

#include "Global.h"
#include <string>
#include <iostream>
#include <vector>
#include<map>

using namespace std;
class Gene
{
public:
    Gene();

    Gene(string head, string tail);
    ~Gene()
    {
        this->headLength=0;
        this->tailLength=0;
        this->ORF=0;
        this->_head.clear();
        this-> _tail.clear();
        this-> geneString.clear();
        this->value=0;
        NodeMap->clear();
        delete this->NodeMap;
    }


    int headLength;
    int tailLength;
    int ORF;
    string _head;
    string _tail;
    string geneString;
    void getORF_Size(map<char,FunctionStructure> *localMap);
    void GenerateNodeMap(map<char,FunctionStructure> *localMap);
    void simplifyNodeMap();
    double value;
    map<int,Node> *NodeMap;
    string toString();
    Node *simplifyNodeMap(Node *workNode);
    Gene  *copyGene();

//    Gene& operator=(const Gene& other) {
//        head_ = other.head_;
//        tail_ = other.tail_;
//        return *this;
//        }
private:

};

#endif // GENE_H
