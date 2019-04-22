#ifndef FORTESTS_H
#define FORTESTS_H


class ForTests
{
public:
    ForTests();
    virtual ~ForTests();
protected:
private:
};
static void forTests()
{
    /*
           sorted_vector<Chromosome *>ChromList;
        Gene *gene;
        gene=new Gene("+a*+","aaaa");
        gene->GenerateNodeMap(localMap);
        Node *node=&(gene->NodeMap[0]);
        double val=computeTreeValue(node,valueA[0],localMap);


        return 0;
        */
//just for output of initial genes
//    for (int j=0;j<ChromList.size();j++){
//        Chromosome *crom;
//        crom=ChromList.at(j);
//        cout<<j<<endl;
//        for(int i=0;i<crom->numberOfGenes;i++)
//        { cout<<crom->Genes->at(i).toString()<<endl;
//        //cout<<crom->Genes->at(0).ORF;
//        }
//        cout<<endl;
//    }
//  map<char,FunctionStructure> localMap;
//    Global *g=new Global;
//      g->mapCharactersWithFunctions();
//    localMap=g->functionDictionary;
//
//Gene *gene=new Gene("*aa*","aaaa");
//gene->getORF_Size(localMap);
//cout<<gene->ORF;
    /*
      sorted_vector<Chromosome *>SelectedChromosomes;
        SelectedChromosomes=gep->SelectChromosomes(ChromList);
        SelectedChromosomes=SelectChromosomes(SelectedChromosomes);
     vector<string> temp;

            for(unsigned int i=0; i<ChromList.size(); i++)
            {
                temp.push_back(SelectedChromosomes.at(i)->Genes->at(0)->geneString);
                SelectedChromosomes.at(i)->Mutate();
            }

            cout<<"------------------"<<endl;

            for(unsigned int i=0; i<ChromList.size(); i++)
            {
                cout<<"         "<<endl;
                cout<<temp.at(i)<<endl;
                cout<<SelectedChromosomes.at(i)->Genes->at(0)->geneString<<endl;
            }

            cout<<"------------------"<<endl;
            */
}
#endif // FORTESTS_H
