#include <iostream>
#include <cstring>
#include <stdlib.h>     
#include <time.h>       
#include <array>
#include <vector>
#include <math.h> 
#include <functional>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <fstream>
#include <string>
#include "parameters.hpp"

using namespace std; 

const float OUT_OF_RANGE                    = numeric_limits<float>::max();
const int   elite_size                      = community_size*elite_percentual/100;
const int   new_generation_size             = community_size - elite_size;
const int   mating_pool_size                = community_size*mating_pool_percentual/100;
const int   number_of_parents               = 2;
const int   number_of_couples               = new_generation_size/number_of_parents;


class Gene
{
    public:
        Gene();
 
        Gene(float a);
        
        void setGene(Gene* newGene);

        void setGene_Float(float number);

        float getGene_Float(void);

        unsigned int * getGene_Binary(void);

        void setBit(unsigned int a,unsigned int position);
        
        void flipBit(unsigned int position);
 
        unsigned int getBit(unsigned int position);

        void printGene_Binary(void);
  
        void printGene_Float(void);
    
    private:
        void updateFloatGene(void);

        double floatGene;
        unsigned int binaryGene[integer_bits+decimal_bits];

};

class Individual
{
    public:

        Individual(Gene a[]);

        Individual(void);

        void printChromosome_Binary(void);
        
        void printChromosome_Float(void);
        
        void printCoast(void);

        array<Gene*, number_of_genes> getChromosome(void);

        void coastUpdate(float (*coastFunction)(array<float, number_of_genes> ));

        float getCoast(void);

        void setChromosome(array<Gene, number_of_genes>& newGene);

        void resetChromosome(void);
        
        unsigned int getGeneBit(unsigned int geneNumber, unsigned int bitLocation);
    

    private:
        float randomFloat(void);

        array<Gene*, number_of_genes> genes;
        float coast;
        
};

class GeneticAlgorithm{
public:

    GeneticAlgorithm(float (*cfun)(array<float, number_of_genes> )): myCoastFunction(cfun)

    
    {
        hammelThreshold                 = ((decimal_bits+integer_bits)*number_of_genes);
        showResultCtrl                  = false;
        showCommunityCtrl               = false;
        for(int i =0; i < number_of_genes; i++) solutionFloat[i] = 0;
        solutionCoast = numeric_limits<float>::max();
        
    }
    ~GeneticAlgorithm(){
    }
    float* run(void){

        srand (time(NULL));

        for(int i = 0; i < executions; i++){
            ofstream myFileSave;

            if(saveCtrl){
                myFileSave.open (myFileSaveName+to_string(i)+".csv", ofstream::trunc | ofstream::in);
                while(myFileSave.is_open() == false){
                    myFileSave.open (myFileSaveName+to_string(i)+".csv", ofstream::trunc | ofstream::in);
                }
            }
            // Set random new gene values 
            resetCommunity();

            unsigned int iterations = 0;
            
            // Update the coast of each solution of the community and sort it in descending order
            CoastFunctionEvaluation ();
            
            // Select couples from mating pool, the approach is defined by user
            Selection               ();

            while(iterations < max_iterations && PrecisionAchieved() == false){

                if(showCommunityCtrl) showCommunity(iterations);
                
                // Save some convergence measurement
                if(saveCtrl)          coastSave(ref(myFileSave));
    
                // Create offspring from selected parents
                Mating                  ();

                // Modify new solutions 
                Mutation                ();

                // Prepare the news solutions adding information from the offsprings genes
                NextGeneration          ();

                // Update the coast of each solution of the community
                CoastFunctionEvaluation ();

                // Separe the elite members from the rest
                Selection               ();
                
                iterations++;
                
            }
            
            // Show solution 
            if(showResultCtrl) showSolution    ();
            // Save final answers
            if(saveCtrl) myFileSave.close     ();

            saveSolution();
        }


        return(getSolution());
    }

    void enableCommunityDisplay(void){showCommunityCtrl = true;}
    void disableCommunityDisplay(void){showCommunityCtrl = false;}

    void enableResultDisplay(void){showResultCtrl = true;}
    void disableResultDisplay(void){showResultCtrl = false;}

    void enableCoastHistory(string name){
        myFileSaveName = name;
        saveCtrl = true;
    }
    void disableCoastHistory(void){
        saveCtrl = false;
    }
    
private:
    void resetCommunity();
    bool PrecisionAchieved();
    void CoastFunctionEvaluation();
    bool isIncest(Individual& , Individual& );
    void Selection  ();
    void parentsCrossover(Individual *, Individual *, array<array<Gene,number_of_genes>,2>& );
    void Mating();
    void Mutation();
    void NextGeneration();
    void coastSave(ofstream& );
    void genesSave(ofstream& );
    void saveSolution();
    float coastsAverage();
    void showCommunity(int );
    void showSolution();
    float* getSolution();

    array  < Individual                           , community_size       >     community_members;
    array  < array<Individual*,number_of_parents> , number_of_couples    >     couples; 
    array  < array<Gene,number_of_genes>          , new_generation_size  >     offsprings_genes;
   
    array<float,number_of_genes>                                                solutionFloat;
    float                                                                       solutionCoast;

    unsigned int    hammelThreshold;

    string          myFileSaveName;

    bool            showCommunityCtrl;
    bool            showResultCtrl;
    bool            saveCtrl;

    float (*myCoastFunction)(array<float, number_of_genes> );

};





