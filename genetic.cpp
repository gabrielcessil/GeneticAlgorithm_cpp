#include "genetic.hpp"
Gene::Gene(){
    setGene_Float(0);
}

Gene::Gene(float a){
    setGene_Float(a);
}

void Gene::setGene(Gene* newGene){
    for(int i = 0; i < integer_bits+decimal_bits; i++){
        binaryGene[i] = newGene->getBit(i);
    }
    updateFloatGene();
}

void Gene::setGene_Float(float number){
    floatGene = number;
   
    float aux_fracPart, aux_intPart;
    aux_fracPart = modf (number , &aux_intPart);   
    
    
    memset(binaryGene, 0, sizeof(binaryGene));

    int i = integer_bits -1;
    while (aux_intPart > 0 && i >= 0)
    { 
        binaryGene[i] = (int)aux_intPart % 2; 
        aux_intPart = (int)aux_intPart / 2; 
        i--; 
    }

    i = integer_bits;
    while(aux_fracPart > 0 && i < integer_bits + decimal_bits)
    {                 
        aux_fracPart = modf ((aux_fracPart*2) , &aux_intPart);
        binaryGene[i] = (int)(aux_intPart);
        i++;
    }
}

float Gene::getGene_Float(void){
    return(floatGene);
}

unsigned int * Gene::getGene_Binary(void){
    return(binaryGene);
}

void Gene::setBit(unsigned int value,unsigned int position){
    binaryGene[position] = value;
    updateFloatGene();

}

void Gene::flipBit(unsigned int position){
    binaryGene[position] = (binaryGene[position]+1)%2;
    updateFloatGene();
}

unsigned int Gene::getBit(unsigned int position){
    return(binaryGene[position]);
}

void Gene::printGene_Binary(void){
    for (int i = 0; i < integer_bits + decimal_bits; i++){
        if(i == integer_bits) cout << ".";
        cout << binaryGene[i]; 
    }
} 

void Gene::printGene_Float(void){
    cout<<setprecision(4) <<floatGene;
} 

void Gene::updateFloatGene(void)
{
    float sum = 0;
    float base = pow(2, integer_bits-1 );
    for(int i = 0; i < integer_bits + decimal_bits; i++)
    {
        sum += binaryGene[i]*base;
        base = base / 2;
        
    }
    floatGene = sum;
}

Individual::Individual(Gene a[]){
    coast = 0;
    for(int i = 0; i < number_of_genes; i++){
        genes[i] = new Gene(a[i]);
    }   
}

Individual::Individual(void){
    coast = 0;
    for(int i=0; i<number_of_genes; i++){
        genes[i] = new Gene(randomFloat());
    }
}

void Individual::printChromosome_Binary(void){
    for(int i=0; i < number_of_genes; i++){
        cout << " | ";
        genes[i]->printGene_Binary();
    }
    cout << " | ";
}

void Individual::printChromosome_Float(void){
    for(int i=0; i < number_of_genes; i++){
        cout << " | ";
        genes[i]->printGene_Float();
    }
    cout << " | ";
}

void Individual::printCoast(void){
    cout << " Coast: " <<setprecision(4) << coast;
}

array<Gene*, number_of_genes> Individual::getChromosome(void){
    // Create and return an array of integer
    return(genes);
}

void Individual::coastUpdate(float (*coastFunction)(array<float, number_of_genes> ))
{
    array<float,number_of_genes> floatGenes;
    for(int i=0; i < number_of_genes; i++) 
        floatGenes[i] = genes[i]->getGene_Float();
    coast = coastFunction(floatGenes);
}

float Individual::getCoast(void){
    return(coast);
}

void Individual::setChromosome(array<Gene, number_of_genes>& newGene)
{
    for(int i = 0; i < number_of_genes; i++){
        genes[i]->setGene(&newGene[i]);
    }
}
void Individual::resetChromosome(void){
    for(int i = 0; i < number_of_genes; i++){
        genes[i]->setGene_Float(randomFloat());
    }
}
unsigned int Individual::getGeneBit(unsigned int geneNumber, unsigned int bitLocation){
    return(genes[geneNumber]->getBit(bitLocation));
}

float Individual::randomFloat(void){
    return(
        min_gene_value + static_cast <float> (rand()) 
            / ( static_cast <float> (RAND_MAX/(max_gene_value-min_gene_value)))
    );
}

void GeneticAlgorithm::resetCommunity(){

    for(int i=0; i<community_size; i++){
        community_members[i].resetChromosome();
    }
}

bool GeneticAlgorithm::PrecisionAchieved(){
    return(community_members[0].getCoast() < desired_coast);
}

void GeneticAlgorithm::CoastFunctionEvaluation( )
{
    //Change coast from each individual
    for(int i =0; i < community_size; i++)
    {
        community_members[i].coastUpdate(myCoastFunction);
    }
    //Sort in descending order of coast value
    sort(community_members.begin(), community_members.end(), [](Individual a, Individual b) {
        return a.getCoast() < b.getCoast();
    });
}

bool GeneticAlgorithm::isIncest(Individual& male, Individual& female){

    // Calculating the Hammel distance: count the different bits 
    unsigned int distance = 0;
    for(int i = 0; i < number_of_genes; i++)
        for(int j = 0; j < decimal_bits + integer_bits; j++)
            if( male.getGeneBit(i,j) != female.getGeneBit(i,j)) distance++;
    // Compare the Hammel distance between them with the threshold
    return(distance < hammelThreshold);

}

void GeneticAlgorithm::Selection()
{
    int male =0, female =1;
    
    int coupleIndex = 0;
    // Randomic
    if(AlgorithmChoice == 0){
        while(coupleIndex < number_of_couples){
            couples[coupleIndex][male]   = &community_members[rand()%mating_pool_size];
            couples[coupleIndex][female] = &community_members[rand()%mating_pool_size];
            coupleIndex++;
        }
    }
    // Top-down 
    else if(AlgorithmChoice == 1){
        while(coupleIndex < number_of_couples){
            couples[coupleIndex][male]   = &community_members[coupleIndex];
            couples[coupleIndex][female] = &community_members[(coupleIndex+1)%mating_pool_size];
            coupleIndex++;
        }
    }
    // Tournament
    else if(AlgorithmChoice == 2){
        while(coupleIndex < number_of_couples){
            // Tournment to define the father and the mother 
            float auxCoastFather = numeric_limits<float>::max();
            float auxCoastMother = numeric_limits<float>::max();
            for(int i = 0; i < 3; i++){
                Individual* auxIndiv = &(community_members[rand() % mating_pool_size]);
                if(auxIndiv->getCoast() < auxCoastFather){
                    couples[coupleIndex][male] = auxIndiv;
                    auxCoastFather = auxIndiv->getCoast();
                }
                auxIndiv = &(community_members[rand() % mating_pool_size]);
                if(auxIndiv->getCoast() < auxCoastMother){
                    couples[coupleIndex][female] = auxIndiv;
                    auxCoastMother = auxIndiv->getCoast();
                }
            }

            coupleIndex++;
        }
    }

    // Lottery
    else if(AlgorithmChoice == 3){

        double sum = 0;
        float rangeTicket[mating_pool_size];
        float maxTicket = 0;

        // Coast summation 
        for(int i=0; i < mating_pool_size; i++){
            sum += community_members[i].getCoast();
        }

        // Define ticket ranges for each mating pool community
        for(int i=0; i < mating_pool_size; i++){
            rangeTicket[i] = 1 - community_members[i].getCoast()/sum;
            maxTicket += rangeTicket[i];
        }
        // Define couples
        while(coupleIndex < number_of_couples){

            float maleTicket    = static_cast <float> (rand())   / ( static_cast <float> (RAND_MAX/(maxTicket)));
            float femaleTicket  = static_cast <float> (rand())   / ( static_cast <float> (RAND_MAX/(maxTicket)));            

            float auxSum = 0;
            // Find who owns the male ticket
            for(int i=0; i < mating_pool_size; i++){
                auxSum  +=rangeTicket[i];
                // If owns the ticket 
                if(auxSum >= maleTicket){
                    couples[coupleIndex][male] = &community_members[i];
                    break;
                }
            }
            auxSum = 0;
            // Find who owns the female ticket
            for(int i=0; i < mating_pool_size; i++){
                auxSum  +=rangeTicket[i];
                // If owns the ticket 
                if(auxSum >= femaleTicket){
                    couples[coupleIndex][female] = &community_members[i];
                    break;
                }
            }

            // Prepare to define next couple
            coupleIndex++;
        }
    }
    // Top-down with incest prevention
    else if(AlgorithmChoice == 4){
        while(coupleIndex < number_of_couples){
            couples[coupleIndex][male] = &community_members[coupleIndex];
            bool found  = false;
            // Search for other individue to make a couple: not circular approach
            for(int i = 0; i < mating_pool_size; i++)
                if(  !isIncest(community_members[coupleIndex], ref(community_members[i])) ){
                    couples[coupleIndex][female] = &community_members[i];
                    found = true;
                }
            
            if(found == false){
                coupleIndex = 0;
                hammelThreshold--;
            }
            else coupleIndex++;
        } 

    }
    
}
    

void GeneticAlgorithm::parentsCrossover(Individual *ma, Individual *fe, array<array<Gene,number_of_genes>,2>& childsChromosome){
    // A two-point crossover aproach:
    
    int count = 0;

     // Define who's the parents for which new choromossome
    int crossoverPoint_Start = rand() % ( (integer_bits + decimal_bits)* number_of_genes);
    int crossoverPoint_End   = crossoverPoint_Start + rand() % ( (integer_bits + decimal_bits)* number_of_genes - crossoverPoint_Start );
            
    for(int i = 0; i < number_of_genes; i++){
        for(int j = 0; j < integer_bits + decimal_bits; j++){
           
            
            if(count > crossoverPoint_Start && count < crossoverPoint_End){
                childsChromosome[i][0].setBit(ma->getGeneBit(i,j) ,j);
                childsChromosome[i][1].setBit(fe->getGeneBit(i,j) ,j);
            }
            else{
                childsChromosome[i][0].setBit(fe->getGeneBit(i,j) ,j);
                childsChromosome[i][1].setBit(ma->getGeneBit(i,j) ,j);
            }
            crossoverPoint_Start++;
        }
    }
}

void GeneticAlgorithm::Mating()
{
    int male = 0, female = 1;
    
    // Random aproach for mating couples
    for(int i =0; i < number_of_couples; i++){
        array<array<Gene,number_of_genes>,2 > child_offsprings_genes = {offsprings_genes[i*2],offsprings_genes[i*2+1]};
        parentsCrossover(couples[i][male], couples[i][female], child_offsprings_genes);
    }
}

void GeneticAlgorithm::Mutation(){
    int nChanges        = (mutation_percentual/100) * (integer_bits+decimal_bits) * number_of_genes * new_generation_size;

    int  whichOffspring = 0;
    int  whichGene      = 0;
    int  whichBit       = 0;
    for(int i =0; i < nChanges; i++){
        whichOffspring = rand() % new_generation_size;               
        whichGene      = rand() % number_of_genes;
        whichBit       = rand() % (integer_bits + decimal_bits); 

        offsprings_genes[whichOffspring][whichGene].flipBit(whichBit);

    }
}

void GeneticAlgorithm::NextGeneration(){
    for(int i = 0; i < new_generation_size; i++){
        community_members[i+elite_size].setChromosome(offsprings_genes[i]);
    }
}

void GeneticAlgorithm::coastSave(ofstream& myfileCoast){
    myfileCoast << coastsAverage() <<"\n";
}
void GeneticAlgorithm::genesSave(ofstream & myfile){
    array<Gene*,number_of_genes> genes = community_members[0].getChromosome();
    for(int j =0; j < number_of_genes; j ++) myfile << genes[j]->getGene_Float() << ",";
    myfile << community_members[0].getCoast() << "\n";
}

float GeneticAlgorithm::coastsAverage(){
    float av = 0;
    for(int i =0; i< elite_size; i++) av += community_members[i].getCoast();
    return(av/elite_size);
}
void GeneticAlgorithm::showCommunity(int iteration){
    cout << "\nCommunity [" << iteration << "] :" << endl;
    for(int i = 0; i < community_size; i++){
        if(i == elite_size) 
        cout << "-----------------------------------------------------------------------------" << endl;
        if(i == mating_pool_size)
        cout << "*****************************************************************************" << endl;

        community_members[i].printChromosome_Binary();
        community_members[i].printChromosome_Float();
        community_members[i].printCoast();
        cout << endl;
    }
}
void GeneticAlgorithm::showSolution(){
        // Display solution for this execution
        cout << "\n\nSolution: " << endl;
        community_members[0].printChromosome_Binary();
        community_members[0].printChromosome_Float();
        community_members[0].printCoast();
        cout << endl;
}
float* GeneticAlgorithm::getSolution(void){
    float * ret = new float[number_of_genes]; 
    array< Gene*, number_of_genes> gene = community_members[0].getChromosome();
    for(int i =0; i < number_of_genes; i++){
        ret[i] = gene[i]->getGene_Float();
    }
    return(ret);
}

void GeneticAlgorithm::saveSolution(){
    if(solutionCoast > community_members[0].getCoast()){
        solutionCoast = community_members[0].getCoast();

        array<Gene*, number_of_genes> gene = community_members[0].getChromosome();
        for(int i =0; i< number_of_genes; i++) solutionFloat[i] = gene[i]->getGene_Float();
    }
}
