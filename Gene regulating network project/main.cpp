//
//  main.cpp
//  Gene regulating network project
//
//  Created by Jelle Molenkamp on 17/12/2018.
//  Copyright © 2018 Jelle Molenkamp. All rights reserved.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <fstream>

const int genomelength = 100;
const int N = 100;
const int Maxgeneration = 10;
const double mutationrate = 0.1;
const double meanoffspring = 2.5;
const double k = 1.0;

struct Individual {
    std::vector<double> genome;
    std::vector<std::vector<double> > interaction;
};


void reprotwo(std::vector<Individual> &population){
    
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    std::poisson_distribution<int> poisson(meanoffspring);
    
    std::vector<Individual> newpopulation;
    
    //check for every individual how many offspring they will produce.
    for(int i = 0; i<population.size();++i){
        
        int event;
        for(;;){
            event = poisson(rng);
            if(event<4)break;
        }
        
        if(event>0){
            Individual newindividual = {population[i].genome,population[i].interaction};
            for(int j = 0; j<event;++j){
                newpopulation.push_back(newindividual);
            }
        } else if (event < 0){
            throw std::logic_error("value from poisson distribution is negative. \n");
        }
    }
    //update population
    population = newpopulation;
}








void reproduction(std::vector<std::vector<int> > &genome){
    
    //create a new matrix to store the new values.
    std::vector<std::vector<int> > newgenome;
    
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    std::poisson_distribution<int> poisson(meanoffspring);
    
    //check for every individual how many offspring they will produce.
    for(int i = 0; i<genome.size();++i){
        
        int event = poisson(rng);
        
        //put new offspring to new newgenome matrix to store them.
        if(event>0){
            for(int j = 0; j<event;++j){
                newgenome.push_back(genome[i]);
            }
        } else if (event < 0 ){
            throw std::logic_error("value from poisson distribution is negative. \n");
        }
        
    }
    //update genome matrix with newgenome matrix
    genome = newgenome;
}


void mutatwo(std::vector<Individual> &population){
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    // flip value if mutation happens.
    for(int i = 0; i<population.size();++i){
        for(int j = 0; j<genomelength; ++j){
            std::bernoulli_distribution biasedCoinFlip(mutationrate);
            
            if(biasedCoinFlip(rng)){
                if(population[i].genome[j]==1){
                    population[i].genome[j] = 0;
                } else if (population[i].genome[j]==0){
                    population[i].genome[j] = 1;
                } else{
                    throw std::logic_error("state of gene is neither 0 or 1. \n");
                }
            }
            
        }
    }
}

void mutation(std::vector<std::vector<int> > &genome){
    
    
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    // flip value if mutation happens.
    for(int i = 0; i<genome.size();++i){
        for(int j = 0; j<genomelength; ++j){
            std::bernoulli_distribution biasedCoinFlip(mutationrate);
            
            if(biasedCoinFlip(rng)){
                if(genome[i][j]==1){
                    genome[i][j] = 0;
                } else if (genome[i][j]==0){
                    genome[i][j] = 1;
                } else{
                    throw std::logic_error("state of gene is neither 0 or 1. \n");
                }
            }
            
        }
    }
    
}

void addinteraction(std::vector<std::vector<double> > &interaction){
    
    // create variable to count number of genes.
    double numbgenes = 1.0;
    //create vector to store degree of nodes.
    std::vector<double> nodeDegree(1,0.0);
    // create vector to store the attachment probabilities.
    std::vector<double> attachprob;
    
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    
    while(numbgenes<genomelength){
        
        // calculate the probabilities *********************************************
        double totaldegree = std::accumulate(nodeDegree.begin(), nodeDegree.end(), 0.0);
        
        //update prob vec with new node(s)
        attachprob.clear();
        for(int i =0; i<nodeDegree.size();++i){
            attachprob.push_back(pow((nodeDegree[i]+1)/(totaldegree+nodeDegree.size()), k));
        }
        
        // create distribution for pulling node with who to connect.
        std::discrete_distribution<int> chooseNode(attachprob.begin(),attachprob.end());
        //create distribution to determine how many connections the new gene will make.
        std::binomial_distribution<int> numbconnect(4,0.3);
        
        //calc number of connections for new node.
        int numb_of_connections;
        
        for(;;){
            numb_of_connections = numbconnect(rng);
            if(numb_of_connections<=numbgenes)
                break;
        }
        if(numb_of_connections<0)
            throw std::logic_error("number of connectios is negative \n");
        
        // ADD THE NEW INTERACTIONS TO THE MODEL ************************************
        std::vector<double> tmpinteract(4,-1);
        
        for(int i = 0; i<numb_of_connections; ++i){
            
            double attachnode;
            
            //determine the node with who the new node will connect. Only if it hasn't connect with this node before.
            for(;;){
                attachnode = chooseNode(rng);
                if(attachnode != tmpinteract[0] && attachnode != tmpinteract[1] && attachnode != tmpinteract[2] && attachnode != tmpinteract[3])
                    break;
            }
            if(attachnode>=nodeDegree.size() || attachnode<0)
                throw std::logic_error("id-value of connection node is invalid \n");
            
            // add interaction to the matrix
            std::vector<double> newinteraction = {numbgenes,attachnode};
            interaction.push_back(newinteraction);
            
            //update degree vec.
            nodeDegree[attachnode] += 1.0;
            
            //add value to tmpinteract to check if don't make connection with the same node.
            tmpinteract.push_back(attachnode);
        }
        
        //update degree vec.
        nodeDegree.push_back(numb_of_connections);
        
        ++numbgenes;
    }
}

void nondInteraction(std::vector<Individual> &population, std::vector<double> &phenotype){
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    std::vector<double> additive(population.size(),0);
    std::vector<double> epistasis(population.size(),0);
    
    
    //Calculate the additive part of phenotype the calculation
    for(int i = 0; i<population.size();++i){
        for(int j = 0; j<population[i].genome.size();++j){
            additive[i] += population[i].genome[j];
        }
    }
    std::cout<<"additive value for every individual is  =  \n";
    for(int i = 0; i<epistasis.size();++i){
        std::cout<<"individual  "<<i<<"  additive value = "<<additive[i]<<std::endl;
    }
    std::cout<<"\n";
    
    //Calculate the epistasis part of the phenotype calculation
    for(int i = 0; i<population.size();++i){
        for(int j = 0; j<population[i].interaction.size();++j){
            int geneOne = population[i].interaction[j][0];
            int geneTwo = population[i].interaction[j][1];
            //BELOW IS OPTION TO CHOOSE WHICH RULE YOU WANT TO APPLY. !!!!!!!!!!!!!!!!!!! ONLY SELECT ONE.
            epistasis[i] += (population[i].genome[geneOne]*population[i].genome[geneTwo]);              //multiplication
            //epistasis[i] += ((population[i].genome[geneOne]+population[i].genome[geneTwo])/2);          // average
            //if(population[i].genome[geneOne] == 1 || population[i].genome[geneTwo] == 1) epistasis[i]+=1.0;                                                             // at least one 1
        }
    }
    std::cout<<"epistasis value for every individual is  =  \n";
    for(int i = 0; i<epistasis.size();++i){
        std::cout<<"individual  "<<i<<"  epistasis = "<<epistasis[i]<<std::endl;
    }
    std::cout<<"\n";
    
    //Calculate the environmental noise and calculate the phenotype
    std::normal_distribution<double> environment(0.0,1.0);
    
    
    for(int i =0; i<phenotype.size();++i){
        double noise = environment(rng);
        std::cout<<"individual "<<i<<"  noise = "<<noise<<std::endl;
        phenotype[i] = additive[i] + epistasis[i] + noise;
    }
    std::cout<<"\n";
}

void direcInteraction(std::vector<Individual> &population,std::vector<double> &phenotype){
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    std::vector<double> addxepiVec(population.size(),0.0);
   

    
    // Calculate the output after the additive epistasis interactions.
    
    for(int i = 0; i<population.size();++i){
        // walk through the population and copy the genome of every individual. With the copy we are able to calc the output without changing the genome.
        std::vector<double> outputepi(population[i].genome.size(),0.0);
        for(int j = 0; j<population[i].genome.size();++j){
            outputepi = population[i].genome;
        }

        for(int j = 0; j<genomelength;++j){
            std::vector<double> inputID;
            
            //figure out by who gene j is connected and store them in a vector in order to do calculations...
            for(int q = 0; q<population[i].interaction.size();++q){
                if(population[i].interaction[q][1]== j){
                    int inputIDval =population[i].interaction[q][0];
                    inputID.push_back(population[i].genome[inputIDval]);
                }
            }
            //calculate the output of every gene after epistasis for individiual i :
            if(inputID.size()>0){
            double tmpaverage = (std::accumulate(inputID.begin(), inputID.end(), 0.0))/(inputID.size());
            //outputepi[j] = outputepi[j]*tmpaverage;                     //for recessive inhibitory epistasis
            outputepi[j] = outputepi[j]*(1.0-tmpaverage);                 //for dominant inhibitory epistasis
            }
        }
        
        addxepiVec[i] = std::accumulate(outputepi.begin(), outputepi.end(), 0.0);
    }

    //Calculate the environmental noise and calculate the phenotype
    std::normal_distribution<double> environment(0.0,1.0);
    
    
   for(int i =0; i<phenotype.size();++i){
        double noise = environment(rng);
 
        phenotype[i] = addxepiVec[i]+ noise;
   }
}

void additivemodel(std::vector<Individual> &population, std::vector<double> &phenotypeAdd){
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    std::vector<double> additive(population.size(),0.0);
    
    //Calculate the additive part of phenotype the calculation
    for(int i = 0; i<population.size();++i){
        for(int j = 0; j<population[i].genome.size();++j){
            additive[i] += population[i].genome[j];
        }
    }
    
    //Calculate the environmental noise and calculate the phenotype
    std::normal_distribution<double> environment(0.0,1.0);
    
    
    for(int i =0; i<phenotypeAdd.size();++i){
        double noise = environment(rng);
        
        phenotypeAdd[i] = additive[i] + noise;
    }
    
    std::cout<<"Phenotypeaddition vec values are  "<<std::endl;
    for(int i = 0; i<population.size();++i){
        std::cout<<phenotypeAdd[i]<<std::endl;
    }
};




int main() {
    try{
        //obtain seed from system clock
        std::chrono::high_resolution_clock::time_point tp =
        std::chrono::high_resolution_clock::now();
        unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
        
        // create and seed pseudo-random number generator
        std::mt19937_64 rng;
        std::clog<<"random seed : "<<seed<<"\n";
        rng.seed(seed);
        
        std::vector<Individual> population;
        
        // create matrix to store the genomes. Every row is the genome of a individual.
        std::vector<std::vector<double> > genome(N,std::vector<double>(genomelength,0));
        
        
        // RANDOMLY ASSIGN THE VALUE (O OR 1) TO THE GENES. ***********************************************
        for(int i = 0; i<genome.size();++i){
            for(int j = 0; j<genomelength;++j){
                //give random value 1 to genes with uniform distribution.
                const double p = 0.5; // 50% change to have value 1.
                std::bernoulli_distribution biasedCoinFlip(p);
                
                if(biasedCoinFlip(rng))
                    genome[i][j] = 1;
                
            }
        }
        
        //DETERMINE EDGE ARCHITECTURE *********************************************************************
        // create matrix to store the interactions in
        std::vector<std::vector<double> > interaction;
        
        addinteraction(interaction);
        
        // PUT THE GENOME VECTORS AND INTERACTIONS TO POPULATION STRUCT ************************************
        for(int i = 0; i<genome.size();++i){
            Individual dataindividual = {genome[i],interaction};
            population.push_back(dataindividual);
        }
        
        
        // LET TE POPULATION EVOLVE BY CREATING NEW GENERATIONS ********************************************
        //To get to a gaussian distribution, we need to let the population evolve...
        
        for(int genCount = 0; genCount<Maxgeneration; ++genCount){
            
            //create new offspring
            //reproduction(genome);
            reprotwo(population);
            
            //look for mutations.
            //mutation(genome);
            mutatwo(population);
            
            std::cout<<genCount<<" is done "<<std::endl;
        }
        std::cout<<"generation calc is done"<<std::endl;
        std::cout<<"population size = "<<population.size()<<std::endl;
        // CALCULATE THE FREQUENCIES OF EXPRESSION VALUES ***************************************************
        
        // create two vectors. One for calculating the phenotype with epistasis and one for calculating the phenotype without epistasis.
        std::vector<double> phenotype(population.size(),0.0);
        std::vector<double> phenotypeAdd(population.size(),0.0);
      
        
        
        
        //nondInteraction(population, phenotype);
        direcInteraction(population, phenotype);
        additivemodel(population, phenotypeAdd);
        
        
        // PRINT TO COMPUTER TO CHECK THE DATA **********************************************************
        for(int i = 0; i<population.size();++i){
            std::cout<<"individual  =  "<<i<<std::endl;
            std::cout<<"genome = "<<std::endl;
            for(int j = 0;j<genomelength;++j){
                std::cout<<population[i].genome[j]<<", ";
            }
            std::cout<<"\n";
            std::cout<<"\n";
        }
        std::cout<<"the interaction matrix is ="<<std::endl;;
        for(int i = 0; i<population[0].interaction.size();++i){
            for(int j = 0; j<2;++j){
                std::cout<<population[0].interaction[i][j]<<"   ";
            }
            std::cout<<"\n";
        }
        
        std::cout<<"\n";
        std::cout<<"phenotype of every individual is  =  \n";
        for(int i = 0; i<phenotype.size();++i){
            std::cout<<"individual  "<<i<<"  phenotype = "<<phenotype[i]<<std::endl;
        }
        
        
        
        
        
        // OUTPUT TO FILE **********************************************************
        //open a output file.
        std::ofstream ofs;
        ofs.open("output.csv");
        if(!ofs.is_open()){
            throw std::runtime_error("unable to open file output.csv \n");
        }
        
        //put the phenotype values to excel in order to make
        for(int i = 0; i<population.size();++i){
            ofs<<i+1<<","<<phenotypeAdd[i]<<","<<phenotype[i]<<"\n";
        }
        
        
    }
    catch(std::exception &fatalException){
        std::cerr <<"fatal error: "<<fatalException.what();
        exit(EXIT_FAILURE);
    }
    
    
    return 0;
}
