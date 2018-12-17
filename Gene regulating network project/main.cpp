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
        
        // create matrix to store the genomes. Every row is the genome of a individual.
        std::vector<std::vector<int> > genome(N,std::vector<int>(genomelength,0));
        
    
        
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
        
        // LET TE POPULATION EVOLVE BY CREATING NEW GENERATIONS ********************************************
        //To get to a gaussian distribution, we need to let the population evolve...
        
        for(int genCount = 0; genCount<Maxgeneration; ++genCount){
            
            //create new offspring
            reproduction(genome);
            
            //look for mutations.
            mutation(genome);
        }
        
        
        // CALCULATE THE FREQUENCIES OF EXPRESSION VALUES ***************************************************
        
        std::vector<int> sumgenome(genome.size(),0);
        
        //count the total value of the genome of every individual and store in the sumgenome vector.
        for(int i = 0; i<genome.size();++i){
            for(int j = 0; j<genomelength;++j){
                sumgenome[i] += genome[i][j];
            }
        }
        
        
        // OUTPUT TO FILE **********************************************************
        //open a output file.
        std::ofstream ofs;
        ofs.open("output.csv");
        if(!ofs.is_open()){
            throw std::runtime_error("unable to open file output.csv \n");
        }
        
        //calculate the frequency of every trait expression and store it in the output file
        for(int i = 0; i<genomelength;++i){
            ofs<<i+1<<", "<<std::count(sumgenome.begin(),sumgenome.end(),i+1)<<"\n";
        }
        
    }
    catch(std::exception &fatalException){
        std::cerr <<"fatal error: "<<fatalException.what();
        exit(EXIT_FAILURE);
    }
    
    
    return 0;
}
