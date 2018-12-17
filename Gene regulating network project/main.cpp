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

void addinteraction(std::vector<std::vector<int> > &interaction){
    
    // create variable to count number of genes.
    int numbgenes = 1;
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
        
        // ADD THE NEW INTERACTIONS TO THE MODEL ************************************
        std::vector<int> tmpinteract(4,-1);
        
        for(int i = 0; i<numb_of_connections; ++i){
            
            int attachnode;
            
            //determine the node with who the new node will connect. Only if it hasn't connect with this node before.
            for(;;){
                attachnode = chooseNode(rng);
                if(attachnode != tmpinteract[0] && attachnode != tmpinteract[1] && attachnode != tmpinteract[2] && attachnode != tmpinteract[3])
                    break;
            }
            
            // add interaction to the matrix
            std::vector<int> newinteraction = {numbgenes,attachnode};
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
        
        //DETERMINE EDGE ARCHITECTURE *********************************************************************
        // create matrix to store the interactions in
        std::vector<std::vector<int> > interaction;
        
        addinteraction(interaction);
        
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
        
        for(int i = 0; i<interaction.size();++i){
            for(int j = 0; j<2;++j){
                std::cout<<interaction[i][j]<<", ";
            }
            std::cout<<"\n "; 
        }
        
    }
    catch(std::exception &fatalException){
        std::cerr <<"fatal error: "<<fatalException.what();
        exit(EXIT_FAILURE);
    }
    
    
    return 0;
}
