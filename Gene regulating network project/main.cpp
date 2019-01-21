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

// GENERAL PARAMETERS
const int genomelength = 100;
const int N =1000;
const int Maxgeneration = 120;
const double meanoffspring = 2.0;
const double k =    0.0;

// CALC FITNESS PARAMETERS
const double Popt = 50.0;
const double omega = 10.0;
const double metacosts = 0.002;

// PHENOTYPE CALC PARAMETER
const int maxphenotype = 100;

//TYPES OF EPISTASIS
const int direction = 1;
//1 - non-directional
//2 - directional
const int typeInteraction = 1;
//NON-DIRECTIONAL                           DIRECTIONAL
//1 - Recessive Epistasis                   1 - recessive inhibitory epistasis
//2 - Dominant Epistasis                    2 - Dominant inhibitory epistasis
//3 - Duplicate Recessive Epistasis


//RATES MUTATIONS
const double mutationrate = 0.1;
const double rpoint = 0.15;
const double rdelete = 0.3;
const double rrecruit = 0.15;
const double rdupli = 0.15;
const double rrecruitepi = 0.3;
const double rdeleteepi = 0.15;


struct Individual {
    std::vector<double> genome;
    std::vector<std::vector<int> > interaction;
    std::vector<double> nodedegree;
    long double fitness;
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
    
    
    std::vector<Individual> newpopulation;
    
    //check for every individual how many offspring they will produce.
    for(int i = 0; i<population.size();++i){
        long double meanpoisson = meanoffspring*population[i].fitness;
        std::poisson_distribution<int> poisson(meanpoisson);
        int event;
        for(;;){
            event = poisson(rng);
            if(event<3)break;
        }
        
        if(event>0){
            std::vector<double> genomeoff = population[i].genome;
            std::vector<std::vector<int> > interactoff = population[i].interaction;
            std::vector<double> nodeoff = population[i].nodedegree;
            long double fitoff = population[i].fitness;
            
            Individual newindividual = {genomeoff,interactoff,nodeoff,fitoff};
            
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



inline void mutapoint(std::vector<Individual> &population, int i, int j){
    if(population[i].genome[j]==1){
        population[i].genome[j] = 0;
    } else if (population[i].genome[j]==0){
        population[i].genome[j] = 1;
    } else{
        std::cout<<"allele = "<<population[i].genome[j]<<std::endl; 
        throw std::logic_error("state of gene is neither 0 or 1. \n");
    }
    
    
}

void mutadelete(std::vector<Individual> &population, int i, int j){
    //delete gene from genome vector and nodedegree vector.
    population[i].genome.erase(population[i].genome.begin()+j);
    population[i].nodedegree.erase(population[i].nodedegree.begin()+j);
    
    //erase interaction if it is connected with the lost gene and lower by one if is higher than lost gene --> resync interaction matrix.

    
    for(int q = 0; q<population[i].interaction.size();++q){
        if(population[i].interaction[q][0]== j || population[i].interaction[q][1]== j){
            population[i].interaction.erase(population[i].interaction.begin()+q);
        }
    }
    for(int q = 0; q<population[i].interaction.size();++q){
        if(population[i].interaction[q][0]>j) {
            --population[i].interaction[q][0];
        }
        if (population[i].interaction[q][1]>j) {
            --population[i].interaction[q][1];
        }
    }
 
    
}

inline void mutarecruit(std::vector<Individual> &population, int i){
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    
    std::bernoulli_distribution chooseallele(0.5);
    int newallele = chooseallele(rng);
    
    population[i].genome.push_back(newallele);
    population[i].nodedegree.push_back(0.0);
    
}


inline void mutadupli(std::vector<Individual> &population, int i, int j){
    population[i].genome.push_back(population[i].genome[j]);
    population[i].nodedegree.push_back(population[i].nodedegree[j]);
    
}



void mutarecruitepi(std::vector<Individual> &population, int i, int j){

    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    // calculate the probabilities *********************************************
    double totaldegree = std::accumulate(population[i].nodedegree.begin(), population[i].nodedegree.end(), 0.0);


    std::vector<double> attachprob;
    
    //update prob vec with new node(s)
    for(int q =0; q<population[i].nodedegree.size();++q){
        
        attachprob.push_back(pow((population[i].nodedegree[q]+1)/(totaldegree+population[i].nodedegree.size()), k));
    }

    
    // create distribution for pulling node with who to connect.
    std::discrete_distribution<int> chooseNode(attachprob.begin(),attachprob.end());
    
  
    
    //check if interaction is already present.
    std::vector<int> checkexist;
    for(int q =0; q<population[i].interaction.size();++q){
        if(population[i].interaction[q][0] == j){
            checkexist.push_back(population[i].interaction[q][1]);
        } else if (population[i].interaction[q][1] == j){
            checkexist.push_back(population[i].interaction[q][0]);
        }
    }

    
    
    // choose node but not the same as itself.
    int attachnode;
    for(;;){
        bool exist = false;
        attachnode = chooseNode(rng);
        //check if interaction is already present.
        for(int q = 0; q<checkexist.size();++q){
            if(checkexist[q]==attachnode) exist = true;
        }
        if(attachnode != j && exist == false)break;
    }
    if(attachnode == j)
        throw std::logic_error("attached node is same as connecting gene \n");

    
    
    //double dj = static_cast<double>(j);
    //double dattachnode = static_cast<double>(attachnode);
    
    std::vector<int> newinteraction = {j,attachnode};
    population[i].interaction.push_back(newinteraction);
    
    
    // update node degreee vec
    
    population[i].nodedegree[j] += 1;
    population[i].nodedegree[attachnode] += 1;
    
    
}



void mutadeleteepi(std::vector<Individual> &population, int i, int j){
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    std::vector<double> epiweights;
    std::vector<int> epistore;
    
    for(int q = 0; q<population[i].interaction.size();++q){
        if(population[i].interaction[q][0] == j || population[i].interaction[q][1] == j){
            epistore.push_back(q);
            epiweights.push_back(1.0);
        }
    }
    
    
    if(epistore.size()!=0){
        std::discrete_distribution<int> deleteepi(epiweights.begin(),epiweights.end());
        int loseID = deleteepi(rng)
        ;
        population[i].nodedegree[population[i].interaction[epistore[loseID]][0]] -= 1;
        population[i].nodedegree[population[i].interaction[epistore[loseID]][1]] -= 1;
        population[i].interaction.erase(population[i].interaction.begin()+(epistore[loseID]));
    }

}



void nondInteraction(std::vector<Individual> &population, std::vector<double> &phenotype,std::vector<int> &frequency){
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
   // std::cout<<"additive value for every individual is  =  \n";
   // for(int i = 0; i<epistasis.size();++i){
   //     std::cout<<"individual  "<<i<<"  additive value = "<<additive[i]<<std::endl;
   // }
   // std::cout<<"\n";
    
    //Calculate the epistasis part of the phenotype calculation
    for(int i = 0; i<population.size();++i){
        for(int j = 0; j<population[i].interaction.size();++j){
            int geneOne = population[i].interaction[j][0];
            int geneTwo = population[i].interaction[j][1];
            switch (typeInteraction) {
                case 1:
                    epistasis[i] += ((population[i].genome[geneOne]+population[i].genome[geneTwo])/2);
                    break;
                case 2:
                    if(population[i].genome[geneOne] == 1 || population[i].genome[geneTwo] == 1) epistasis[i]+=1.0;
                    break;
                case 3:
                    epistasis[i] += (population[i].genome[geneOne]*population[i].genome[geneTwo]);
                    break;
                default:
                    throw std::logic_error("integer for type of interaction is incorrect \n");
                    break;
            }
        }
    }
   // std::cout<<"epistasis value for every individual is  =  \n";
   // for(int i = 0; i<epistasis.size();++i){
   //     std::cout<<"individual  "<<i<<"  epistasis = "<<epistasis[i]<<std::endl;
   // }
   // std::cout<<"\n";
    
    //Calculate the environmental noise and calculate the phenotype
    std::normal_distribution<double> environment(0.0,1.0);
    
    
   for(int i =0; i<phenotype.size();++i){
        double noise = environment(rng);
        //std::cout<<"individual "<<i<<"  noise = "<<noise<<std::endl;
        phenotype[i] = additive[i] + epistasis[i] + noise;
    }
    //std::cout<<"\n";
    
    for(int i = 0;i<maxphenotype;++i){
        for(int j = 0; j<phenotype.size();++j){
            if(phenotype[j]>=i && phenotype[j]<i+1){
                frequency[i] += 1;
            }
        }
    }
    
}

void direcInteraction(std::vector<Individual> &population,std::vector<double> &phenotype,std::vector<int> &frequency){
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
                switch (typeInteraction) {
                    case 1:
                        outputepi[j] = outputepi[j]*tmpaverage;
                        break;
                    case 2:
                        outputepi[j] = outputepi[j]*(1.0-tmpaverage);
                        break;
                    default:
                        throw std::logic_error("integer for type of interaction is incorrect \n");
                        break;
                }
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
    
    for(int i = 0;i<maxphenotype;++i){
        for(int j = 0; j<phenotype.size();++j){
            if(phenotype[j]>=i && phenotype[j]<i+1){
                frequency[i] += 1;
            }
        }
    }
    
}

void additivemodel(std::vector<Individual> &population, std::vector<double> &phenotypeAdd,std::vector<int> &frequencyAdd){
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
    
    for(int i = 0;i<maxphenotype;++i){
        for(int j = 0; j<phenotypeAdd.size();++j){
            if(phenotypeAdd[j]>=i && phenotypeAdd[j]<i+1){
                frequencyAdd[i] += 1;
            }
        }
    }
    
   // std::cout<<"Phenotypeaddition vec values are  "<<std::endl;
   // for(int i = 0; i<population.size();++i){
  //      std::cout<<phenotypeAdd[i]<<std::endl;
  //  }
    
}

void mutation(std::vector<Individual> &population){
    //obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
    std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    
    // create and seed pseudo-random number generator
    std::mt19937_64 rng;
    // std::clog<<"random seed : "<<seed<<"\n";
    rng.seed(seed);
    
    std::vector<double> vecRates = {rpoint,rdelete,rrecruit,rdupli,rrecruitepi,rdeleteepi};
    std::discrete_distribution<int> determinemutation(vecRates.begin(),vecRates.end());
    std::bernoulli_distribution mutationprob(mutationrate);
    
    for(int i = 0; i<population.size();++i){
        for(int j =0; j<population[i].genome.size();++j){
            if(mutationprob(rng)){
                int mutationtype = determinemutation(rng);
                
                
                switch (mutationtype) {
                    case 0: //point mutation
                        mutapoint(population,i,j);
                        break;
                    case 1: //gene deletion
                        mutadelete(population,i,j);
                        break;
                    case 2: // gene recruitment
                        mutarecruit(population, i);
                        break;
                    case 3: // gene duplication
                        mutadupli(population, i, j);
                        break;
                    case 4: // interaction recruitment
                        mutarecruitepi(population, i, j);
                        break;
                    case 5: // interaction deletion
                        mutadeleteepi(population, i, j);
                        break;
                    default:
                        throw std::logic_error("mutation type unknown.  \n");
                        break;
                }
            }
            
        }
    }

}


inline void calcfitness(std::vector<Individual> &population, const std::vector<double> &phenotype){
    for(int i = 0; i<population.size();++i){
        double L = std::count(population[i].genome.begin(), population[i].genome.end(), 0.0);
        double fitness = exp(-(pow((phenotype[i]-Popt), 2.0)/(2*pow(omega, 2.0)))) - L*metacosts;
        if(fitness<0.0){
            std::cout<<"fitness below 0!"<<std::endl;
            population[i].fitness = 0.0;
        } else if (fitness>1.0){
            std::cout<<"fitness above 1! "<<std::endl;
            population[i].fitness = 1.0;
        } else {
            population[i].fitness = fitness;
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
        std::vector<std::vector<int> > interaction(0);
        std::vector<double> nodedegree(genomelength,0.0);
        double W = 1.0;
        
        
        
        // PUT THE GENOME VECTORS AND INTERACTIONS TO POPULATION STRUCT ************************************
        for(int i = 0; i<genome.size();++i){
            Individual dataindividual = {genome[i],interaction,nodedegree,W};
            population.push_back(dataindividual);
        }
        
        
        
        // OUTPUT TO FILE OPEN**********************************************************
        //open a output file.
        std::ofstream ofs;
        ofs.open("output.csv");
        if(!ofs.is_open()){
            throw std::runtime_error("unable to open file output.csv \n");
        }
        ofs<<"generation,,population size,,average genes,,average active genes,,average interactions \n";
        
        
        // LET TE POPULATION EVOLVE BY CREATING NEW GENERATIONS ********************************************
        //To get to a gaussian distribution, we need to let the population evolve...
        
        
        
        for(int genCount = 0; genCount<Maxgeneration; ++genCount){
      
            //create new offspring
            reprotwo(population);
            
            //look for mutations.
            mutation(population);
            
            // create four vectors. two for calculating the phenotype with epistasis and two for calculating the phenotype without epistasis.
            std::vector<double> phenotype(population.size(),0.0);
            //std::vector<double> phenotypeAdd(population.size(),0.0);
            std::vector<int>    frequency(maxphenotype,0);
            //std::vector<int>    frequencyAdd(maxphenotype,0);
            
            switch (direction) {
                case 1:
                    nondInteraction(population, phenotype,frequency);
                    break;
                case 2:
                    direcInteraction(population, phenotype,frequency);
                    break;
                default:
                    throw std::logic_error("integer for type of direction is incorrect \n");
                    break;
            }
            
            //additivemodel(population, phenotypeAdd,frequencyAdd);
            std::vector<double> genomeactive(population.size(),0.0);
            std::vector<double> genometotal(population.size(),0.0);
            std::vector<double> interacttotal(population.size(),0.0);
            
            
            
            for(int i = 0; i <population.size();++i){
                genometotal[i] = population[i].genome.size();
                genomeactive[i] = std::count(population[i].genome.begin(), population[i].genome.end(), 1);
                interacttotal[i] = population[i].interaction.size();
            }
            
            for(int i = 0; i<population.size();++i){
            }
            
            double averagetotal = std::accumulate(genometotal.begin(), genometotal.end(), 0.0)/genometotal.size();
           
            double averageactive = std::accumulate(genomeactive.begin(), genomeactive.end(), 0.0)/genomeactive.size();
            
            double  interactionaverage = std::accumulate(interacttotal.begin(), interacttotal.end(), 0.0)/interacttotal.size();
            
            calcfitness(population, phenotype);
       
            //for(int i = 0; i<population.size();++i){
            //    std::cout<<"individual  =  "<<i<<std::endl;
            //    std::cout<<"genome = "<<std::endl;
            //    for(int j = 0;j<population[i].genome.size();++j){
            //        std::cout<<population[i].genome[j]<<", ";
            //    }
            //    std::cout<<"\n";
            //    std::cout<<"interaction matrix = "<<std::endl;
            //    for(int j = 0; j<population[i].interaction.size();++j){
            //        if(population[i].interaction[j].size() != 2){
            //            std::cout<<"empty vec in gen control "<<std::endl;
            //        } else {
            //            std::cout<<population[i].interaction[j][0]<<"   "<<population[i].interaction[j][1]<<std::endl;
            //        }
            //   }
                
            // }
           // std::cout<<"\n";
            //std::cout<<"\n";
            
            ofs<<genCount<<",,"<<population.size()<<",,"<<averagetotal<<",,"<<averageactive<<",,"<<interactionaverage<<"\n";
            
          
            if(genCount== (Maxgeneration-1)){
                ofs<<"\n";
                ofs<<"\n";
                ofs<<"\n";
                for(int i = 0; i<frequency.size();++i){
                    ofs<<i+1<<",,"<<frequency[i]<<"\n";
                }
            }
            
          
            // for(int i = 0;i<population.size();++i){
            //    std::cout<<"fitness of "<<i<<" =  "<<population[i].fitness<<std::endl;
            //    std::cout<<"interact size of "<<i<<" =  "<<population[i].interaction.size()<<std::endl;
            //    std::cout<<"genome size of "<<i<<" =  "<<population[i].genome.size()<<std::endl;
            //    std::cout<<"nodedegree size of "<<i<<" =  "<<population[i].nodedegree.size()<<std::endl;
            //}
            
            
            std::cout<<genCount<<" is done "<<std::endl;
            std::cout<<"average genes = "<<averagetotal<<std::endl;
            std::cout<<"average active genes = "<<averageactive<<std::endl;
            std::cout<<"average interactions = "<<interactionaverage<<std::endl;
            std::cout<<"population size = "<<population.size()<<std::endl;
            std::cout<<"\n";

        }
        std::cout<<"generation calc is done"<<std::endl;
        
        
        
        //put the phenotype values to excel in order to make
        //for(int i = 0; i<population.size();++i){
        //    ofs<<i+1<<","<<phenotypeAdd[i]<<","<<phenotype[i]<<"\n";
        //}
        
        // PRINT TO COMPUTER TO CHECK THE DATA **********************************************************
       // for(int i = 0; i<population.size();++i){
       //     std::cout<<"individual  =  "<<i<<std::endl;
       //     std::cout<<"genome = "<<std::endl;
       //     for(int j = 0;j<genomelength;++j){
       //         std::cout<<population[i].genome[j]<<", ";
       //     }
       //     std::cout<<"\n";
       //     std::cout<<"\n";
       // }
       // std::cout<<"the interaction matrix is ="<<std::endl;;
       // for(int i = 0; i<population[0].interaction.size();++i){
       //     for(int j = 0; j<2;++j){
       //         std::cout<<population[0].interaction[i][j]<<"   ";
       //     }
       //     std::cout<<"\n";
        //}
        
       // std::cout<<"\n";
       // std::cout<<"phenotype of every individual is  =  \n";
       // for(int i = 0; i<phenotype.size();++i){
       //     std::cout<<"individual  "<<i<<"  phenotype = "<<phenotype[i]<<std::endl;
       // }
        
    }
    catch(std::exception &fatalException){
        std::cerr <<"fatal error: "<<fatalException.what();
        exit(EXIT_FAILURE);
    }
    
    
    return 0;
}
