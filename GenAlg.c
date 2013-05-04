#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Individual.h"



struct Pair{
	struct Individual first;
	struct Individual second;
};

struct TimetableGA{
	int popsize;
	int ngen;
	double pcross;
	double pmut;

	struct Individual *population;
	struct Individual (*crossover)(struct Pair, double pcross);
	struct Pair (*selection)(struct Individual *population, int popsize);	
};

void sortPopByFitness(struct Individual *population, int size){
	qsort(population, size, sizeof(struct Individual), compare);
}

void initGA(struct TimetableGA *ga, int popsize, int ngen, double pcross, double pmut){
	ga -> popsize = popsize;
	ga -> ngen = ngen;
	ga -> pcross = pcross;
	ga -> pmut = pmut;
	
	ga -> population = (struct Individual*)malloc(popsize * sizeof(struct Individual));
}

void setCrossover(struct TimetableGA *ga, struct Individual (*crossover)(struct Pair, double pcross)){
	ga -> crossover = crossover;
}

void setSelection(struct TimetableGA *ga, struct Pair (*selection)(struct Individual *population, int popsize)){
	ga -> selection = selection;
}

void finalizeGA(struct TimetableGA *ga){
	free(ga -> population);
}

struct Pair tournamentSelection(struct Individual *population, int popsize){
	int k;
	struct Pair pair;
	for(k = 0; k < 2; ++k){
		int i;
		struct Individual fighters[8];
		double fit_sum = 0;
		for(i = 0; i < popsize; ++i){
			fit_sum += population[i].fitness;
		}
		for(i = 0; i < 8; ++i){
			double rd = (rand()%100000)/100000.0 * fit_sum;
			double cur_sum = population[0].fitness;
			int j = 0;
			while(rd > cur_sum && j < popsize){
				++j;
				cur_sum += population[j].fitness;
			}
			fighters[i] = population[j];
		}
		
		sortPopByFitness(fighters, 8);
		if(k == 0)
			pair.first = fighters[0];
		else
			pair.second = fighters[0];
	}
	return pair;
}

int main(){
	struct TimetableGA ga;
	initGA(&ga, 100,100,1,1);
	printf("dupa");
	finalizeGA(&ga);
}
