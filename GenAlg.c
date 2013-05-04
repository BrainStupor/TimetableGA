#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Individual.h"

struct Pair{
	struct Individual *first;
	struct Individual *second;
};

struct TimetableGA{
	int popsize;
	int ngen;
	double pcross;
	double pmut;

	struct Individual *population;
	struct Individual (*crossover)(struct Pair);
	struct Pair (*selection)(struct Individual *population, int popsize, double pcross);	
};


void initGA(struct TimetableGA *ga, int popsize, int ngen, double pcross, double pmut){
	ga -> popsize = popsize;
	ga -> ngen = ngen;
	ga -> pcross = pcross;
	ga -> pmut = pmut;
	
	ga -> population = (struct Individual*)malloc(popsize * sizeof(struct Individual));
}

void setCrossover(struct TimetableGA *ga, void (*crossover)(struct Individual *population, int popsize)){
	ga -> crossover = crossover;
}

void setSelection(struct TimetableGA *ga, struct Pair (*selection)(struct Individual *population, int popsize, double pcross)){
	ga -> selection = selection;
}

void finalizeGA(struct TimetableGA *ga){
	free(ga -> population);
}

int main(){
	struct TimetableGA ga;
	initGA(&ga, 100,100,1,1);
	printf("dupa");
	finalizeGA(&ga);
}
