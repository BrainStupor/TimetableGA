#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Individual.h"

#define FRAND() ((double)rand()/(double)(RAND_MAX+1))

 
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
	struct Individual (*crossover)(struct Pair);
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
	int i;
	for(i = 0; i < popsize; ++i){
		Individual_Init(&(ga -> population[i]));
	}
}

void setCrossover(struct TimetableGA *ga, struct Individual (*crossover)(struct Pair)){
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

struct Individual CrossOver(struct Pair pair)
{
    struct Individual Child;
    int fit=pair.first.fitness+pair.second.fitness;
    int i,j;
    for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS; ++j){
            double rd=(rand()%fit);
            if (rd<pair.first.fitness)
            {
                Child.genotype[i][j]=pair.first.genotype[i][j];
            }
            else
            {
                Child.genotype[i][j]=pair.second.genotype[i][j];
            }

		}
    }

    return Child;
};

void Mutate(struct Individual * first, double pmut)
{
	int i,j;
	for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS; ++j){
			if(FRAND() <= pmut){
				int a,b,c,d;
				a=i;
				b=rand()%(DAYS*PERIODS_PER_DAY);
				c=j;
				d=rand()%CLASSROOMS;			
				struct Tuple tmp;
				tmp=first->genotype[b][d];
				first->genotype[b][d]=first->genotype[a][c];
				first->genotype[a][c]=tmp;
			}
		}
	}
	for(i = 0; i < CLASSES; ++i){
		for(j = 0; j < SUBJECTS; ++j){
			if(FRAND() <= pmut){
				int random_teacher = randomTeacher(j);
				class_list.teachers[i][j] = random_teacher;
				int k,l;
				for(k = 0; k < DAYS*PERIODS_PER_DAY; ++k){
					for(l = 0; l < CLASSROOMS; ++l){
						if(first -> genotype[k][l].class_id == i && first -> genotype[k][l].subject_id == j)
							first ->genotype[k][l].teacher_id = random_teacher;
					}
				}
			}
		}
	}
};

void nextGen(struct TimetableGA *ga){
	int i;
	for(i = 0; i < (ga -> popsize); ++i){
		fitness(&(ga -> population[i]));	
	}
	struct Individual *newpop = (struct Individual*)malloc(ga->popsize * sizeof(struct Individual));
	for(i = 0; i < (ga -> popsize); ++i){
		struct Pair selected = ga->selection(ga->population, ga->popsize);
		if(FRAND() <= (ga->pcross)){
			newpop[i] = ga->crossover(selected);
		}
		else{
			newpop[i] = (selected.first.fitness > selected.second.fitness) ? selected.first : selected.second;
		}
		Mutate(&newpop[i], ga->pmut);
	}
	for(i = 0; i < (ga -> popsize); ++i){
		ga -> population[i] = newpop[i];
	}
	
}

int main(){
	
	
	
	
	
	struct TimetableGA ga;
	
	
	initSubjectList(&subject_list);
	
	subject_list.hours[0] = 2;
	subject_list.hours[1] = 5;
	subject_list.hours[2] = 7;
	subject_list.hours[3] = 8;	
	subject_list.hours[4] = 6;
	subject_list.hours[5] = 8;
	
	initTeacherList(&teacher_list);
	teacher_list.subjects[0][0] = 3;
	teacher_list.subjects[1][0] = 5;
	teacher_list.subjects[2][0] = 0;
	teacher_list.subjects[2][1] = 2;
	teacher_list.subjects[3][0] = 2;
	teacher_list.subjects[3][1] = 4;
	teacher_list.subjects[4][0] = 0;
	teacher_list.subjects[4][1] = 5;
	teacher_list.subjects[5][0] = 1;
	teacher_list.subjects[5][1] = 3;
	teacher_list.subjects[6][0] = 2;
	teacher_list.subjects[6][1] = 4;
	teacher_list.subjects[7][0] = 0;
	teacher_list.subjects[7][1] = 1;
	teacher_list.subjects[7][2] = 4;
	
	initClassList(&class_list);
	
	initGA(&ga, 100,100,1,1);
	
	finalizeGA(&ga);
	freeClassList(&class_list);
	freeSubjectList(&subject_list);
	freeTeacherList(&teacher_list);
}
