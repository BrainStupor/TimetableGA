#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "GenAlg.h"
#include "mpi.h"
#include <stddef.h>


int main(int argc, char *argv[]){

	MPI_Init(&argc,&argv);
 	int numprocs, myid, server, totalin, totalout, workerid;    	

	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    	MPI_Comm_rank(MPI_COMM_WORLD,&myid);	

	srand(time(0) + myid);
	struct TimetableGA ga;	
	
	initSubjectList(&subject_list); 
	loadSubjectParams(SUBJECTS_FNAME, &subject_list);
	
	initTeacherList(&teacher_list);
	loadTeacherParams(TEACHERS_FNAME, &teacher_list);
	initClassList(&class_list);
	
	int popsize;
	int ngen;
	double pcross;
	double pmut;
	loadGAParams(CONFIG_FNAME, &popsize, &ngen, &pcross, &pmut);
	
	if(myid == 1)
		pmut = 0.5;		//specjalna populacja losującas
	
	initGA(&ga, popsize,ngen,pcross,pmut);
	
	setSelection(&ga, &tournamentSelection);
	setCrossover(&ga, &singlePointCrossover);
	
	int gen;
	double prevfit=0,nextfit=0;
	for(gen = 0; gen < ga.ngen; ++gen){
		nextGen(&ga);
		
		nextfit = ga.population[0].fitness;
		if(nextfit == MAX_FIT)
			break;	
		
		prevfit = nextfit; 
		if(myid == MAIN_PROCESS_ID)
			printf("gen = %d\tbest = %f\tconf = %d\n", gen, ga.population[0].fitness, conflicts(&(ga.population[0])));
		
	}
	
	/*const int nitems=4;
    	int          blocklengths[4] = {1,1,1,1};
	
	MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Datatype indv;
	MPI_Aint     offsets[4];
	

	offsets[0] = offsetof(tuple, teacher_id);
	offsets[1] = offsetof(tuple, class_id);
	offsets[2] = offsetof(tuple, room_id);
	offsets[3] = offsetof(tuple, subject_id);


	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &indv);
    MPI_Type_commit(&indv);
	
	
	struct Individual *final_results;
	int i;
	if(myid != MAIN_PROCESS_ID){
		MPI_Send(&(ga.population[0].fitness),  1, MPI_DOUBLE, MAIN_PROCESS_ID, 1, MPI_COMM_WORLD);
		MPI_Send(ga.population[0].genotype,  1, indv, MAIN_PROCESS_ID, 1, MPI_COMM_WORLD);
	}	
	if(myid == MAIN_PROCESS_ID){
		final_results = (struct Individual*)malloc(numprocs * sizeof(struct Individual));
		final_results[0] = ga.population[0];
		for(i = 1; i < numprocs; ++i){
			MPI_Status status;
			MPI_Recv(&(final_results[i].fitness),   1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(final_results[i].genotype,   1, indv, i, 1, MPI_COMM_WORLD, &status);
		}
	}
	MPI_Type_free(&indv);
	
	sortPopByFitness(final_results, numprocs);	
	
	
	printIndividual("najlepszy.txt", &(final_results[0]));
	
	
	free(final_results);*/
	if(myid == MAIN_PROCESS_ID)
		printIndividual("najlepszy.txt", &(ga.population[0]));
	finalizeGA(&ga);
	freeClassList(&class_list);
	freeSubjectList(&subject_list);
	freeTeacherList(&teacher_list);


 	MPI_Finalize();
}

