#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include "mpi.h"
////	DEFINICJE WARTOSCI OPISUJACYCH NASZA 'SZKOLE'	////
#define DAYS 5		//ILOSC DNI W TYGODNIU
#define PERIODS_PER_DAY 8		//ILOSC GODZIN LEKCYJNYCH, PRZEZ KTORE SZKOLA ZOSTAJE OTWARTA
#define CLASSROOMS 10		//ILOSC SAL
#define CLASSES 5		//ILOSC KLAS
#define TEACHERS 8		//ILOSC NAUCZYCIELI
#define MAX_SUBJECTS_PER_TEACHER 3		//MAKSYMALNA ILOSC PRZEDMIOTOW JAKICH MOZE UCZYC NAUCZYCIEL
#define SUBJECTS 6		//ILOSC PRZEDMIOTOW KTORYCH UCZY SIE KAZDA KLASA
#define PERIODS_PER_TEACHER 20		//DOCELOWA ILOSC GODZIN W TYGODNIU UCZONYCH PRZEZ NAUCZYCIELA


////	DEFINICJE WLASNOSCI ALGORYTMU	////
#define MAX_FIT 5000		//MAKSYMALNA WARTOSC DOSTOSOWANIA, MUSI BYC WIEKSZA NIZ MAKSYMALNY KOSZT OSOBNIKA
#define CONFLICT_COST 25	//KARA PRZYZNAWANA ZA KAZDY KONFLIKT KLASY/NAUCZYCIELA W ROZKLADZIE
#define PERIOD_COST 5		//KARA PRZYZNAWANA ZA KAZDA GODZINE ZA MALO/ZA DUZO UCZONA PRZEZ NAUCZYCIELA
#define WINDOW_COST 1		//KARA PRZYZNAWANA ZA KAZDE 'OKIENKO'
#define MIGRATION_COEFFICIENT 0.05  //JAKI PROCENT OSOBNIKOW MIGRUJE

////	OZNACZENIA	AKCJI MPI	////
#define START_BCAST 0
#define END_BCAST 1
#define START_ALLRED 2
#define END_ALLRED 3
#define START_RECV 4
#define END_RECV 5
#define START_SEND 6
#define END_SEND 7

#define REQUEST  1
#define REPLY    2

#define MAIN_PROCESS_ID 0

////	FUNKCJE LOSUJACE	////
#define FRAND() (rand()%10000/10000.0)
#define RANDOM_PERIOD() (rand()%(DAYS*PERIODS_PER_DAY))

////	DEFINICJA POMOCNICZA, NIE PRZYDZIELONO ZADNEGO NAUCZYCIELA/KLASY/SALI/PRZEDMIOTU	////
#define NONE -1		

//// 	DEFINICJE NAZW PLIKOW KONFIGURACYJNYCH
#define CONFIG_FNAME "config.txt"
#define SUBJECTS_FNAME "subjects.txt"
#define TEACHERS_FNAME "teachers.txt"



////	STRUKTURY UZYWANE W PROGRAMIE	////
typedef struct Tuple{		//PRZEDSTAWIA JEDNA GODZINE LEKCYJNA ODBYWAJACA SIE W JEDNEJ SALI
	int teacher_id;
	int class_id;
	int room_id;
	int subject_id;
}tuple;

struct Individual{			//OSOBNIK, POTENCJALNY ROZKLAD ZAJEC
	double fitness;
	struct Tuple genotype[DAYS*PERIODS_PER_DAY][CLASSROOMS];	
};

struct Pair{				//PARA OSOBNIKOW, UZYWANA DO KRZYZOWANIA
	struct Individual first;
	struct Individual second;
};

struct TimetableGA{			//STRUKTURA PRZEDSTAWIAJACA NASZ ALGORYTM
	int popsize;
	int ngen;
	double pcross;
	double pmut;

	struct Individual *population;
	struct Individual (*crossover)(struct Pair);
	struct Pair (*selection)(struct Individual *population, int popsize);
};

struct SubjectList{			//USTAWIENIA ILOSCI GODZIN Z DANEGO PRZEDMIOTU
	int *hours;
}subject_list;

struct TeacherList{			//USTAWIENIA JAKIE PRZEDMIOTY SA UCZONE PRZEZ KTORYCH NAUCZYCIELI
	int **subjects;	//[id nauczyciela][lista przedmiotow]
}teacher_list;

struct ClassList{			//USTAWIENIA KTORE KLASY MAJA KTORYCH NAUCZYCIELI Z KTORYCH PRZEDMIOTOW
	int **teachers;	//[id klasy][id przedmiotu] -> id nauczyciela tego przedmiotu
}class_list;
	

//// 	FUNKCJE UZYWANIE W PROGRAMIE	////

//ALGORYTM:
void initGA(struct TimetableGA *ga, int popsize, int ngen, double pcross, double pmut);		//INICJALIZUJE ALGORYTM Z PARAMETRAMI: WIELKOSC POPULACJI, ILOSC POKOLEN, PRAWDOPODOBIENSTWO KRZYZOWANIA, PRAWDOPODOBIENSTWO MUTACJI
void finalizeGA(struct TimetableGA *ga);		//ZWALNIA PAMIEC ZAALOKOWANA PRZEZ ALGORYTM
void setCrossover(struct TimetableGA *ga, struct Individual (*crossover)(struct Pair));		//USTAWIA METODE KRZYZOWANIA
void setSelection(struct TimetableGA *ga, struct Pair (*selection)(struct Individual *population, int popsize));	//USTAWIA METODE SELEKCJI
void nextGen(struct TimetableGA *ga);		//POSUWA ALGORYTM JEDNO POKOLENIE DO PRZODU, WYKONUJAC WSZYSTKIE NIEZBEDNE OPERACJE
void adjustFitness(struct TimetableGA *ga);		//ODPOWIEDNIO ZMIENIA WARTOSC DOSTOSOWANIA WSZYSTKICH OSOBNIKOW TAK ZEBY ZADNA NIE BYLA UJEMNA

//METODY STOSOWANE W PROGRAMIE:
struct Pair tournamentSelection(struct Individual *population, int popsize);		//SELEKCJA TURNIEJOWA
struct Individual singlePointCrossover(struct Pair pair);		//KRZYZOWANIE JEDNOPUNKTOWE
void Mutate(struct Individual * first, double pmut);		//MUTACJA

//USTAWIENIA ALGORYTMU:
void initClassList(struct ClassList *cl);
void freeClassList(struct ClassList *cl);
void initTeacherList(struct TeacherList *cl);
void freeTeacherList(struct TeacherList *cl);
void initSubjectList(struct SubjectList *cl);
void freeSubjectList(struct SubjectList *cl);
void loadGAParams(char *filename, int *popsize, int *ngen, double *pcross, double *pmut);

//OSOBNIK
void Individual_Init(struct Individual *ind);	//LOSOWO INICJALIZUJE OSOBNIKA
void printIndividual(char *filename, struct Individual *ind);	//WYPISUJE UZYSKANY PLAN DO PLIKU

//POZOSTALE:
void sortPopByFitness(struct Individual *population, int size);		//SORTUJE POPULACJE
int compare(const void *a, const void *b);		//KOMPARATOR UZYWANY DO SORTOWANIA
int randomTeacher(int subject_id);		//ZWRACA IDENTYFIKATOR LOSOWEGO NAUCZYCIELA UCZACEGO DANEGO PRZEDMIOTU
int firstFreeSlotID(struct Individual *ind, int period_id);		//ZWRACA INDEKS WOLNEJ SALI W DANEJ GODZINIE LEKCYJNEJ, ZWRACA NONE JESLI WSZYSTKIE SA ZAJETE


//FUNKCJA DOSTOSOWANIA:
void fitness(struct Individual *ind);		//USTAWIA W POLU fitness OSOBNIKA JEGO WARTOSC DOSTOSOWANIA
int teachersPeriods(struct Individual *ind, int teacher_id);	//ZWRACA ILOSC GODZIN UCZONYCH PRZEZ NAUCZYCIELA O DANYM ID
int classSubjectPeriods(struct Individual *ind, int class_id, int subject_id);		//ZWRACA ILOSC GODZIN DANA KLASA MA Z DANEGO PRZEDMIOTU W PLANIE
int classHasPeriod(struct Individual *ind, int class_id, int period_id); 	//ZWRACA WARTOSC LOGICZNA CZY KLASA MA ZAJECIA W TRAKCIE DANEJ GODZINY LEKCYJNEJ
int windows(struct Individual *ind, int class_id);		//OBLICZA ILOSC OKIENIEK JAKIE MA W PLANIE DANA KLASA
int conflicts(struct Individual *ind);		//OBLICZA ILOSC KONFLIKTOW W PLANIE


/////		IMPLEMENTACJE:		/////


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
		fitness(&(ga -> population[i]));
	}
	sortPopByFitness(ga -> population, ga -> popsize);
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

void adjustFitness(struct TimetableGA *ga){
	int i;
	double min = 999999;
	for(i = 0; i < ga->popsize; ++i){
		min = ga->population[i].fitness < min ? ga->population[i].fitness : min;
	}
	if(min < 0){
		for(i = 0; i < ga->popsize; ++i){
			ga->population[i].fitness -= min;
		}
	}
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

void loadGAParams(char *filename, int *popsize, int *ngen, double *pcross, double *pmut){
	printf("Laduje parametry algorytmu\n");
	FILE *in = fopen (filename,"r");
	fscanf(in, "%d %d %f %f", popsize, ngen, pcross, pmut);
	fclose(in);
}

void loadSubjectParams(char *filename, struct SubjectList *sl){
	printf("Laduje ustawienia przedmiotow\n");
	FILE *in = fopen(filename,"r");
	int i;
	for(i = 0; i < SUBJECTS-1; ++i){
		fscanf(in, "%d ", &(sl->hours[i]));
	}
	fscanf(in, "%d", &(sl->hours[SUBJECTS-1]));
	fclose(in);
}

void loadTeacherParams(char *filename, struct TeacherList *tl){
	printf("Laduje ustawienia nauczycieli\n");
	FILE *in = fopen(filename,"r");
	int i,j;
	
	for(i = 0; i < TEACHERS-1; ++i){
		for(j = 0; j < MAX_SUBJECTS_PER_TEACHER -1; ++j){
			fscanf(in, "%d ", &(tl->subjects[i][j]));
		}
		fscanf(in, "%d\n", &(tl->subjects[i][MAX_SUBJECTS_PER_TEACHER-1]));
	}
	
	for(j = 0; j < MAX_SUBJECTS_PER_TEACHER -1; ++j){
		fscanf(in, "%d ", &(tl->subjects[TEACHERS-1][j]));
	}
	fscanf(in, "%d", &(tl->subjects[TEACHERS-1][MAX_SUBJECTS_PER_TEACHER-1]));
	fclose(in);
}

struct Individual singlePointCrossover(struct Pair pair){
    struct Individual Child;
    double fit=pair.first.fitness+pair.second.fitness;
    int i,j;
    for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		
		int split_pt = rand()%CLASSROOMS;
		for(j = 0; j < CLASSROOMS; ++j){
			if(j <= split_pt)
				Child.genotype[i][j]=pair.first.genotype[i][j];
			else
				Child.genotype[i][j]=pair.second.genotype[i][j];
		}		
    }

    return Child;
};

void Mutate(struct Individual * first, double pmut){
	int i,j,k,l;
	for(i = 0; i < CLASSES; ++i){
		for(j = 0; j < SUBJECTS; ++j){
			int total = 0;
			for(k = 0; k < DAYS * PERIODS_PER_DAY; ++k){
				for(l = 0; l < CLASSROOMS; ++l){
					if(first -> genotype[k][l].subject_id == j){
						if(total < subject_list.hours[j]){
							++total;
						}
						else if(total > subject_list.hours[j]){
							first -> genotype[k][l].subject_id = NONE;
							first -> genotype[k][l].class_id = NONE;
							first -> genotype[k][l].teacher_id = NONE;
						}
						
					}
				}	
			}
			while(total < subject_list.hours[j]){
				int rand_period = RANDOM_PERIOD();
				int slot_id = firstFreeSlotID(first, rand_period);
				if(slot_id == NONE){
					continue;
				}
				else{
					first -> genotype[rand_period][slot_id].teacher_id = class_list.teachers[i][j];
					first -> genotype[rand_period][slot_id].class_id = i;
					first -> genotype[rand_period][slot_id].room_id = slot_id;
					first -> genotype[rand_period][slot_id].subject_id = j;
					++total;
				}
			}
		}
	}
	
	
	//mutacja wlasciwa
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
	for(i = 0; i < (ga -> popsize); ++i){
		fitness(&(ga -> population[i]));	
	}
	sortPopByFitness(ga -> population, ga -> popsize);
	adjustFitness(ga);


	//Rownoleglosc
	int dest, src;
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 	const int nitems=4;
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

	if(rank==size-1)
	{ 		
		dest=0;
	}
	else
	{
		dest= (rank+1);
	}
	
	
	for(i=0;i<MIGRATION_COEFFICIENT*ga->popsize;++i)
	{
		MPI_Send(&(ga->population[i].fitness),  1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
		MPI_Send(ga->population[i].genotype,  1, indv, dest, 1, MPI_COMM_WORLD);
	}

	if (rank==0)
	{
		
		src=size-1; 
	}
	else
	{
		src=rank-1; 
	}
	
	for(i=0;i<MIGRATION_COEFFICIENT*ga->popsize;++i)
	{
		MPI_Status status;
		MPI_Recv(&(ga->population[ga->popsize-i].fitness),   1, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(ga->population[ga->popsize-i].genotype,   1, indv, src, 1, MPI_COMM_WORLD, &status);
	}
	
	for(i = 0; i < ga->popsize; ++i){
			fitness(&(ga -> population[i]));	
	}
	adjustFitness(ga);
		
	MPI_Type_free(&indv);
	
}

void printIndividual(char *filename, struct Individual *ind){
	FILE *out = fopen(filename, "w");
	fprintf(out, "Oznaczenia: (id_przedmiotu, id_klasy, id_nauczyciela)\n");
	int i,j,k;
	for(i = 0; i < DAYS; ++i){
		fprintf(out, "Dzien %d:\n", i+1);
		for(j = 0; j < PERIODS_PER_DAY; ++j){
			fprintf(out, "Lekcja #%d:\t", j+1);
			for(k = 0; k < CLASSROOMS; ++k){
				if(ind->genotype[i*PERIODS_PER_DAY + j][k].class_id != NONE)
					fprintf(out, "(%d,%d,%d)\t", ind->genotype[i*PERIODS_PER_DAY + j][k].subject_id, ind->genotype[i*PERIODS_PER_DAY + j][k].class_id, ind->genotype[i*PERIODS_PER_DAY + j][k].teacher_id);
				else
					fprintf(out, "(WOLNE)\t");
			}
			fprintf(out, "\n");
		}
		fprintf(out, "\n");
	}
	fclose(out);
}


int randomTeacher(int subject_id){
	int teachers[TEACHERS];
	int total = 0;
	int i,j;
	for(i = 0; i < TEACHERS; ++i){
		teachers[i] = NONE;
	}
	for(i = 0; i < TEACHERS; ++i){
		for(j = 0; j < MAX_SUBJECTS_PER_TEACHER; ++j){
			if(teacher_list.subjects[i][j] == subject_id){
				teachers[total++] = i;
				break;
			}
		}
	}
	return teachers[rand()%total];
}

void initClassList(struct ClassList *cl){
	cl -> teachers = (int**)malloc(CLASSES * sizeof(int*));
	int i;
	for(i = 0; i < CLASSES; ++i){
		cl -> teachers[i] = (int*)malloc(SUBJECTS * sizeof(int));
		int j;
		for(j = 0; j < SUBJECTS; ++j){
			cl -> teachers[i][j] = NONE;
		}
	}
	
	for(i = 0; i < CLASSES; ++i){
		int j;
		for(j = 0; j < SUBJECTS; ++j){
			cl -> teachers[i][j] = randomTeacher(j);
		}
	}
	
}

void freeClassList(struct ClassList *cl){
	int i;
	for(i = 0; i < CLASSES; ++i){
		free(cl -> teachers[i]);
	}
	free(cl -> teachers);
}


void initSubjectList(struct SubjectList *sl){	
	sl -> hours = (int*)malloc(SUBJECTS * sizeof(int));
	int i;
	for(i = 0; i < SUBJECTS; ++i){
		sl -> hours[i] = 0;
	}
}

void freeSubjectList(struct SubjectList *sl){
	free(sl -> hours);
}

void initTeacherList(struct TeacherList *tl){
	tl -> subjects = (int**)malloc(TEACHERS * sizeof(int*));
	int i;
	for(i = 0; i < TEACHERS; ++i){
		tl -> subjects[i] = (int*)malloc(MAX_SUBJECTS_PER_TEACHER * sizeof(int));
		int j;
		for(j = 0; j < MAX_SUBJECTS_PER_TEACHER; ++j){
			tl -> subjects[i][j] = NONE;
		}
	}
	
}

void freeTeacherList(struct TeacherList *tl){
	int i;
	for(i = 0; i < TEACHERS; ++i){
		free(tl -> subjects[i]);
	}
	free(tl -> subjects);
}





int compare(const void *a, const void *b){
	struct Individual *ind_a = (struct Individual*)a;
	struct Individual *ind_b = (struct Individual*)b;
	
	return (ind_a -> fitness) > (ind_b -> fitness) ? -1 : 1;
}

int firstFreeSlotID(struct Individual *ind, int period_id){
	int i;
	for(i = 0; i < CLASSROOMS; ++i){
		if(ind -> genotype[period_id][i].room_id == NONE){
			return i;
		}
	}
	return NONE;
}

void Individual_Init(struct Individual *ind){
	int i,j,k;
	
	for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS; ++j){
			ind -> genotype[i][j].teacher_id = NONE;
			ind -> genotype[i][j].class_id = NONE;
			ind -> genotype[i][j].room_id = NONE;
			ind -> genotype[i][j].subject_id = NONE;
		}
	}
	for(i = 0; i < CLASSES; ++i){
		for(j = 0; j < SUBJECTS; ++j){
			for(k = 0; k < subject_list.hours[j];){
				int rand_period = RANDOM_PERIOD();
				int slot_id = firstFreeSlotID(ind, rand_period);
				if(slot_id == NONE){
					continue;
				}
				else{
					ind -> genotype[rand_period][slot_id].teacher_id = class_list.teachers[i][j];
					ind -> genotype[rand_period][slot_id].class_id = i;
					ind -> genotype[rand_period][slot_id].room_id = slot_id;
					ind -> genotype[rand_period][slot_id].subject_id = j;
					++k;
				}
			}
		}
	}
}

int teachersPeriods(struct Individual *ind, int teacher_id){
	int total = 0;
	int i,j;
	for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS; ++j){
			if(ind -> genotype[i][j].teacher_id == teacher_id)
				total++;
		}
	}
	return total;
}

int classSubjectPeriods(struct Individual *ind, int class_id, int subject_id){
	int total = 0;
	int i,j;
	for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS; ++j){
			if(ind -> genotype[i][j].class_id == class_id && ind -> genotype[i][j].subject_id == subject_id)
				total++;
		}
	}
	return total;
}

int classHasPeriod(struct Individual *ind, int class_id, int period_id){
	int i;
	for(i = 0; i < CLASSROOMS; ++i){
		if(ind -> genotype[period_id][i].class_id == class_id)
			return 1;
	}
	return 0;
}

int windows(struct Individual *ind, int class_id){
	int i,j,first_period, last_period;
	int total = 0;
	for(i = 0; i < DAYS; ++i){
		for(first_period = i * PERIODS_PER_DAY; !classHasPeriod(ind, class_id, first_period); ++first_period)
			if(first_period == ((i+1) * PERIODS_PER_DAY) - 1)
				break;
				
		if(first_period == ((i+1) * PERIODS_PER_DAY) - 1)
				continue;
				
		for(last_period = ((i+1) * PERIODS_PER_DAY) -1; !classHasPeriod(ind, class_id, last_period); --last_period);
		
		for(j = first_period+1; j < last_period-1; ++j){
			if(!classHasPeriod(ind,class_id,j))
				total++;
		}
	}
	return total;
}

int conflicts(struct Individual *ind){
	int total = 0;
	int i,j,k;
	for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS-1; ++j){
			int teacher_id = ind -> genotype[i][j].teacher_id;
			int class_id = ind -> genotype[i][j].class_id;
			for(k = j + 1 ; k < CLASSROOMS; ++k){
				if(ind -> genotype[i][k].teacher_id == teacher_id && teacher_id != NONE && ind -> genotype[i][k].teacher_id != NONE)
						++total;
				if(ind -> genotype[i][k].class_id == class_id && class_id != NONE && ind -> genotype[i][k].class_id != NONE)
						++total;
			}
		}
	}
	return total;
}

void fitness(struct Individual *ind){
	double cost = 0;
	int i,j,k;
	//Iloœæ konfliktów:
	cost += CONFLICT_COST * conflicts(ind);
	//printf("konf: %f\n" ,cost / CONFLICT_COST);
	//Ilosc godzin nauczycieli
	for(i = 0; i < TEACHERS; ++i){
		cost += PERIOD_COST * abs(PERIODS_PER_TEACHER - teachersPeriods(ind, i));
	}
	///////////////////////
	
	//okienka:
	for(i = 0; i < CLASSES; ++i){
		cost += WINDOW_COST * windows(ind, i);
	}
	
	//ilosc lekcji z odpowiednich przedmiotow:
	for(i = 0; i < CLASSES; ++i){
		for(j = 0; j < SUBJECTS; ++j){
			cost += PERIOD_COST * abs(subject_list.hours[j] - classSubjectPeriods(ind, i, j));
		}
	}
	
	ind -> fitness = MAX_FIT - cost;
	//printf("%f ", cost);
}


