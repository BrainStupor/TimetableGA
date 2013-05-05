#include <time.h>
#include <stdlib.h>
#define DAYS 5
#define PERIODS_PER_DAY 8
#define CLASSROOMS 10
#define CLASSES 5
#define TEACHERS 8
#define PERIODS_PER_CLASS 30
#define RANDOM_PERIOD() (rand()%(DAYS*PERIODS_PER_DAY))
#define RANDOM_TEACHER() (rand()%TEACHERS)
#define NONE -1
#define MAX_SUBJECTS_PER_TEACHER 3
#define SUBJECTS 6
#define PERIODS_PER_TEACHER 20

#define CONFLICT_COST 3
#define PERIOD_COST 2
#define WINDOW_COST 1



struct SubjectList{
	int *hours;
}subject_list;

struct TeacherList{
	int **subjects;	//[id nauczyciela][lista przedmiotow]
}teacher_list;

struct ClassList{	
	int **teachers;	//[id klasy][id przedmiotu] -> id nauczyciela tego przedmiotu
}class_list;


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

freeTeacherList(struct TeacherList *tl){
	int i;
	for(i = 0; i < TEACHERS; ++i){
		free(tl -> subjects[i]);
	}
	free(tl -> subjects);
}

struct Tuple{
	int teacher_id;
	int class_id;
	int room_id;
	int subject_id;
};

struct Individual{
	double fitness;
	struct Tuple genotype[DAYS*PERIODS_PER_DAY][CLASSROOMS];	
};

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

void fitness(struct Individual *ind){
	double cost = 0;
	int i,j,k;
	//Iloœæ konfliktów:
	for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS-1; ++j){
			int teacher_id = ind -> genotype[i][j].teacher_id;
			int class_id = ind -> genotype[i][j].class_id;
			for(k = j + 1 ; k < CLASSROOMS; ++k){
				if(ind -> genotype[i][k].teacher_id == teacher_id)
						cost += CONFLICT_COST;
				if(ind -> genotype[i][k].class_id == class_id)
						cost += CONFLICT_COST;
			}
		}
	}
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
	
	ind -> fitness = 5000 - cost;
}


