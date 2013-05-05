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



struct SubjectList{
	int *hours;
};

struct TeacherList{
	int **subjects;	//[id nauczyciela][lista przedmiotow]
};

struct ClassList{	
	int **teachers;	//[id klasy][id przedmiotu] -> id nauczyciela tego przedmiotu
};

struct SubjectList subject_list;
struct TeacherList teacher_list;
struct ClassList class_list;

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

void fitness(struct Individual *ind){
	double score = 0;
	int i,j,k;
	for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS-1; ++j){
			int teacher_id = ind -> genotype[i][j].teacher_id;
			int class_id = ind -> genotype[i][j].class_id;
			for(k = j + 1 ; k < CLASSROOMS; ++k){
				if(ind -> genotype[i][k].teacher_id == teacher_id)
						score += 10;
				if(ind -> genotype[i][k].class_id == class_id)
						score += 10;
			}
		}
	}
	ind -> fitness = score;
}


