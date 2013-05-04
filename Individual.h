#include <time.h>
#include <stdlib.h>
#define DAYS 5
#define PERIODS_PER_DAY 8
#define CLASSROOMS 10
#define CLASSES 5
#define TEACHERS 5
#define PERIODS_PER_CLASS 30
#define RANDOM_PERIOD() (rand()%(DAYS*PERIODS_PER_DAY))
#define RANDOM_TEACHER() (rand()%TEACHERS)
#define NONE -1
#define MAX_SUBJECTS_PER_TEACHER 3
#define SUBJECTS 5



struct SubjectList{
	int *hours;
};

struct TeacherList{
	int **subjects;	//[id nauczyciela][lista przedmiotow]
};



void initSubjectList(struct SubjectList *sl){
	int * hours = sl -> hours;
	hours = (int*)malloc(SUBJECTS * sizeof(int));
	int i;
	for(i = 0; i < SUBJECTS; ++i){
		hours[i] = 0;
	}
}

void freeSubjectList(struct SubjectList *sl){
	free(sl -> hours);
}

void initTeacherList(struct TeacherList *tl, int number_of_teachers){
	int **subjects = tl -> subjects;
	subjects = (int**)malloc(number_of_teachers * sizeof(int*));
	int i;
	for(i = 0; i < number_of_teachers; ++i){
		subjects[i] = (int*)malloc(MAX_SUBJECTS_PER_TEACHER * sizeof(int));
		int j;
		for(j = 0; j < MAX_SUBJECTS_PER_TEACHER; ++j){
			subjects[i][j] = NONE;
		}
	}
	
}

freeTeacherList(struct TeacherList *tl){
	int i;
	for(i = 0; i < MAX_SUBJECTS_PER_TEACHER; ++i){
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
	int i,j;
	
	for(i = 0; i < DAYS*PERIODS_PER_DAY; ++i){
		for(j = 0; j < CLASSROOMS; ++j){
			ind -> genotype[i][j].teacher_id = NONE;
			ind -> genotype[i][j].class_id = NONE;
			ind -> genotype[i][j].room_id = NONE;
		}
	}
	
	for(i = 0; i < CLASSES; ++i){
		for(j = 0; j < PERIODS_PER_CLASS;){
			int rand_period = RANDOM_PERIOD();
			int slot_id = firstFreeSlotID(ind, rand_period);
			if(slot_id == NONE){
				continue;
			}
			else{
				ind -> genotype[rand_period][slot_id].teacher_id = RANDOM_TEACHER();
				ind -> genotype[rand_period][slot_id].class_id = i;
				ind -> genotype[rand_period][slot_id].room_id = slot_id;
				++j;
			}
		}
	}
}

double fitness(struct Individual *ind){
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
}


