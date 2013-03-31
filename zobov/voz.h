#define MAXVERVER 100000
#define NGUARD 84 /*Actually, the number of SPACES between guard points
##define NGUARD 42 /*Actually, the number of SPACES between guard points
		    in each dim */

typedef int pid_t;

typedef struct Partadj {
  pid_t nadj;
  pid_t *adj;
} PARTADJ;
