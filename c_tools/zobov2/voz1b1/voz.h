#ifndef __VOZ_H
#define __VOZ_H

#define MAXVERVER 100000
#define NGUARD 84 /*Actually, the number of SPACES between guard points
##define NGUARD 42 /*Actually, the number of SPACES between guard points
		    in each dim */

typedef int pid_t;

typedef struct Partadj {
  pid_t nadj;
  pid_t *adj;
} PARTADJ;

#ifdef __cplusplus
extern "C" {
#endif

int delaunadj (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs);
int vorvol (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, float *vol);
int posread(char *posfile, float ***p, float fact);

#ifdef __cplusplus
}
#endif


#endif
