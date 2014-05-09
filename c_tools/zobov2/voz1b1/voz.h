/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/zobov2/voz1b1/voz.h
    Copyright (C) 2010-2014 Guilhem Lavaux
    Copyright (C) 2011-2014 P. M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/
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
