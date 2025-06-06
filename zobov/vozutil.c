#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libqhull/qhull_a.h"
#include "voz.h"

#define FOREACHvertex2_(vertices) FOREACHsetelement_(vertexT, vertices2,vertex2)

typedef struct {
  double x, y;
  double angle; // 角度を格納
} Point;

/* Shoelaceの公式で凸包の面積を計算 */
double polygon_area(Point *points, int n) {
  double area = 0.0;
  for (int i = 0; i < n; i++) {
      int j = (i + 1) % n;
      area += points[i].x * points[j].y - points[i].y * points[j].x;
  }
  return fabs(area) / 2.0;
}

/* 角度に基づいてソートするための比較関数 */
int compare_angles(const void *a, const void *b) {
  double angle_a = ((Point *)a)->angle;
  double angle_b = ((Point *)b)->angle;
  return (angle_a > angle_b) - (angle_a < angle_b);  // 昇順ソート
}

int compar(const void * n1, const void * n2) {
  int i1,i2;

  i1 = *(int *)n1;
  i2 = *(int *)n2;
  return 2*(i1 > i2) - 1 + (i1 == i2);
}


/* Finds Delaunay adjacencies of a set of points */
int delaunadj (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs) {

  int dim= 3;	            /* dimension of points */
  boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
  char flags[250];          /* option flags for qhull, see qh_opt.htm */
  FILE *outfile= stdout;    /* output from qh_produce_output()
			       use NULL to skip qh_produce_output() */
  FILE *errfile= stderr;    /* error messages from qhull code */
  int exitcode;             /* 0 if no error from qhull */
  int curlong, totlong;	    /* memory remaining after qh_memfreeshort */
  int i, ver, count;
  int numfacets, numsimplicial, numridges, totneighbors, numneighbors,
    numcoplanars, numtricoplanars;
  setT *vertices, *vertices2, *vertex_points, *coplanar_points;
  vertexT *vertex, **vertexp;
  vertexT *vertex2, **vertex2p;
  int vertex_i, vertex_n;
  facetT *facet, **facetp, *neighbor, **neighborp;
  pointT *point, **pointp;
  int numdiv;

  PARTADJ adjst;

  int errorReported = 0;

  adjst.adj = (int *)malloc(MAXVERVER*sizeof(int));
  if (adjst.adj == NULL) {
    printf("Unable to allocate adjst.adj\n");
    exit(0);
  }

  /* Delaunay triangulation*/
  sprintf (flags, "qhull s d Qt");
  exitcode= qh_new_qhull (dim, nvpall, points, ismalloc,
                      flags, outfile, errfile);


  if (!exitcode) {                  /* if no error */
    /* 'qh facet_list' contains the convex hull */

    /* From qh_printvneighbors */
    qh_countfacets(qh facet_list, NULL, 0, &numfacets, &numsimplicial,
		   &totneighbors, &numridges, &numcoplanars, &numtricoplanars);
    qh_vertexneighbors();
    vertices= qh_facetvertices (qh facet_list, NULL, 0);
    vertex_points= qh_settemp (nvpall);
    coplanar_points= qh_settemp (nvpall);
    qh_setzero (vertex_points, 0, nvpall);
    qh_setzero (coplanar_points, 0, nvpall);
    FOREACHvertex_(vertices)
      qh_point_add (vertex_points, vertex->point, vertex);
    FORALLfacet_(qh facet_list) {
      FOREACHpoint_(facet->coplanarset)
	qh_point_add (coplanar_points, point, facet);
    }
    ver = 0;
    FOREACHvertex_i_(vertex_points) {
      (*adjs)[ver].nadj = 0;
      if (vertex) {
	/* Count the neighboring vertices, check that all are real
	   neighbors */
	adjst.nadj = 0;
	FOREACHneighbor_(vertex) {
	  if ((*adjs)[ver].nadj > -1) {
	    if (neighbor->visitid) {
	      vertices2 = neighbor->vertices;
	      FOREACHvertex2_(vertices2) {
		if (ver != qh_pointid(vertex2->point)) {
		  adjst.adj[adjst.nadj] = (int)qh_pointid(vertex2->point);
		  adjst.nadj ++;
		  if (adjst.nadj > MAXVERVER) {
		    printf("Increase MAXVERVER to at least %d!\n",adjst.nadj);
		    exit(0);
		  }
		}
	      }
	    } else {
	      printf(" %d",ver);
	      (*adjs)[ver].nadj = -1; /* There are unreal vertices here */
	    }
	  }
	}
      } else (*adjs)[ver].nadj = -2;

      /* Enumerate the unique adjacencies*/
      if (adjst.nadj >= 4) {
	qsort((void *)adjst.adj, adjst.nadj, sizeof(int), &compar);
	count = 1;

	for (i=1; i<adjst.nadj; i++)
	  if (adjst.adj[i] != adjst.adj[i-1]) {
	    if (adjst.adj[i] >= nvpbuf && !errorReported) {
        errorReported = 1;
	      printf("Guard point encountered.  Increase border and/or nguard.\n");
	      printf("P:(%f,%f,%f), G: (%f,%f,%f)\n",points[3*ver],points[3*ver+1],points[3*ver+2],
		     points[3*adjst.adj[i]],points[3*adjst.adj[i]+1],points[3*adjst.adj[i]+2]);
	    }
	    count++;
	  }
	(*adjs)[ver].adj = (int *)malloc(count*sizeof(int));
	if ((*adjs)[ver].adj == NULL) {
	  printf("Unable to allocate (*adjs)[ver].adj\n");
	  exit(0);
	}
	(*adjs)[ver].adj[0] = adjst.adj[0];
	count = 1;
	for (i=1; i<adjst.nadj; i++)
	  if (adjst.adj[i] != adjst.adj[i-1]) {
	    (*adjs)[ver].adj[count] = adjst.adj[i];
	    count++;
	  }
	(*adjs)[ver].nadj = count;
      } else {
	printf("Number of adjacencies %d < 4, particle %d -> %d\n",adjst.nadj,ver,ver);
	exit(0);
      }
      ver++;
      if (ver == nvp) break;
    }
    qh_settempfree (&coplanar_points);
    qh_settempfree (&vertex_points);
    qh_settempfree (&vertices);
  }

   /* 隣接点とその座標を表示 */
  printf("\n--- Delaunay Neighbors ---\n");
  for (int ver = 0; ver < nvp; ver++) {
    printf("Point %d (%.6f, %.6f): ", ver, points[2*ver], points[2*ver+1]);
    for (int i = 0; i < (*adjs)[ver].nadj; i++) {
      int adj_index = (*adjs)[ver].adj[i];
      printf("[%d: (%.6f, %.6f)] ", adj_index, points[2*adj_index], points[2*adj_index+1]);
    }
    printf("\n");
  }

  qh_freeqhull(!qh_ALL);                 /* free long memory */
  qh_memfreeshort (&curlong, &totlong);  /* free short memory and memory allocator */
  if (curlong || totlong)
    fprintf (errfile, "qhull internal warning (delaunadj): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
  free(adjst.adj);
  return exitcode;
}

/* Finds Delaunay adjacencies of a set of points */
int delaunadj_2D (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs) {

  int dim= 2;	            /* dimension of points */
  boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
  char flags[250];          /* option flags for qhull, see qh_opt.htm */
  FILE *outfile= stdout;    /* output from qh_produce_output()
			       use NULL to skip qh_produce_output() */
  FILE *errfile= stderr;    /* error messages from qhull code */
  int exitcode;             /* 0 if no error from qhull */
  int curlong, totlong;	    /* memory remaining after qh_memfreeshort */
  int i, ver, count;
  int numfacets, numsimplicial, numridges, totneighbors, numneighbors,
    numcoplanars, numtricoplanars;
  setT *vertices, *vertices2, *vertex_points, *coplanar_points;
  vertexT *vertex, **vertexp;
  vertexT *vertex2, **vertex2p;
  int vertex_i, vertex_n;
  facetT *facet, **facetp, *neighbor, **neighborp;
  pointT *point, **pointp;
  int numdiv;

  PARTADJ adjst;

  int errorReported = 0;

  adjst.adj = (int *)malloc(MAXVERVER*sizeof(int));
  if (adjst.adj == NULL) {
    printf("Unable to allocate adjst.adj\n");
    exit(0);
  }

  /* Delaunay triangulation*/
  sprintf (flags, "qhull d Qz");
  exitcode= qh_new_qhull (dim, nvpall, points, ismalloc,
                      flags, outfile, errfile);


  if (!exitcode) {                  /* if no error */
    /* 'qh facet_list' contains the convex hull */

    /* From qh_printvneighbors */
    qh_countfacets(qh facet_list, NULL, 0, &numfacets, &numsimplicial,
		   &totneighbors, &numridges, &numcoplanars, &numtricoplanars);
    qh_vertexneighbors();
    vertices= qh_facetvertices (qh facet_list, NULL, 0);
    vertex_points= qh_settemp (nvpall);
    coplanar_points= qh_settemp (nvpall);
    qh_setzero (vertex_points, 0, nvpall);
    qh_setzero (coplanar_points, 0, nvpall);
    FOREACHvertex_(vertices)
      qh_point_add (vertex_points, vertex->point, vertex);
    FORALLfacet_(qh facet_list) {
      FOREACHpoint_(facet->coplanarset)
	qh_point_add (coplanar_points, point, facet);
    }
    ver = 0;
    FOREACHvertex_i_(vertex_points) {
      (*adjs)[ver].nadj = 0;
      if (vertex) {
	/* Count the neighboring vertices, check that all are real
	   neighbors */
	adjst.nadj = 0;
	FOREACHneighbor_(vertex) {
	  if ((*adjs)[ver].nadj > -1) {
	    if (neighbor->visitid) {
	      vertices2 = neighbor->vertices;
	      FOREACHvertex2_(vertices2) {
		if (ver != qh_pointid(vertex2->point)) {
		  adjst.adj[adjst.nadj] = (int)qh_pointid(vertex2->point);
		  adjst.nadj ++;
		  if (adjst.nadj > MAXVERVER) {
		    printf("Increase MAXVERVER to at least %d!\n",adjst.nadj);
		    exit(0);
		  }
		}
	      }
	    } else {
	      printf(" %d",ver);
	      (*adjs)[ver].nadj = -1; /* There are unreal vertices here */
	    }
	  }
	}
      } else (*adjs)[ver].nadj = -2;

      /* Enumerate the unique adjacencies*/
      if (adjst.nadj >= 3) {
	qsort((void *)adjst.adj, adjst.nadj, sizeof(int), &compar);
	count = 1;

	for (i=1; i<adjst.nadj; i++)
	  if (adjst.adj[i] != adjst.adj[i-1]) {
	    if (adjst.adj[i] >= nvpbuf && !errorReported) {
        errorReported = 1;
	      printf("Guard point encountered.  Increase border and/or nguard.\n");
	      printf("P:(%f,%f), G: (%f,%f)\n",points[2*ver],points[2*ver+1],
		     points[2*adjst.adj[i]],points[2*adjst.adj[i]+1]);
	    }
	    count++;
	  }
	(*adjs)[ver].adj = (int *)malloc(count*sizeof(int));
	if ((*adjs)[ver].adj == NULL) {
	  printf("Unable to allocate (*adjs)[ver].adj\n");
	  exit(0);
	}
	(*adjs)[ver].adj[0] = adjst.adj[0];
	count = 1;
	for (i=1; i<adjst.nadj; i++)
	  if (adjst.adj[i] != adjst.adj[i-1]) {
	    (*adjs)[ver].adj[count] = adjst.adj[i];
	    count++;
	  }
	(*adjs)[ver].nadj = count;
      } else {
	printf("Number of adjacencies %d < 4, particle %d -> %d\n",adjst.nadj,ver,ver);
	exit(0);
      }
      ver++;
      if (ver == nvp) break;
    }
    qh_settempfree (&coplanar_points);
    qh_settempfree (&vertex_points);
    qh_settempfree (&vertices);
  }

  /* 隣接点とその座標を表示 */
  printf("\n--- Delaunay Neighbors ---\n");
  for (int ver = 0; ver < nvp; ver++) {
    if(ver%100 == 0){
      printf("Point %d (%.6f, %.6f): ", ver, points[2*ver], points[2*ver+1]);
      for (int i = 0; i < (*adjs)[ver].nadj; i++) {
        int adj_index = (*adjs)[ver].adj[i];
        printf("[%d: (%.6f, %.6f)] ", adj_index, points[2*adj_index], points[2*adj_index+1]);
      }
      printf("\n");
    }
  }

  qh_freeqhull(!qh_ALL);                 /* free long memory */
  qh_memfreeshort (&curlong, &totlong);  /* free short memory and memory allocator */
  if (curlong || totlong)
    fprintf (errfile, "qhull internal warning (delaunadj): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
  free(adjst.adj);
  return exitcode;
}

/* Calculates the Voronoi volume from a set of Delaunay adjacencies */
int vorvol (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, float *vol) {
  int dim= 3;	            /* dimension of points */
  boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
  char flags[250];          /* option flags for qhull, see qh_opt.htm */
  FILE *outfile= NULL;    /* output from qh_produce_output()
			       use NULL to skip qh_produce_output() */
  FILE *errfile= stderr;    /* error messages from qhull code */
  int exitcode;             /* 0 if no error from qhull */
  facetT *facet;	    /* set by FORALLfacets */
  int curlong, totlong;	    /* memory remaining after qh_memfreeshort */

  coordT *point, *normp, *coordp, *feasiblep, *deladj;
  int i, j, k;
  boolT zerodiv;
  float runsum;
  char region;
  /*coordT *points;
    pointT *intpoints;*/

  /* make point array from adjacency coordinates (add offset)*/
  /*points = (coordT *)malloc(4*numpoints*sizeof(coordT));
  if (points == NULL) {
    printf("Unable to allocate points\n");
    exit(0);
    }*/
  for (i=0; i<numpoints; i++) {
    runsum = 0.;
    deladj = deladjs + 3*i;
    point = points + 4*i;
    for (j=0;j<3;j++) {
      runsum += deladj[j]*deladj[j];
      point[j] = deladj[j];
    }
    point[3] = -0.5*runsum;
  }
  sprintf (flags, "qhull H0");

  exitcode= qh_new_qhull (4, numpoints, points, ismalloc,
			  flags, outfile, errfile);

  numpoints = 0;
  if (!exitcode) {                  /* if no error */
    FORALLfacets {
      numpoints++;
    }
    /* Now we know how many points */
    /*intpoints = (pointT *)malloc(dim*numpoints*sizeof(pointT));
    if (intpoints == NULL) {
    printf("Unable to allocate intpoints\n");
    exit(0);
    }*/

    j = 0;
    FORALLfacets {
      if (!qh feasible_point) {
	fprintf (stdout, "qhull input error (qh_printafacet): option 'Fp' needs qh feasible_point\n");
	qh_errexit( qh_ERRinput, NULL, NULL);
      }
      point= coordp= intpoints + j*3;
      j++;
      normp= facet->normal;
      feasiblep= qh feasible_point;
      if (facet->offset < -qh MINdenom) {
	for (k= qh hull_dim; k--; )
	  *(coordp++)= (*(normp++) / - facet->offset) + *(feasiblep++);
    //printf("nx %f, ny %f, nz %f, nw %f, offset %f\n" , facet->normal[0], facet->normal[1], facet->normal[2], facet->normal[3], -facet->offset);
      }else {
	for (k= qh hull_dim; k--; ) {
	  *(coordp++)= qh_divzero (*(normp++), facet->offset, qh MINdenom_1,
				   &zerodiv) + *(feasiblep++);
	  if (zerodiv) {
	    qh_memfree (point, qh normal_size);
	    printf("LABELprintinfinite\n");
	    exit(0);
	  }
	}
      }
    }
  }
  qh_freeqhull (!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);

  /* Now we calculate the volume */
  sprintf (flags, "qhull FA");
  exitcode= qh_new_qhull (dim, numpoints, intpoints, ismalloc,
			  flags, outfile, errfile);

  qh_getarea(qh facet_list);
  *vol = qh totvol;

  qh_freeqhull (!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  if (curlong || totlong)
    fprintf (errfile, "qhull internal warning (vorvol): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
  /*free(points); free(intpoints);*/

  return exitcode;
}

int vorvol_2D (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, float *area) {
  int dim = 2;              /* 2D */
  boolT ismalloc = False;
  char flags[250];
  FILE *outfile = NULL;
  FILE *errfile = stderr;
  int exitcode;
  facetT *facet;
  int curlong, totlong;
  
  coordT *point, *deladj, *coordp, *feasiblep, *normp;
  boolT zerodiv;
  int i, j, k;
  float runsum;

  /* Adjacency coordinates を2次元に変換 */
  for (i = 0; i < numpoints; i++) {
    runsum = 0.;
    deladj = deladjs + 2 * i;
    point = points + 3 * i;  /* Homogeneous座標系を考慮 */
    for (j = 0; j < 2; j++) {
      runsum += (deladj[j]*100.0) * (deladj[j]*100.0);
      point[j] = deladj[j]*100.0;
    }
    point[2] = -0.5 * runsum; /* Homogeneous座標 */
  }

  //printf("-------------print points--------------------------\n");

  //for (i = 0; i < numpoints; i++) {
    //printf("[%f, %f],\n", points[3*i], points[3*i+1]);
  //}

  sprintf(flags, "qhull H0");

  exitcode = qh_new_qhull(3, numpoints, points, ismalloc, flags, outfile, errfile);

  numpoints = 0;
  if (!exitcode) {
    FORALLfacets {
      numpoints++;
    }

    j = 0;
    //printf("--------------normal--------------------------\n");
    FORALLfacets {
      if (!qh feasible_point) {
        fprintf (stdout, "qhull input error (qh_printafacet): option 'Fp' needs qh feasible_point\n");
        qh_errexit( qh_ERRinput, NULL, NULL);
            }
            point= coordp= intpoints + j*2;
            j++;
            normp= facet->normal;
            feasiblep= qh feasible_point;
            if (facet->offset < -qh MINdenom) {
        for (k= qh hull_dim; k--; )
          *(coordp++)= (*(normp++)/ - facet->offset) + *(feasiblep++);
            }else {
        for (k= qh hull_dim; k--; ) {
          *(coordp++)= qh_divzero (*(normp++), facet->offset, qh MINdenom_1,
                 &zerodiv) + *(feasiblep++);
          if (zerodiv) {
            qh_memfree (point, qh normal_size);
            printf("LABELprintinfinite\n");
            exit(0);
          }
        }
            }
      //printf("[%f, %f],\n" , intpoints[2*(j-1)+0], intpoints[2*(j-1)+1]);
    }
  //}
  }

  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);

  /* 面積計算 */
  sprintf(flags, "qhull FA");
  exitcode = qh_new_qhull(dim, numpoints, intpoints, ismalloc, flags, outfile, errfile);

  //qh_getarea(qh facet_list);
  //*area = qh totarea/10000.0/2.0;

  /* 凸包の頂点を取得 */
  int num_hull_points = qh num_vertices;
  Point *hull_points = malloc(num_hull_points * sizeof(Point));
  double centroid_x = 0, centroid_y = 0;
  
  vertexT *vertex;
  i=0;
  FORALLvertices {
      hull_points[i].x = vertex->point[0];
      hull_points[i].y = vertex->point[1];
      centroid_x += hull_points[i].x;
      centroid_y += hull_points[i].y;
      i++;
  }

  /* 重心を求める */
  centroid_x /= num_hull_points;
  centroid_y /= num_hull_points;

  /* 各点の角度を計算 */
  for (int i = 0; i < num_hull_points; i++) {
      hull_points[i].angle = atan2(hull_points[i].y - centroid_y, hull_points[i].x - centroid_x);
  }

  /* 角度でソート（反時計回りに並び替え） */
  qsort(hull_points, num_hull_points, sizeof(Point), compare_angles);

  /* 並び替えた凸包の点を表示 */
  //printf("Sorted Convex Hull Points:\n");
  //for (int i = 0; i < num_hull_points; i++) {
      //printf("(%f, %f)\n", hull_points[i].x, hull_points[i].y);
  //}

  /* 面積計算 */

  *area = polygon_area(hull_points, num_hull_points)/10000.0;
  printf("volume %f\n" , *area*10000.0);

   /* メモリ解放 */
  free(hull_points);
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    fprintf(errfile, "qhull internal warning (vorarea): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);

  return exitcode;
}

int old_vorvol_2D(coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, float *area) {
  int dim = 2; // 2Dに変更
  boolT ismalloc = False;
  char flags[250];
  FILE *outfile = NULL;
  FILE *errfile = stderr;
  int exitcode;
  facetT *facet;
  int curlong, totlong;
  coordT *point, *deladj;
  int i, j;
  float runsum;
  int num_facets = 0;  // 修正点

  // Adjacency coordinates を2次元に変換
  for (i = 0; i < numpoints; i++) {
    runsum = 0.;
    deladj = deladjs + 2 * i;
    point = points + 2 * i;  // 3Dではなく2D

    for (j = 0; j < 2; j++) {
      runsum += deladj[j] * deladj[j];
      point[j] = deladj[j];
    }
  }

  // Qhull 実行
  sprintf(flags, "qhull v Qbb Qx FA");
  exitcode = qh_new_qhull(dim, numpoints, points, ismalloc, flags, outfile, errfile);
  
  if (!exitcode) {
    FORALLfacets {
      num_facets++;
    }
    j = 0;
    FORALLfacets {
      point = intpoints + j * 2;
      j++;
      for (i = 0; i < 2; i++) {
        point[i] = facet->normal[i] / -facet->offset;
      }
    }
  }

  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);

  // 面積計算
  if (num_facets == 0) {
    fprintf(errfile, "qhull error: no Voronoi facets found\n");
    return -1;
  }

  sprintf(flags, "qhull v Qx Qz FA");
  exitcode = qh_new_qhull(dim, num_facets, intpoints, ismalloc, flags, outfile, errfile);

  if (!exitcode) {
    qh_getarea(qh facet_list);
    *area = qh totarea;
  }

  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);

  if (curlong || totlong)
    fprintf(errfile, "qhull internal warning (vorarea): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);

  return exitcode;
}

