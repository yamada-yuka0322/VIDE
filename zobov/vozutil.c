#include "libqhull/qhull_a.h"
#include "voz.h"

#define FOREACHvertex2_(vertices) FOREACHsetelement_(vertexT, vertices2,vertex2)

/*char qh_version[] = "user_eg 3.1 2001/10/04";  */

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
  qh_freeqhull(!qh_ALL);                 /* free long memory */
  qh_memfreeshort (&curlong, &totlong);  /* free short memory and memory allocator */
  if (curlong || totlong)
    fprintf (errfile, "qhull internal warning (delaunadj): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
  free(adjst.adj);
  return exitcode;
}

/* Finds Delaunay adjacencies for 2D point distribution */
int delaunadj_2D (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs) {
  int dim = 2;               /* 2Dに変更 */
  boolT ismalloc = False;
  char flags[250];
  FILE *outfile = stdout;
  FILE *errfile = stderr;
  int exitcode;
  int curlong, totlong;
  int i, ver, count;
  int numfacets, numsimplicial, numridges, totneighbors, numneighbors,
      numcoplanars, numtricoplanars;
  setT *vertices, *vertex_points, *coplanar_points;
  vertexT *vertex, **vertexp;
  facetT *facet, **facetp, *neighbor, **neighborp;
  pointT *point, **pointp;
  PARTADJ adjst;
  int errorReported = 0;

  adjst.adj = (int *)malloc(MAXVERVER * sizeof(int));
  if (adjst.adj == NULL) {
    printf("Unable to allocate adjst.adj\n");
    exit(0);
  }

  /* 2D Delaunay triangulation */
  sprintf(flags, "qhull s Qt");  // 2D用に修正
  exitcode = qh_new_qhull(dim, nvpall, points, ismalloc, flags, outfile, errfile);

  if (!exitcode) {
    qh_countfacets(qh facet_list, NULL, 0, &numfacets, &numsimplicial,
                   &totneighbors, &numridges, &numcoplanars, &numtricoplanars);
    qh_vertexneighbors();
    vertices = qh_facetvertices(qh facet_list, NULL, 0);
    vertex_points = qh_settemp(nvpall);
    coplanar_points = qh_settemp(nvpall);
    qh_setzero(vertex_points, 0, nvpall);
    qh_setzero(coplanar_points, 0, nvpall);
    FOREACHvertex_(vertices)
      qh_point_add(vertex_points, vertex->point, vertex);
    FORALLfacet_(qh facet_list) {
      FOREACHpoint_(facet->coplanarset)
        qh_point_add(coplanar_points, point, facet);
    }
    ver = 0;
    FOREACHvertex_i_(vertex_points) {
      (*adjs)[ver].nadj = 0;
      if (vertex) {
        adjst.nadj = 0;
        FOREACHneighbor_(vertex) {
          if ((*adjs)[ver].nadj > -1) {
            if (neighbor->visitid) {
              vertices2 = neighbor->vertices;
              FOREACHvertex2_(vertices2) {
                if (ver != qh_pointid(vertex2->point)) {
                  adjst.adj[adjst.nadj] = (int)qh_pointid(vertex2->point);
                  adjst.nadj++;
                  if (adjst.nadj > MAXVERVER) {
                    printf("Increase MAXVERVER to at least %d!\n", adjst.nadj);
                    exit(0);
                  }
                }
              }
            } else {
              printf(" %d", ver);
              (*adjs)[ver].nadj = -1;
            }
          }
        }
      } else (*adjs)[ver].nadj = -2;

      if (adjst.nadj >= 3) { // 2Dでは3以上に変更
        qsort((void *)adjst.adj, adjst.nadj, sizeof(int), &compar);
        count = 1;
        for (i = 1; i < adjst.nadj; i++)
          if (adjst.adj[i] != adjst.adj[i - 1]) {
            if (adjst.adj[i] >= nvpbuf && !errorReported) {
              errorReported = 1;
              printf("Guard point encountered. Increase border and/or nguard.\n");
              printf("P:(%f,%f), G: (%f,%f)\n", points[2 * ver], points[2 * ver + 1],
                     points[2 * adjst.adj[i]], points[2 * adjst.adj[i] + 1]);
            }
            count++;
          }
        (*adjs)[ver].adj = (int *)malloc(count * sizeof(int));
        if ((*adjs)[ver].adj == NULL) {
          printf("Unable to allocate (*adjs)[ver].adj\n");
          exit(0);
        }
        (*adjs)[ver].adj[0] = adjst.adj[0];
        count = 1;
        for (i = 1; i < adjst.nadj; i++)
          if (adjst.adj[i] != adjst.adj[i - 1]) {
            (*adjs)[ver].adj[count] = adjst.adj[i];
            count++;
          }
        (*adjs)[ver].nadj = count;
      } else {
        printf("Number of adjacencies %d < 3, particle %d -> %d\n", adjst.nadj, ver, ver);
        exit(0);
      }
      ver++;
      if (ver == nvp) break;
    }
    qh_settempfree(&coplanar_points);
    qh_settempfree(&vertex_points);
    qh_settempfree(&vertices);
  }
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    fprintf(errfile, "qhull internal warning (delaunadj_2D): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
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

int vorvol_2D (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, float *vol) {
  int dim= 2;               /* 2D の場合の次元 */
  boolT ismalloc= False;    /* qhull がポイントを解放するか */
  char flags[250];          /* qhull のオプション */
  FILE *outfile= NULL;    /* qh_produce_output() の出力 */
  FILE *errfile= stderr;    /* エラーメッセージ用ファイル */
  int exitcode;             /* qhull の実行結果 */
  facetT *facet;            /* FORALLfacets で使用 */
  int curlong, totlong;     /* メモリ解放チェック */

  coordT *point, *normp, *coordp, *feasiblep, *deladj;
  int i, j, k;
  boolT zerodiv;
  float runsum;

  /* points 配列の作成 (2D 用) */
  for (i=0; i<numpoints; i++) {
    runsum = 0.;
    deladj = deladjs + 2*i;
    point = points + 3*i;
    for (j=0; j<2; j++) {
      runsum += deladj[j]*deladj[j];
      point[j] = deladj[j];
    }
    point[2] = -0.5*runsum;
  }
  sprintf (flags, "qhull QJ");

  exitcode= qh_new_qhull (3, numpoints, points, ismalloc, flags, outfile, errfile);

  numpoints = 0;
  if (!exitcode) {
    FORALLfacets {
      numpoints++;
    }

    j = 0;
    FORALLfacets {
      if (!qh feasible_point) {
        fprintf (stdout, "qhull input error: option 'Fp' needs qh feasible_point\n");
        qh_errexit(qh_ERRinput, NULL, NULL);
      }
      point= coordp= intpoints + j*2;
      j++;
      normp= facet->normal;
      feasiblep= qh feasible_point;
      if (facet->offset < -qh MINdenom) {
        for (k= qh hull_dim; k--; )
          *(coordp++)= (*(normp++) / - facet->offset) + *(feasiblep++);
      } else {
        for (k= qh hull_dim; k--; ) {
          *(coordp++)= qh_divzero (*(normp++), facet->offset, qh MINdenom_1, &zerodiv) + *(feasiblep++);
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

  /* 面積計算（体積の代わり） */
  sprintf (flags, "qhull FA");
  exitcode= qh_new_qhull (dim, numpoints, intpoints, ismalloc, flags, outfile, errfile);

  qh_getarea(qh facet_list);
  *vol = qh totarea;

  qh_freeqhull (!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  if (curlong || totlong)
    fprintf (errfile, "qhull internal warning (vorvol): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);

  return exitcode;
}

