/* Convert SDF to ascii
   Usage:
     SDFcvt [-s skip][-n nrecout][-h hdrfile] [-bfv] sdffile name ...

   Not too fancy, but it should handle most SDF files, including
   non-unit arrcnt.  By using -s and -n is is possible to skip
   to a place in the file and examine pieces of files
   that would otherwise be too big to handle. We do NOT read the
   whole file into memory.  We process output lines one at a time,
   reading one value of each of the names into memory at a time.
   Thrashing is possible.  We do nothing to prevent it.

   With -v  it prints a potentially useful pair of header lines
   starting with '#'.

   CHANGELOG:
   6/28/94:  Added the FORMAT_NATIVE option, and made it the default.
       This will write SDF_INTS as ints, SDF_SHORTS as shorts, etc.  It's
       useful for post-processing by programs that can handle mixed types,
       e.g., perl, but not AVS.

   10/18/93: Added the -S (Swapping) flag to write binary data swapped.
   Use limits.h/INT_MAX instead of SunOS values.h/MAXINT.

   8/1/93: add the AT_A_TIME construction.  We no longer process lines
   	one at a time.  Instead, we read up to AT_A_TIME to beat the
	excessive overhead in each SDFrdvecs call.  One of these days
	I'll add a hash table to SDF's symbol-lookup functions.  Until
	then, reading one element at a time will be extremely slow.

   8/20/93: added the option to produce output in binary.  This is MUCH
   	faster.  With a little luck, perl (and presumably other programs)
	can read it MUCH faster	too (perl uses read and unpack).

*/
#define AT_A_TIME 512		/* read this many values at a time */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "SDF.h"
#include "protos.h"
#include "byteswap.h"
#include "error.h"
#include "Malloc.h"

#ifndef __SUN5__
extern int getopt(int, char**, const char *);
extern char *optarg;
extern int optind, opterr;
#endif

#define FORMAT_FLOAT 1
#define FORMAT_DOUBLE 2
#define FORMAT_INT 3
#define FORMAT_NATIVE 4

int binary = 0;
int format = FORMAT_NATIVE;
int swap_out = 0;

int fwrite_maybe_swap(void *p, size_t sz, size_t n, FILE *fp){
    if( swap_out ) 
	Byteswap(sz, n, p, p);
    return fwrite(p, sz, n, fp);
}

void output_char(short c){
    float f;
    double d;
    int i;
    
    if( binary ){
	switch(format){
	case FORMAT_DOUBLE:
	    d = c;
	    fwrite_maybe_swap(&d, sizeof(d), 1, stdout);
	    break;
	case FORMAT_FLOAT:
	    f = c;
	    fwrite_maybe_swap(&f, sizeof(f), 1, stdout);
	    break;
	case FORMAT_INT:
	    i = c;
	    fwrite_maybe_swap(&i, sizeof(i), 1, stdout);
	    break;
	case FORMAT_NATIVE:
	    fwrite_maybe_swap(&c, sizeof(c), 1, stdout);
	    break;
	}
    }else{
	printf("%d", c);
    }
}

void output_short(short s){
    float f;
    double d;
    int i;
    
    if( binary ){
	switch(format){
	case FORMAT_DOUBLE:
	    d = s;
	    fwrite_maybe_swap(&d, sizeof(d), 1, stdout);
	    break;
	case FORMAT_FLOAT:
	    f = s;
	    fwrite_maybe_swap(&f, sizeof(f), 1, stdout);
	    break;
	case FORMAT_INT:
	    i = s;
	    fwrite_maybe_swap(&i, sizeof(i), 1, stdout);
	    break;
	case FORMAT_NATIVE:
	    fwrite_maybe_swap(&s, sizeof(s), 1, stdout);
	    break;
	}
    }else{
	printf("%d", s);
    }
}

void output_int(int i){
    float f;
    double d;

    if( binary ){
	switch(format){
	case FORMAT_DOUBLE:
	    d = i;
	    fwrite_maybe_swap(&d, sizeof(d), 1, stdout);
	    break;
	case FORMAT_FLOAT:
	    f = i;
	    fwrite_maybe_swap(&f, sizeof(f), 1, stdout);
	    break;
	case FORMAT_NATIVE:
	case FORMAT_INT:
	    fwrite_maybe_swap(&i, sizeof(i), 1, stdout);
	    break;
	}
    }else{
	printf("%d", i);
    }
}

void output_long(long l){
    float f;
    double d;

    if( binary ){
	switch(format){
	case FORMAT_DOUBLE:
	    d = l;
	    fwrite_maybe_swap(&d, sizeof(d), 1, stdout);
	    break;
	case FORMAT_FLOAT:
	    f = l;
	    fwrite_maybe_swap(&f, sizeof(f), 1, stdout);
	    break;
	case FORMAT_NATIVE:
	case FORMAT_INT:
	    fwrite_maybe_swap(&l, sizeof(l), 1, stdout);
	    break;
	}
    }else{
	printf("%ld", l);
    }
}

void output_int64(int64_t i64){
    float f;
    double d;
    
    if( binary ){
	switch(format){
	case FORMAT_DOUBLE:
	    d = i64;
	    fwrite_maybe_swap(&d, sizeof(d), 1, stdout);
	    break;
	case FORMAT_FLOAT:
	    f = i64;
	    fwrite_maybe_swap(&f, sizeof(f), 1, stdout);
	    break;
	case FORMAT_NATIVE:
	case FORMAT_INT:
	    fwrite_maybe_swap(&i64, sizeof(i64), 1, stdout);
	    break;
	}
    }else{
#if __WORDSIZE==64
	printf("%ld", i64);
#else
	printf("%lld", i64);
#endif
    }
}



void output_float(float f){
    double d;
    int i;

    if( binary ){
	switch(format){
	case FORMAT_DOUBLE:
	    d = f;
	    fwrite_maybe_swap(&d, sizeof(d), 1, stdout);
	    break;
	case FORMAT_FLOAT:
	case FORMAT_NATIVE:
	    fwrite_maybe_swap(&f, sizeof(f), 1, stdout);
	    break;
	case FORMAT_INT:
	    i = f;
	    fwrite_maybe_swap(&i, sizeof(i), 1, stdout);
	    break;
	}
    }else{
	switch(format){
	case FORMAT_DOUBLE:
	    printf("%.15e", f);
	    break;
	case FORMAT_FLOAT:
	case FORMAT_NATIVE:
	    printf("%.6e", f);
	    break;
	case FORMAT_INT:
	    printf("%.0f", f);
	    break;
	}
    }
}

void output_double(double d){
    float f;
    int i;

    if( binary ){
	switch(format){
	case FORMAT_NATIVE:
	case FORMAT_DOUBLE:
	    fwrite_maybe_swap(&d, sizeof(d), 1, stdout);
	    break;
	case FORMAT_FLOAT:
	    f = d;
	    fwrite_maybe_swap(&f, sizeof(f), 1, stdout);
	    break;
	case FORMAT_INT:
	    i = d;
	    fwrite_maybe_swap(&i, sizeof(i), 1, stdout);
	    break;
	}
    }else{
	switch(format){
	case FORMAT_DOUBLE:
	case FORMAT_NATIVE:
	    printf("%.15e", d);
	    break;
	case FORMAT_FLOAT:
	    printf("%.6e", d);
	    break;
	case FORMAT_INT:
	    printf("%.0f", d);
	    break;
	}
    }
}

void output_string(const char *s){
    fwrite(s, strlen(s), 1, stdout);
}

void output_separator(void){
    if( !binary )
	putchar(' ');
}

void output_record_terminator(void){
    if( !binary )
	putchar('\n');
}

int main(int argc, char **argv){
    int c;
    char *fname, *hdr;
    SDF *sdfp;
    int64_t n, start;
    int nnames, nn, i, j, k, kk;
    int at_a_time;
    int64_t nrecs, minnrecs, maxnrecs;
    char *newname;
    char **names;
    void **addrs;
    int *nread;
    int *acnt;
    int *strides;
    enum SDF_type_enum *types, t;
    int verbose;

    hdr = NULL;
    start = 0;
    n = 0;
    verbose = 0;
    at_a_time = 0;
    while((c=getopt(argc, argv, "a:h:s:n:vbfdiS")) != -1)
	switch(c){
	case 'h':
	    hdr = optarg;
	    break;
	case 's':
	    start = atoll(optarg);
	    break;
	case 'n':
	    n = atoll(optarg);
	    break;
	case 'a':
	    at_a_time = atoi(optarg);
	    break;
	case 'v':
	    verbose = 1;
	    break;
	case 'b':
	    binary = 1;
	    break;
	case 'f':
	    format = FORMAT_FLOAT;
	    break;
	case 'S':
	    swap_out = 1;
	    break;
	case 'd':
	    format = FORMAT_DOUBLE;
	    break;
	case 'i':
	    format = FORMAT_INT;
	    break;
	case '?':
	    fprintf(stderr, "Usage: %s [-s start][-n number][-a at-a-time][-v<erbose>][-h hdrfile][-b<inary>][-f<loat>][-d<ouble>][-i<nt>][-S<wap>] sdffile name ...\n", argv[0]);
	}

    fname = argv[optind];
    optind++;
    sdfp = SDFopen(hdr, fname);
    if( sdfp == NULL ){
	Error("Could not SDFopen(%s, %s): %s\n", 
	      hdr?hdr:"<null>", 
	      fname?fname:"<null>", 
	      SDFerrstring);
    }
    if( at_a_time == 0 )
	at_a_time = AT_A_TIME;
    
    nnames = argc - optind;
    names = Malloc(nnames*sizeof(char*));
    addrs = Malloc(nnames*sizeof(void*));
    nread = Malloc(nnames*sizeof(int));
    types = Malloc(nnames*sizeof(enum SDF_type_enum));
    strides = Malloc(nnames*sizeof(int));
    acnt = Malloc(nnames*sizeof(int));
    
    nn = 0;
    minnrecs = INT64_MAX;
    maxnrecs = 0;
    for(; optind < argc; optind++){
	newname = argv[optind];
	if( !SDFhasname(newname, sdfp) ){
	    Error("%s not in %s\n", newname, fname);
	}
	t = SDFtype(newname, sdfp);
	acnt[nn] = SDFarrcnt(newname, sdfp);
	addrs[nn] = Malloc(SDFtype_sizes[t] * acnt[nn] * at_a_time);
	strides[nn] = SDFtype_sizes[t]*acnt[nn];
	names[nn] = newname;
	types[nn] = t;
	nrecs = SDFnrecs(newname, sdfp);
	if( nrecs < minnrecs )
	    minnrecs = nrecs;
	if( nrecs > maxnrecs )
	    maxnrecs = nrecs;
	nread[nn] = at_a_time;
	nn++;
    }
    if( n==0 ){
	n = maxnrecs - start;
    }
    if( start+n > minnrecs ){
	Error("Not enough records to start at %ld\n", start);
    }

    if( start ){
	for(j=0; j<nnames; j++){
	    SDFseek(names[j], start, SDF_SEEK_SET, sdfp);
	}
    }

    if( verbose ){
	putchar('#');
	for(i=0; i<argc; i++){
	    printf("%s%c", argv[i], (i<argc-1)?' ':'\n');
	}
	putchar('#');
	for(j=0; j<nnames; j++){
	    printf("%s%c", names[j], (j<nnames-1)?' ':'\n');
	}
    }
	    
    while( n>0 ){
	if( n < at_a_time ){
	    at_a_time = n;
	    for(j=0; j<nnames; j++){
		nread[j] = at_a_time;
	    }
	}
	if( SDFrdvecsarr(sdfp, nnames, names, nread, addrs, strides) )
	    Error("SDFrdvecsarr failed n=%ld: %s\n", n, SDFerrstring); 
	n -= at_a_time;
	for(i=0; i<at_a_time; i++){
	    for(j=0; j<nnames; j++){
		for(k=0; k<acnt[j]; k++){
		    kk = acnt[j]*i + k;
		    if(j||k)
			output_separator();
		    switch(types[j]){
		    case SDF_CHAR:
			output_char( ((char *)addrs[j])[kk] );
			break;
		    case SDF_SHORT:
			output_short( ((short *)addrs[j])[kk] );
			break;
		    case SDF_INT:
			output_int( ((int *)addrs[j])[kk] );
			break;
		    case SDF_LONG:
			output_long( ((long *)addrs[j])[kk] );
			break;
		    case SDF_INT64:
			output_int64( ((int64_t *)addrs[j])[kk] );
			break;
		    case SDF_FLOAT:
			output_float( ((float *)addrs[j])[kk] );
			break;
		    case SDF_DOUBLE:
			output_double( ((double *)addrs[j])[kk] );
			break;
		    case SDF_STRING:
			output_string((char *)addrs[j]);
			break;
		    default:
			Error("Unrecognized type\n");
		    }
		}
	    }
	    output_record_terminator();
	}
    }
    exit(0);
}

