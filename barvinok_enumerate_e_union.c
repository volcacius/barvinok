#include <unistd.h>
#include <sys/times.h>
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
#include <util.h>
#include <barvinok.h>
#include "config.h"
#include <time.h>

/* The input of this example program is a polytope in combined
 * data and parameter space followed by two lines indicating
 * the number of existential variables and parameters respectively.
 * The first lines starts with "E ", followed by a number.
 * The second lines starts with "P ", followed by a number.
 * These two lines are (optionally) followed by the names of the parameters.
 * The polytope is in PolyLib notation.
 */

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    0
#else
#define MAXRAYS  600
#endif

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "pip",   no_argument,  0,  'p' },
    { "convert",   no_argument,  0,  'c' },
    { "floor",     no_argument,  0,  'f' },
    { "range",	no_argument,	0,  'r' },
    { "omega_disjoint", no_argument, 0, 'o' },
    { 0, 0, 0, 0 }
};
#endif

typedef struct _Polyhedron_union { 
	Enumeration *pt;
        struct _Polyhedron_union *next;} Polyhedron_union;

Polyhedron *DMUnion(Enumeration *en, unsigned MR)  {
        Enumeration *e1;
        Polyhedron *d;
        e1=en;
        d=e1->ValidityDomain;
         for (e1=en->next; e1; e1=e1->next) {
           d= DomainUnion( d, e1->ValidityDomain, MR);
         }
         return d;
  }


int main(int argc, char **argv)
{
    Polyhedron *A=0;
    Polyhedron *d1=0, *d2=0;
    Polyhedron *d=0;
    Matrix *M=0;
    char **param_name=0;
    int exist, nparam, npol;
    char s[128];
    evalue *EP=0, *_en=0;
    Enumeration *sen=0, *pr=0, *e=0, *en=0, *res=0, *en1=0, *en2=0, *tmp=0;
    Polyhedron_union * Polun=0;
    Polyhedron_union * pu=0;
    int c, ind = 0;
    int range = 0;
    int convert = 0;
    int pip = 0;
    int floor = 0;
    int omega_disjoint = 0;
    int p;
    Polyhedron *lp=0, *lp1=0, *lp1next=0;
    int nlp;
    clock_t time;

    while ((c = getopt_long(argc, argv, "opfcr", options, &ind)) != -1) {
	switch (c) {
	case 'p':
	    pip = 1;
	    break;
	case 'f':
	    floor = 1;
	    break;
	case 'c':
	    convert = 1;
	    break;
	case 'r':
	    range = 1;
	    break;
        case 'o':
            omega_disjoint = 1;
	    break;
	}
    }

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "N %d", &npol)<1))
	fgets(s, 128, stdin);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "E %d", &exist)<1))
	fgets(s, 128, stdin);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "P %d", &nparam)<1))
	fgets(s, 128, stdin);
    param_name = Read_ParamNames(stdin, nparam);

    time = clock();
    for(p=0; p<npol; p++) {
      Polyhedron* read_p;
      M = Matrix_Read();
      read_p = Constraints2Polyhedron(M, MAXRAYS);
      read_p->next = A;
      A = read_p;
      Matrix_Free(M);

      //Polyhedron_Print(stdout, P_VALUE_FMT, A);
    }
    fprintf(stderr, "Read in %d polytopes in %g seconds.\n",npol, ((double)clock()-time)/CLOCKS_PER_SEC);
    /* compute disjoint union */
    fprintf(stderr, "Starting to compute disjoint union of %d polytopes.\n", npol);
    fflush(stderr);
    time = clock();
    if (!omega_disjoint)
      lp = Disjoint_Domain(A, 0, MAXRAYS);
    else
      lp = A;
    nlp=0;
    for(lp1=lp; lp1; lp1=lp1->next)
      nlp++;
    fprintf(stderr, "The disjoint union is calculated in %g seconds and consists of %d polytopes.\n",
            ((double)clock()-time)/CLOCKS_PER_SEC, nlp);
    fflush(stderr);
   /* 
    printf("exist: %d, nparam: %d\n", exist, nparam);
    if (pip && exist > 0)
	EP = barvinok_enumerate_pip(A, exist, nparam, MAXRAYS);
    else
	EP = barvinok_enumerate_e(A, exist, nparam, MAXRAYS);
*/
 
    nlp=0; 
   	for (lp1=lp ; lp1; lp1=lp1->next)
	{
        fprintf(stderr, "Started enumerating disjoint polytope %d.\n",++nlp);
		fflush(stderr);
        time = clock();
		lp1next = lp1->next;
		lp1->next = NULL;
		//en= Polyhedron_Enumerate(lp1, C, MAXRAYS,NULL);
                if (pip && exist > 0)
                  _en = barvinok_enumerate_pip(lp1, exist, nparam, MAXRAYS);
                else
                  _en = barvinok_enumerate_e(lp1, exist, nparam, MAXRAYS);
		reduce_evalue(_en);
		evalue_combine(_en);
		if (range)
		    evalue_range_reduction(_en);
		//print_evalue(stdout, en, param_name);
		if (floor) {
		    fprintf(stderr, "WARNING: floor conversion not supported\n");
		    evalue_frac2floor(_en);
		    //print_evalue(stdout, _en, param_name);
		} else if (convert) {
		    evalue_mod2table(_en, nparam);
		    //print_evalue(stdout, _en, param_name);
		}
 		lp1->next = lp1next;
		en = partition2enumeration(_en);
		sen= NULL;
		//en = (Enumeration*)malloc(sizeof(Enumeration));
		//en 
		for(e=en;e;e=e->next)
			if(!Degenerate(e))
			{
				pr = (Enumeration  *)malloc(sizeof(Enumeration));
				pr->EP=e->EP;     
				pr->ValidityDomain=e->ValidityDomain;
				pr->next=sen;
				sen=pr;
			}

		if(sen!= NULL)
		{
			pu = (Polyhedron_union  *)malloc(sizeof(Polyhedron_union));
			pu->pt=sen;
			pu->next = Polun;
			Polun = pu;
		}
        fprintf(stderr, "Enumerated disjoint polytope %d in %g seconds.\n",nlp, ((double)clock()-time)/CLOCKS_PER_SEC);
		fflush(stderr);
	}
	if(!Polun)
	{
#ifdef UE_DEBUG
		fprintf(stdout,"Empty Polun\n");
#endif
		return 0;
	}
      
	while(Polun->next != NULL)  {
		res=NULL;
		en1=Polun->pt;
		en2=(Polun->next)->pt;

		d1=DMUnion(en1, MAXRAYS);
		d2=DMUnion(en2, MAXRAYS);

		for(en1=Polun->pt;en1;en1=en1->next)
		{

			for(en2=(Polun->next)->pt;en2;en2=en2->next)
			{
				d = DomainIntersection(en1->ValidityDomain,en2->ValidityDomain,MAXRAYS);
				if( d && !emptyQ(d)&&!IncludeInRes(d,res,MAXRAYS))  {
					evalue ev;
					value_init(ev.d);
					value_assign( ev.d, en2->EP.d );
					if(value_zero_p(ev.d))
						ev.x.p=ecopy(en2->EP.x.p);
					else
					{
						value_init(ev.x.n);
						value_assign( ev.x.n, en2->EP.x.n );
					}

					eadd(&en1->EP,&ev);
					tmp = (Enumeration  *)malloc(sizeof(Enumeration));
					tmp->ValidityDomain =d;
					tmp->EP=ev;
					tmp->next= res;
					res=tmp;
				}
			}
			d=DomainDifference(en1->ValidityDomain,d2 ,MAXRAYS);
			if( d && !emptyQ(d)&&!IncludeInRes(d,res,MAXRAYS))
			{
				tmp = (Enumeration  *)malloc(sizeof(Enumeration));
				tmp->ValidityDomain =d;

				tmp->EP=en1->EP;
				tmp->next= res;
				res=tmp;
			}
		}
		for(en2=(Polun->next)->pt; en2; en2= en2->next)
		{
			d= DomainDifference(en2->ValidityDomain,d1,MAXRAYS);
			if( d && !emptyQ(d)&&!IncludeInRes(d,res,MAXRAYS) )
			{
				tmp = (Enumeration  *)malloc(sizeof(Enumeration));
				tmp->ValidityDomain =d;
				tmp->EP=en2->EP;
				tmp->next= res;
				res=tmp;
			}
		}
	    
		Polun->pt=res;
	        		     
		Polun->next= (Polun->next)->next;
	}
	res=Polun->pt;
		
	Remove_RedundantDomains(&res); 
	//return(res);
	// Print result in res;
        //print_evalue(stdout, res, param_name);
 
   	for (tmp=res ; tmp; tmp=tmp->next)
	{
		Print_Domain(stdout, tmp->ValidityDomain, param_name);
		print_evalue(stdout, &tmp->EP, param_name);	
	}

    //free_evalue_refs(EP);
    //free(EP);
    Free_ParamNames(param_name, nparam);
    Polyhedron_Free(A);
    return 0;
}
