#ifndef EV_OPERATIONS_H
#define EV_OPERATIONS_H

#include <polylib/polylibgmp.h>

#define polynomial new_polynomial
#define periodic new_periodic
#define evector new_evector
#define enode_type new_enode_type
#define enode_node new_enode_node
#define _enumeration _new_enumeration
#define Enumeration _new_Enumeration
#define enode _new_enode
#define _enode _new__enode
#define evalue _new_evalue
#define _evalue _new__evalue
#define eadd _new_eadd
#define ecopy _new_ecopy
#define new_enode _new_new_enode
#define free_evalue_refs _new_free_evalue_refs
#define print_evalue _new_print_evalue
#define print_enode _new_print_enode
#define reduce_evalue _new_reduce

typedef enum { polynomial, periodic, evector, modulo } enode_type;

typedef struct _evalue {
  Value d;              /* denominator */
  union {
    Value n;            /* numerator (if denominator != 0) */
    struct _enode *p;	/* pointer   (if denominator == 0) */
  } x;
} evalue;

typedef struct _enode {
  enode_type type;      /* polynomial or periodic or evector */
  int size;             /* number of attached pointers */
  int pos;	        /* parameter position */
  evalue arr[1];        /* array of rational/pointer */
} enode;

typedef struct _enumeration {
  
  Polyhedron *ValidityDomain;    /* contraints on the parameters     */
  evalue EP;                     /* dimension = combined space       */
  struct _enumeration *next;     /* Ehrhart Polynomial, corresponding
	                            to parameter values inside the
                                    domain ValidityDomain below      */
} Enumeration;


void evalue_set_si(evalue *ev, int n, int d);
void evalue_set(evalue *ev, Value n, Value d);
void evalue_copy(evalue *dst, evalue *src);
enode *new_enode(enode_type type,int size,int pos);
enode *ecopy(enode *e);
int eequal(evalue *e1,evalue *e2);
void free_evalue_refs(evalue *e);
void print_evalue(FILE *DST,evalue *e,char **pname);
void print_enode(FILE *DST,enode *p,char **pname);
void reduce_evalue (evalue *e);
void aep_evalue(evalue *e, int *ref);
void addeliminatedparams_evalue(evalue *e,Matrix *CT);
void eadd(evalue *e1,evalue *res);
void emul (evalue *e1, evalue *res );
#endif
