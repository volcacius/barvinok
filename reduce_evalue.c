#include <polylib/polylibgmp.h>

/*------------------------------------------------------------*/
/* int eequal (e1, e2)                                        */
/*      e1, e2 : pointers to evalues                          */
/* returns 1 (true) if they are equal, 0 (false) if not       */
/*------------------------------------------------------------*/
static int eequal(evalue *e1,evalue *e2) { 
 
  int i;
  enode *p1, *p2;
  
  if (value_ne(e1->d,e2->d))
    return 0;
  
  /* e1->d == e2->d */
  if (value_notzero_p(e1->d)) {    
    if (value_ne(e1->x.n,e2->x.n))
      return 0;
    
    /* e1->d == e2->d != 0  AND e1->n == e2->n */
    return 1;
  }
  
  /* e1->d == e2->d == 0 */
  p1 = e1->x.p;
  p2 = e2->x.p;
  if (p1->type != p2->type) return 0;
  if (p1->size != p2->size) return 0;
  if (p1->pos  != p2->pos) return 0;
  for (i=0; i<p1->size; i++)
    if (!eequal(&p1->arr[i], &p2->arr[i]) ) 
      return 0;
  return 1;
} /* eequal */

/*------------------------------------------------------------*/
/* void reduce_evalue ( e )                                   */
/*       e : pointer to an evalue                             */
/*------------------------------------------------------------*/
void reduce_evalue (evalue *e) {
  
  enode *p;
  int i, j, k;
  
  if (value_notzero_p(e->d))
    return;	/* a rational number, its already reduced */
  if(!(p = e->x.p))
    return;	/* hum... an overflow probably occured */
  
  /* First reduce the components of p */
  for (i=0; i<p->size; i++)
    reduce_evalue(&p->arr[i]);

  if (p->type==periodic) {
    
    /* Try to reduce the period */
    for (i=1; i<=(p->size)/2; i++) {
      if ((p->size % i)==0) {
	
	/* Can we reduce the size to i ? */
	for (j=0; j<i; j++)
	  for (k=j+i; k<e->x.p->size; k+=i)
	    if (!eequal(&p->arr[j], &p->arr[k])) goto you_lose;

            /* OK, lets do it */
            for (j=i; j<p->size; j++) free_evalue_refs(&p->arr[j]);
            p->size = i;
            break;

you_lose:   /* OK, lets not do it */
            continue;
         }
      }

      /* Try to reduce its strength */
    if (p->size == 1) {
      memcpy(e,&p->arr[0],sizeof(evalue));
      free(p);
    }
  }
  else if (p->type==polynomial) {
	  
    /* Try to reduce the degree */
    for (i=p->size-1;i>=0;i--) {
      if (value_one_p(p->arr[i].d) && value_zero_p(p->arr[i].x.n))
	
	/* Zero coefficient */
	continue;
      else
	break;
    }
    if (i==-1) p->size = 1;
    else if (i+1<p->size) p->size = i+1;

    /* Try to reduce its strength */
    if (p->size == 1) {
      memcpy(e,&p->arr[0],sizeof(evalue));
      free(p);
    }
  }
} /* reduce_evalue */
