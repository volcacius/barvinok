#include <assert.h>

#include "ev_operations.h"

void evalue_set_si(evalue *ev, int n, int d) {
    value_set_si(ev->d, d);
    value_init(ev->x.n);
    value_set_si(ev->x.n, n);
}

void evalue_set(evalue *ev, Value n, Value d) {
    value_assign(ev->d, d);
    value_init(ev->x.n);
    value_assign(ev->x.n, n);
}

void aep_evalue(evalue *e, int *ref) {
  
    enode *p;
    int i;
  
    if (value_notzero_p(e->d))
        return;	        /* a rational number, its already reduced */
    if(!(p = e->x.p))
        return;	        /* hum... an overflow probably occured */
  
    /* First check the components of p */
    for (i=0;i<p->size;i++)
        aep_evalue(&p->arr[i],ref);
  
    /* Then p itself */
    p->pos = ref[p->pos-1]+1;
    return;
} /* aep_evalue */

/** Comments */
void addeliminatedparams_evalue(evalue *e,Matrix *CT) {
	
    enode *p;
    int i, j;
    int *ref;

    if (value_notzero_p(e->d))
        return;	         /* a rational number, its already reduced */
    if(!(p = e->x.p))
        return;	         /* hum... an overflow probably occured */
  
    /* Compute ref */
    ref = (int *)malloc(sizeof(int)*(CT->NbRows-1));
    for(i=0;i<CT->NbRows-1;i++)
        for(j=0;j<CT->NbColumns;j++)
            if(value_notzero_p(CT->p[i][j])) {
                ref[i] = j;
                break;
            }
  
    /* Transform the references in e, using ref */
    aep_evalue(e,ref);
    free( ref );
    return;
} /* addeliminatedparams_evalue */

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

you_lose:   	/* OK, lets not do it */
                continue;
            }
        }

        /* Try to reduce its strength */
        if (p->size == 1) {
	    value_clear(e->d);
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
	    value_clear(e->d);
            memcpy(e,&p->arr[0],sizeof(evalue));
            free(p);
        }
    }
    else if (p->type==modulo) {
	  
        /* Try to reduce the degree */
        for (i=p->size-1;i>=1;i--) {
            if (value_one_p(p->arr[i].d) && value_zero_p(p->arr[i].x.n))
	
                /* Zero coefficient */
                continue;
            else
                break;
        }
        if (i==0) p->size = 2;
        else if (i+1<p->size) p->size = i+1;

        /* Try to reduce its strength */
        if (p->size == 2) {
	    value_clear(e->d);
            memcpy(e,&p->arr[1],sizeof(evalue));
            free(p);
        }
    }
} /* reduce_evalue */

void print_evalue(FILE *DST,evalue *e,char **pname) {
  
  if(value_notzero_p(e->d)) {    
    if(value_notone_p(e->d)) {
      value_print(DST,VALUE_FMT,e->x.n);
      fprintf(DST,"/");
      value_print(DST,VALUE_FMT,e->d);
    }  
    else {
      value_print(DST,VALUE_FMT,e->x.n);
    }
  }  
  else
    print_enode(DST,e->x.p,pname);
  return;
} /* print_evalue */

void print_enode(FILE *DST,enode *p,char **pname) {
  
  int i;
  
  if (!p) {
    fprintf(DST, "NULL");
    return;
  }
  if (p->type == evector) {
    fprintf(DST, "{ ");
    for (i=0; i<p->size; i++) {
      print_evalue(DST, &p->arr[i], pname);
      if (i!=(p->size-1))
	fprintf(DST, ", ");
    }
    fprintf(DST, " }\n");
  }
  else if (p->type == polynomial) {
    fprintf(DST, "( ");
    for (i=p->size-1; i>=0; i--) {
      print_evalue(DST, &p->arr[i], pname);
      if (i==1) fprintf(DST, " * %s + ", pname[p->pos-1]);
      else if (i>1) 
	fprintf(DST, " * %s^%d + ", pname[p->pos-1], i);
    }
    fprintf(DST, " )\n");
  }
  else if (p->type == periodic) {
    fprintf(DST, "[ ");
    for (i=0; i<p->size; i++) {
      print_evalue(DST, &p->arr[i], pname);
      if (i!=(p->size-1)) fprintf(DST, ", ");
    }
    fprintf(DST," ]_%s", pname[p->pos-1]);
  }
  else if (p->type == modulo) {
    fprintf(DST, "( ");
    for (i=p->size-1; i>=1; i--) {
      print_evalue(DST, &p->arr[i], pname);
      if (i >= 2) {
        fprintf(DST, " * ");
	if (i > 2)
	  fprintf(DST, "(");
	fprintf(DST, "(");
        print_evalue(DST, &p->arr[0], pname);
	fprintf(DST, ") mod %d", p->pos);
	if (i>2) 
	  fprintf(DST, ")^%d + ", i-1);
	else
	  fprintf(DST, " + ", i-1);
      }
    }
    fprintf(DST, " )\n");
  }
  return;
} /* print_enode */ 

static int mod_term_smaller(evalue *e1, evalue *e2)
{
    if (value_notzero_p(e1->d)) {
	if (value_zero_p(e2->d))
	    return 1;
	return value_lt(e1->x.n, e2->x.n);
    }
    if (value_notzero_p(e2->d))
	return 0;
    if (e1->x.p->pos < e2->x.p->pos)
	return 1;
    else if (e1->x.p->pos > e2->x.p->pos)
	return 0;
    else
	return mod_term_smaller(&e1->x.p->arr[e1->x.p->type==modulo],
				&e2->x.p->arr[e2->x.p->type==modulo]);
}

static void eadd_rev(evalue *e1, evalue *res)
{
    evalue ev;
    value_init(ev.d);
    evalue_copy(&ev, e1);
    eadd(res, &ev);
    free_evalue_refs(res);	  
    *res = ev;
}

static void eadd_rev_cst (evalue *e1, evalue *res)
{
    evalue ev;
    value_init(ev.d);
    evalue_copy(&ev, e1);
    eadd(res, &ev.x.p->arr[ev.x.p->type==modulo]);
    free_evalue_refs(res);	  
    *res = ev;
}

void eadd(evalue *e1,evalue *res) {

 int i; 
    if (value_notzero_p(e1->d) && value_notzero_p(res->d)) {
         /* Add two rational numbers */
	 Value g,m1,m2;
	 value_init(g);
	 value_init(m1);
	 value_init(m2);
	     
         value_multiply(m1,e1->x.n,res->d);
         value_multiply(m2,res->x.n,e1->d);
         value_addto(res->x.n,m1,m2);
         value_multiply(res->d,e1->d,res->d);
         Gcd(res->x.n,res->d,&g);
         if (value_notone_p(g)) {
	      value_division(res->d,res->d,g);
              value_division(res->x.n,res->x.n,g);
         }
         value_clear(g); value_clear(m1); value_clear(m2);
         return ;
     }
     else if (value_notzero_p(e1->d) && value_zero_p(res->d)) {
	  switch (res->x.p->type) {
	  case polynomial:
	      /* Add the constant to the constant term of a polynomial*/
	       eadd(e1, &res->x.p->arr[0]);
	       return ;
	  case periodic:
	      /* Add the constant to all elements of a periodic number */
	      for (i=0; i<res->x.p->size; i++) {
		  eadd(e1, &res->x.p->arr[i]);
	      }
	      return ;
	  case evector:
	      fprintf(stderr, "eadd: cannot add const with vector\n");
	      return;
	  case modulo:
	       eadd(e1, &res->x.p->arr[1]);
	       return ;
	  }
     }
     /* add polynomial or periodic to constant 
      * you have to exchange e1 and res, before doing addition */
     
     else if (value_zero_p(e1->d) && value_notzero_p(res->d)) {
	  eadd_rev(e1, res); 
	  return;
     }
     else {   // ((e1->d==0) && (res->d==0)) 
                 if ((e1->x.p->type != res->x.p->type) ) {
		      /* adding to evalues of different type. two cases are possible  
		       * res is periodic and e1 is polynomial, you have to exchange
		       * e1 and res then to add e1 to the constant term of res */
		     if (e1->x.p->type == polynomial) {
			  eadd_rev_cst(e1, res); 
	             }
                     else if (res->x.p->type == polynomial) {
                          /* res is polynomial and e1 is periodic,
		            add e1 to the constant term of res */
			 
			  eadd(e1,&res->x.p->arr[0]);
		     } else
			assert(0);
	                	 
		     return;
	         }
	         else if (e1->x.p->pos != res->x.p->pos ||
			    (res->x.p->type == modulo &&
			     !eequal(&e1->x.p->arr[0], &res->x.p->arr[0]))) { 
	      	 /* adding evalues of different position (i.e function of different unknowns
		  * to case are possible  */
			   
			switch (res->x.p->type) {
			case modulo:
			    if(mod_term_smaller(res, e1))
				eadd(e1,&res->x.p->arr[1]);
			    else
				eadd_rev_cst(e1, res);
			    return;
			case polynomial: //  res and e1 are polynomials
			       //  add e1 to the constant term of res
			       
			    if(res->x.p->pos < e1->x.p->pos)
		               eadd(e1,&res->x.p->arr[0]);
			    else
				eadd_rev_cst(e1, res);
		              // value_clear(g); value_clear(m1); value_clear(m2);
		               return;
		        case periodic:  // res and e1 are pointers to periodic numbers
				  //add e1 to all elements of res 
				   
			    if(res->x.p->pos < e1->x.p->pos)
			          for (i=0;i<res->x.p->size;i++) {
			               eadd(e1,&res->x.p->arr[i]);
			          }
			    else
				eadd_rev(e1, res);
			    return;
			}
	         }  
                 
                
		 //same type , same pos  and same size
                 if (e1->x.p->size == res->x.p->size) {
	              // add any element in e1 to the corresponding element in res 
		      if (res->x.p->type == modulo)
			assert(eequal(&e1->x.p->arr[0], &res->x.p->arr[0]));
		      i = res->x.p->type == modulo ? 1 : 0;
	              for (; i<res->x.p->size; i++) {
                            eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
                      }
                      return ;
                }
                
		/* Sizes are different */
                if (res->x.p->type==polynomial) {
                    /* VIN100: if e1-size > res-size you have to copy e1 in a   */
                    /* new enode and add res to that new node. If you do not do */
                    /* that, you lose the the upper weight part of e1 !         */

                     if(e1->x.p->size > res->x.p->size) {
                          enode *tmp;
                          tmp = ecopy(e1->x.p);
                          for(i=0;i<res->x.p->size;++i) {
                              eadd(&res->x.p->arr[i], &tmp->arr[i]);
                              //  free_evalue_refs(&res->x.p->arr[i]);
                          }
                          res->x.p = tmp;
    	  
                     }
                     else {
	  	     
                        for (i=0; i<e1->x.p->size ; i++) {
                             eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
                        } 
                      		   
                        return ;
                     } 
                }
                
    /* add two periodics of the same pos (unknown) but whith different sizes (periods) */
                else if (res->x.p->type==periodic) {
		      /* you have to create a new evalue 'ne' in whitch size equals to the lcm
		       of the sizes of e1 and res, then to copy res periodicaly in ne, after
		       to add periodicaly elements of e1 to elements of ne, and finaly to 
		       return ne. */
		       int x,y,p;
		       Value ex, ey ,ep;
		       evalue *ne;
	               value_init(ex); value_init(ey);value_init(ep);
		       x=e1->x.p->size;
	               y= res->x.p->size;
		       value_set_si(ex,e1->x.p->size);
		       value_set_si(ey,res->x.p->size);
		       value_assign (ep,*Lcm(ex,ey));
		       p=(int)mpz_get_si(ep);
	               ne= (evalue *) malloc (sizeof(evalue)); 
	               value_init(ne->d);
	               value_set_si( ne->d,0);
	    	           	    
	               ne->x.p=new_enode(res->x.p->type,p, res->x.p->pos);
	               for(i=0;i<p;i++)  {
		              value_assign(ne->x.p->arr[i].d, res->x.p->arr[i%y].d);
		              if (value_notzero_p(ne->x.p->arr[i].d))   {
				value_init(ne->x.p->arr[i].x.n);
				value_assign(ne->x.p->arr[i].x.n, res->x.p->arr[i%y].x.n);
		              }
		              else { 
			          ne->x.p->arr[i].x.p =ecopy(res->x.p->arr[i%y].x.p);
		              }
	               }
	               for(i=0;i<p;i++)  {
	                    eadd(&e1->x.p->arr[i%x], &ne->x.p->arr[i]);
	               }
      
	              value_assign(res->d, ne->d);
		      res->x.p=ne->x.p;
	    	    
                        return ;
                }
                else { /* evector */
                     fprintf(stderr, "eadd: ?cannot add vectors of different length\n");
                     return ;
                }
     }
     return ;
 }/* eadd  */ 
 
static void emul_rev (evalue *e1, evalue *res)
{
    evalue ev;
    value_init(ev.d);
    evalue_copy(&ev, e1);
    emul(res, &ev);
    free_evalue_refs(res);	  
    *res = ev;
}

static void emul_poly (evalue *e1, evalue *res)
{
    int i, j, o = res->x.p->type == modulo;
    evalue tmp;
    int size=(e1->x.p->size + res->x.p->size - o - 1); 
    value_init(tmp.d);
    value_set_si(tmp.d,0);
    tmp.x.p=new_enode(res->x.p->type, size, res->x.p->pos);
    if (o)
	evalue_copy(&tmp.x.p->arr[0], &e1->x.p->arr[0]);
    for (i=o; i < e1->x.p->size; i++) {
	evalue_copy(&tmp.x.p->arr[i], &e1->x.p->arr[i]);
	emul(&res->x.p->arr[o], &tmp.x.p->arr[i]);
    }
    for (; i<size; i++)
	evalue_set_si(&tmp.x.p->arr[i], 0, 1);
    for (i=o+1; i<res->x.p->size; i++)
	for (j=o; j<e1->x.p->size; j++) {
	    evalue ev;
	    value_init(ev.d);
	    evalue_copy(&ev, &e1->x.p->arr[j]);
	    emul(&res->x.p->arr[i], &ev);
	    eadd(&ev, &tmp.x.p->arr[i+j-o]);
	    free_evalue_refs(&ev);
	}
    free_evalue_refs(res);
    *res = tmp;
}

/* Computes the product of two evalues "e1" and "res" and puts the result in "res". you must
 * do a copy of "res" befor calling this function if you nead it after. The vector type of 
 * evalues is not treated here */

void emul (evalue *e1, evalue *res ){
    int i,j;

if((value_zero_p(e1->d)&&e1->x.p->type==evector)||(value_zero_p(res->d)&&(res->x.p->type==evector))) {    
    fprintf(stderr, "emul: do not proced on evector type !\n");
    return;
}
     
   if(value_zero_p(e1->d)&& value_zero_p(res->d)) {
       switch(e1->x.p->type) {
       case polynomial:
	   switch(res->x.p->type) {
	   case polynomial:
	       if(e1->x.p->pos == res->x.p->pos) {
	       /* Product of two polynomials of the same variable */
		    emul_poly(e1, res);
		    return;
	       }
	       else {
		  /* Product of two polynomials of different variables */     
	          
		if(res->x.p->pos < e1->x.p->pos)
		     for( i=0; i<res->x.p->size ; i++) 
		            emul(e1, &res->x.p->arr[i]);
		else
		    emul_rev(e1, res);
			  
		 return ;
	       }      
	   case periodic:
	   case modulo:
	        /* Product of a polynomial and a periodic or modulo */
		emul_rev(e1, res);
		return;
	   }
       case periodic:
	   switch(res->x.p->type) {
	   case periodic:
		if(e1->x.p->pos==res->x.p->pos && e1->x.p->size==res->x.p->size) {
		 /* Product of two periodics of the same parameter and period */	
                   
		     for(i=0; i<res->x.p->size;i++) 
		           emul(&(e1->x.p->arr[i]), &(res->x.p->arr[i]));
	             		     
		     return;
		}
		else{
                  if(e1->x.p->pos==res->x.p->pos && e1->x.p->size!=res->x.p->size) {  
	           /* Product of two periodics of the same parameter and different periods */		  
		    evalue *newp;
	            Value x,y,z;
		    int ix,iy,lcm;
		    value_init(x); value_init(y);value_init(z);
		    ix=e1->x.p->size;
		    iy=res->x.p->size;
		    value_set_si(x,e1->x.p->size);
		    value_set_si(y,res->x.p->size);
		    value_assign (z,*Lcm(x,y));
		    lcm=(int)mpz_get_si(z);
		    newp= (evalue *) malloc (sizeof(evalue));
		    value_init(newp->d);
		    value_set_si( newp->d,0);
		    newp->x.p=new_enode(periodic,lcm, e1->x.p->pos);
		        for(i=0;i<lcm;i++)  {
		           value_assign(newp->x.p->arr[i].d, res->x.p->arr[i%iy].d);
			   if (value_notzero_p(newp->x.p->arr[i].d))   {
			   value_assign(newp->x.p->arr[i].x.n, res->x.p->arr[i%iy].x.n);
			   }
			   else {
			   newp->x.p->arr[i].x.p =ecopy(res->x.p->arr[i%iy].x.p);	   
			   }

			}
			for(i=0;i<lcm;i++)  
		            emul(&e1->x.p->arr[i%ix], &newp->x.p->arr[i]);		
			  
			value_assign(res->d,newp->d);
			res->x.p=newp->x.p;
			
			  value_clear(x); value_clear(y);value_clear(z); 
			  return ;	 
		  }
		  else {
	             /* Product of two periodics of different parameters */
			  
		        for(i=0; i<res->x.p->size; i++)
         	            emul(e1, &(res->x.p->arr[i]));
			
			return;
		  }
		}		       
	   case polynomial:
                  /* Product of a periodic and a polynomial */
			  
		       for(i=0; i<res->x.p->size ; i++)
		            emul(e1, &(res->x.p->arr[i]));    
                 
		       return; 
			       
	   }		   
       case modulo:
	    switch(res->x.p->type) {
	    case polynomial:
	        for(i=0; i<res->x.p->size ; i++)
		    emul(e1, &(res->x.p->arr[i]));    
	        return; 
	    case periodic:
		assert(0);
	    case modulo:
	        if (e1->x.p->pos == res->x.p->pos &&
			    eequal(&e1->x.p->arr[0], &res->x.p->arr[0]))
		    emul_poly(e1, res);
		else {
		    if(mod_term_smaller(res, e1))
			for(i=1; i<res->x.p->size ; i++)
			    emul(e1, &(res->x.p->arr[i]));    
		    else
			emul_rev(e1, res);
		    return; 
		}
	    }
       }		   
   }
   else {
       if (value_notzero_p(e1->d)&& value_notzero_p(res->d)) {
	   /* Product of two rational numbers */
	
	    Value g;
	    value_init(g);
	    value_multiply(res->d,e1->d,res->d);
	    value_multiply(res->x.n,e1->x.n,res->x.n );
	    Gcd(res->x.n, res->d,&g);
	   if (value_notone_p(g)) {
	       value_division(res->d,res->d,g);
	       value_division(res->x.n,res->x.n,g);
	   }
	   value_clear(g);
	   return ;
       }
       else { 
	     if(value_zero_p(e1->d)&& value_notzero_p(res->d)) { 
	      /* Product of an expression (polynomial or peririodic) and a rational number */
		     
		emul_rev(e1, res);
		 return ;
	     }
	     else {
	       /* Product of a rationel number and an expression (polynomial or peririodic) */ 
	         
		   i = res->x.p->type == modulo ? 1 : 0;
		   for (; i<res->x.p->size; i++) 
	              emul(e1, &res->x.p->arr[i]);		 
		  
		 return ;
	     }
       }
   }
   
   return ;
}

void evalue_copy(evalue *dst, evalue *src)
{
    value_assign(dst->d, src->d);
    if(value_notzero_p(src->d)) {
	 value_init(dst->x.n);
	 value_assign(dst->x.n, src->x.n);
    } else
	 dst->x.p = ecopy(src->x.p);
}

enode *new_enode(enode_type type,int size,int pos) {
  
  enode *res;
  int i;
  
  if(size == 0) {
    fprintf(stderr, "Allocating enode of size 0 !\n" );
    return NULL;
  }
  res = (enode *) malloc(sizeof(enode) + (size-1)*sizeof(evalue));
  res->type = type;
  res->size = size;
  res->pos = pos;
  for(i=0; i<size; i++) {
    value_init(res->arr[i].d);
    value_set_si(res->arr[i].d,0);
    res->arr[i].x.p = 0;
  }
  return res;
} /* new_enode */

enode *ecopy(enode *e) {
  
  enode *res;
  int i;
  
  res = new_enode(e->type,e->size,e->pos);
  for(i=0;i<e->size;++i) {
    value_assign(res->arr[i].d,e->arr[i].d);
    if(value_zero_p(res->arr[i].d))
      res->arr[i].x.p = ecopy(e->arr[i].x.p);
    else {
      value_init(res->arr[i].x.n);
      value_assign(res->arr[i].x.n,e->arr[i].x.n);
    }
  }
  return(res);
} /* ecopy */

int eequal(evalue *e1,evalue *e2) { 
 
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

void free_evalue_refs(evalue *e) {
  
  enode *p;
  int i;
  
  if (value_notzero_p(e->d)) {
    
    /* 'e' stores a constant */
    value_clear(e->d);
    value_clear(e->x.n);
    return; 
  }  
  value_clear(e->d);
  p = e->x.p;
  if (!p) return;	/* null pointer */
  for (i=0; i<p->size; i++) {
    free_evalue_refs(&(p->arr[i]));
  }
  free(p);
  return;
} /* free_evalue_refs */

