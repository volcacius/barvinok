
#include "ev_operations.h"
#include <polylib/polylibgmp.h>


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
              if (res->x.p->type==polynomial) {
                  /* Add the constant to the constant term of a polynomial*/
                   eadd(e1, &res->x.p->arr[0]);
                   return ;
              }
              else if (res->x.p->type==periodic) {
                          /* Add the constant to all elements of a periodic number */
                          for (i=0; i<res->x.p->size; i++) {
                              eadd(e1, &res->x.p->arr[i]);
                          }
			     
                          return ;
                    } 
                    else {
                            fprintf(stderr, "eadd: cannot add const with vector\n");
                        
                            return;
                    }
     }
     /* add polynomial or periodic to constant 
      * you have to exchange e1 and res, before doing addition */
     
     else if (value_zero_p(e1->d) && value_notzero_p(res->d)) {
              
	      evalue ev;
	      value_init(ev.d);
	      value_assign(ev.d,res->d);
	      value_init(ev.x.n);
	      value_assign(ev.x.n,res->x.n);
	      value_set_si( res->d, 0 );
	      res->x.p=ecopy(e1->x.p);
	          eadd(&ev,res);
	      free_evalue_refs(&ev);
	     
              return ;
     }
     else {   // ((e1->d==0) && (res->d==0)) 
                 if ((e1->x.p->type != res->x.p->type) ) {
		      /* adding to evalues of different type. two cases are possible  
		       * res is periodic and e1 is polynomial, you have to exchange
		       * e1 and res then to add e1 to the constant term of res */
		     if ((res->x.p->type == periodic)&&(e1->x.p->type == polynomial)) {
	               
		          evalue eval;
		          value_set_si( eval.d, 0 );
		          eval.x.p=ecopy(res->x.p);
	                  res->x.p= ecopy(e1->x.p);
                          eadd(&eval,&res->x.p->arr[0]);
			 		         	     
	             }
                     else if ((res->x.p->type == polynomial)&&(e1->x.p->type == periodic)) {
                          /* res is polynomial and e1 is periodic,
		            add e1 to the constant term of res */
			 
			  eadd(e1,&res->x.p->arr[0]);
		     }
	                	 
		     return;
	         }
	         else if (e1->x.p->pos  != res->x.p->pos ) { 
	      	 /* adding evalues of different position (i.e function of different unknowns
		  * to case are possible  */
			   
			 if (res->x.p->type == polynomial) {//  res and e1 are polynomials
			       //  add e1 to the constant term of res
			       
		               eadd(e1,&res->x.p->arr[0]);
		              // value_clear(g); value_clear(m1); value_clear(m2);
		               return;
		          }
		          else {  // res and e1 are pointers to periodic numbers
				  //add e1 to all elements of res 
				   
			          for (i=0;i<res->x.p->size;i++) {
			               eadd(e1,&res->x.p->arr[i]);
			          }
			          return;
		          }
                 
				          
	         }  
                 
                
		 //same type , same pos  and same size
                 if (e1->x.p->size == res->x.p->size) {
	              // add any element in e1 to the corresponding element in res 
	              for (i=0; i<res->x.p->size; i++) {
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
 

/* Computes the product of tow evalues "e1" and "res" and puts the result in "res". you must
 * do a copy of "res" befor calling this function if you nead it after. The vector type of 
 * evalues is not treated here */

void emul (evalue *e1, evalue *res ){
    int i,j;

if((value_zero_p(e1->d)&&e1->x.p->type==evector)||(value_zero_p(res->d)&&(res->x.p->type==evector))) {    
    fprintf(stderr, "emul: do not proced on evector type !\n");
    return;
}
     
   if(value_zero_p(e1->d)&& value_zero_p(res->d)) {
	   if(e1->x.p->type==polynomial && res->x.p->type==polynomial) {
	       if(e1->x.p->pos == res->x.p->pos) {
	       /* Product of tow polynomials of the same variable */
		 
		    evalue *tmp;
	            int size=(((e1->x.p->size)-1) + ((res->x.p->size)-1) +1); 
		    tmp=(evalue *) malloc (sizeof(evalue));
	            value_init(tmp->d);
		    value_set_si(tmp->d,0);
		    tmp->x.p=new_enode(polynomial, size, res->x.p->pos);
		    for(i=0 ;i<e1->x.p->size; i++) {
			    value_init(tmp->x.p->arr[i].d);
			    value_assign(tmp->x.p->arr[i].d, e1->x.p->arr[i].d);
		        if(value_notzero_p(tmp->x.p->arr[i].d)){
				value_init(tmp->x.p->arr[i].x.n);
				value_assign(tmp->x.p->arr[i].x.n, e1->x.p->arr[i].x.n);
			}
			else {
			        tmp->x.p->arr[i].x.p=ecopy (e1->x.p->arr[i].x.p);
					
			}
		    }
		    for(i;i<size;i++){
		         value_init(tmp->x.p->arr[i].d);
			 value_init(tmp->x.p->arr[i].x.n);
			 value_set_si(tmp->x.p->arr[i].d,1);
			 value_set_si(tmp->x.p->arr[i].x.n,0);
		    }
		    for(i=0;i<e1->x.p->size;i++)
			    emul(&(res->x.p->arr[0]),&(tmp->x.p->arr[i]));
		    for(i=1;i<res->x.p->size;i++)    
			    for(j=0;j<e1->x.p->size;j++) {
				    evalue ev;
				    value_init(ev.d);
				    value_assign(ev.d, e1->x.p->arr[j].d);
				    if(value_notzero_p(ev.d)) {
					 value_init(ev.x.n);
					 value_assign(ev.x.n,e1->x.p->arr[j].x.n);
				    }
				    else{
				         ev.x.p=ecopy(e1->x.p->arr[j].x.p);
				    }
				    emul(&(res->x.p->arr[i]),&ev);
			            eadd(&ev,&(tmp->x.p->arr[i+j]));
				    free_evalue_refs(&ev);
			    }
		value_assign(res->d,tmp->d);
		res->x.p=tmp->x.p;
		return;
	       }
	       else {
		  /* Product of tow polynomials of different variables */     
	          
		 for( i=0; i<res->x.p->size ; i++) 
		            emul(e1, &res->x.p->arr[i]);
			  
		 return ;
	       }      
	   }
	   else{
              if(e1->x.p->type==periodic && res->x.p->type==periodic) {   
		if(e1->x.p->pos==res->x.p->pos && e1->x.p->size==res->x.p->size) {
		 /* Product of tow periodics of the same parameter and period */	
                   
		     for(i=0; i<res->x.p->size;i++) 
		           emul(&(e1->x.p->arr[i]), &(res->x.p->arr[i]));
	             		     
		     return;
		}
		else{
                  if(e1->x.p->pos==res->x.p->pos && e1->x.p->size!=res->x.p->size) {  
	           /* Product of tow periodics of the same parameter and different periods */		  
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
	             /* Product of tow periodics of different parameters */
			  
		        for(i=0; i<res->x.p->size; i++)
         	            emul(e1, &(res->x.p->arr[i]));
			
			return;
		  }
		}		       
	      }
	       else {
		  if(e1->x.p->type==periodic && res->x.p->type==polynomial) { 
                  /* Product of a periodic and a polynomial */
			  
		       for(i=0; i<res->x.p->size ; i++)
		            emul(e1, &(res->x.p->arr[i]));    
                 
		       return; 
			       
		  }
		  else {
	             /* Product of a polynomial and a periodic */		  
		      
		       evalue ev;
                       value_init(ev.d);
		       value_set_si(ev.d,0);
		       ev.x.p=ecopy(res->x.p);
		       res->x.p=ecopy(e1->x.p);
		           emul(&ev,res);
		     
		       free_evalue_refs(&ev);	  
		       return ;
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
		     
		 evalue ev;
		 value_init(ev.d);
		 value_init(ev.x.n);
		 value_assign(ev.d, res->d);
		 value_assign(ev.x.n, res->x.n);
		 value_set_si( res->d, 0 );
		 res->x.p=ecopy(e1->x.p);
		       emul(&ev,res);
		
		 free_evalue_refs(&ev);
		 return ;
	     }
	     else {
	       /* Product of a rationel number and an expression (polynomial or peririodic) */ 
	         
		   for (i=0; i<res->x.p->size; i++) 
	              emul(e1, &res->x.p->arr[i]);		 
		  
		 return ;
	     }
       }
   }
   
   return ;
}


