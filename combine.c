#include "combine.h"

/*                                                                                 
 * Vector p3 is a linear combination of two vectors (p1 and p2) such that
 * p3[pos] is zero. First element of each vector (p1,p2,p3) is a status
 * element and is not changed in p3. The value of 'pos' may be 0 however.
 * The parameter 'length' does not include status element one.
 * 
 * This function was copied from the PolyLib source.
 */                                                                                
void Combine(Value *p1, Value *p2, Value *p3, int pos, unsigned length) {  

  Value a1, a2, gcd;
  Value abs_a1,abs_a2,neg_a1;
                                                                                   
  /* Initialize all the 'Value' variables */                                         
  value_init(a1); value_init(a2); value_init(gcd);
  value_init(abs_a1); value_init(abs_a2); value_init(neg_a1);
                                                                                      
  /* a1 = p1[pos] */
  value_assign(a1,p1[pos]);                                                         

  /* a2 = p2[pos] */
  value_assign(a2,p2[pos]);                                                        

  /* a1_abs = |a1| */
  value_absolute(abs_a1,a1);                                                       

  /* a2_abs = |a2| */
  value_absolute(abs_a2,a2);

  /* gcd  = Gcd(abs(a1), abs(a2)) */
  Gcd(abs_a1,abs_a2,&gcd);

  /* a1 = a1/gcd */
  value_division (a1,a1,gcd);
                                                                                   
  /* a2 = a2/gcd */
  value_division (a2,a2,gcd);

  /* neg_a1 = -(a1) */
  value_oppose(neg_a1,a1);

  Vector_Combine(p1+1,p2+1,p3+1,a2,neg_a1,length);
  Vector_Normalize(p3+1,length);
                                                                                      
  /* Clear all the 'Value' variables */                                            
  value_clear(a1); value_clear(a2); value_clear(gcd);
  value_clear(abs_a1); value_clear(abs_a2); value_clear(neg_a1);
                                                                                    
  return;                                                                          
} /* Combine */                                                                      

