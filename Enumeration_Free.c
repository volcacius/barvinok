#include <polylib/polylibgmp.h>

void Enumeration_Free(Enumeration *en)
{
  Enumeration *ee;

  while( en )
  {
	  free_evalue_refs( &(en->EP) );
	  Polyhedron_Free( en->ValidityDomain );
	  ee = en ->next;
	  free( en );
	  en = ee;
  }
}
