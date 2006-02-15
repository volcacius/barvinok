/* 
 *	Bernstein Expansion
 *
 *	- polylib functions
 */


#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include <polylib/polylibgmp.h>

#include "bernstein.h"

// global to be used on the coefficient max
extern Polyhedron *VD = NULL;	// current validity domain

/* Converts a value (p) to a long long integer */
long long value2longlong(Value p)
{
	char *str; 
	long long l;

	str = mpz_get_str(0,10, p);
	l = atoll(str);
	free(str);

	return l;
}


/* Given a domain, converts the vertices matrix to an long long matrix and then
   sends it to ginac functions */
void doExpansion(Param_Polyhedron *PP, Param_Domain *Q, unsigned int nb_param
		 , char **param_name)
{
	Param_Vertices *V;
	unsigned bx;
	unsigned int nbVertices = 0, nbColumns, nbRows;
	int i,ix;

	// counts the number of vertices and dimensions
	for(i=0,ix=0,bx=MSB,V=PP->V; V && i<PP->nbV; i++,V=V->next) {
		if (Q->F[ix] & bx) {
			nbVertices++;
			nbColumns = V->Vertex->NbColumns;
			nbRows = V->Vertex->NbRows;
		}
		NEXT(ix,bx);
	}
 	long long *matrix = calloc(sizeof(long long), nbVertices * nbRows * nbColumns);

	// convert matrix to long long
	convertMatrix(matrix, PP, Q, param_name, nbColumns, nbRows);

	// convert vertices and send to ginac
	verticesConvert(matrix, nbRows, nbColumns, nbVertices, param_name, nb_param);

	free(matrix);
}


/*
 * Convert the matrix pointed by p to long long format (store it in matrix)
 */
void convertMatrix(long long *matrix, Param_Polyhedron *PP, Param_Domain *Q, char **param_name
		   , unsigned int nbColumns, unsigned int nbRows)
{
	Param_Vertices *V;
	Value *p;
	unsigned int vertexNb = 0;
	unsigned bx;
	unsigned int j,k;
	int i,ix;

#ifdef DEBUG
	fprintf(stdout,"Vertices:\n");
#endif
	for(i=0,ix=0,bx=MSB,V=PP->V; V && i<PP->nbV; i++,V=V->next) {
		int nr,nc;

		if (Q->F[ix] & bx) {
#ifdef DEBUG
			Print_Vertex(stdout,V->Vertex,param_name);
			fprintf(stdout,"\n");
#endif
			p=*(V->Vertex->p);
			nr = V->Vertex->NbRows;
			nc = V->Vertex->NbColumns;

			// copy matrix to long long
			for(j = 0; j < V->Vertex->NbColumns; j++) {
				for(k = 0; k < V->Vertex->NbRows; k++) {
					// TODO: loosing precision to long long
					long long l = value2longlong(*p++);
					 matrix[vertexNb*nbRows*nbColumns+nr*j+k] = l;
				}
			}
			vertexNb++;
		}
		NEXT(ix,bx);
	}
}


/* Prints the domain corresponding to Q */
void printDomain(Param_Domain *Q, Polyhedron *B, char **param_name)
{

	printf("\nDomain: \n");

	if(VD != NULL) {
		Domain_Free(VD);
	}
	VD = DomainSimplify(Q->Domain, B, MAXRAYS);
	Print_Domain(stdout, VD, param_name);
}


/* Reads the polynomial matrix, converts it to long long precision and calls ginac functions */
Matrix *readPolynomial(unsigned int nbVariables, unsigned int nbParams, char **param_name)
{
	Matrix *polynomial;

	polynomial = Matrix_Read();

#ifdef DEBUG
	/* Print the polynomial matrix */
	printf("================================\n");
	printf("Polynomial: \n");
	Matrix_Print(stdout, P_VALUE_FMT, polynomial);
	printf("================================\n");
#endif

	long long *matrix = matrix2longlong(polynomial);
	// parameters in the polynomial coefficients
	if(polynomial->NbColumns == nbVariables) {
		unsigned int nbCoefficients = polynomial->NbRows;;
		unsigned int i;

		// FIXME: free
		long long **llPolynomialCoefficients = (long long **) calloc(sizeof(long *), nbCoefficients);
		Matrix **mPolynomialCoefficients = (Matrix **) calloc(sizeof(Matrix *), nbCoefficients);
		unsigned int *llRows = (unsigned int *) calloc(sizeof(unsigned int), nbCoefficients);
		unsigned int *llColumns = (unsigned int *) calloc(sizeof(unsigned int), nbCoefficients);;

		for(i = 0; i < nbCoefficients; i++) {
			// read the matrix and set rows and columns number.
			mPolynomialCoefficients[i] = Matrix_Read();
			llRows[i] = mPolynomialCoefficients[i]->NbRows;
			llColumns[i] = mPolynomialCoefficients[i]->NbColumns;
#ifdef DEBUG
			/* Print the i coefficient matrix */
			printf("================================\n");
			printf("Coefficient i: \n");
			Matrix_Print(stdout, P_VALUE_FMT, mPolynomialCoefficients[i]);
			printf("================================\n");
#endif
			llPolynomialCoefficients[i] = matrix2longlong(mPolynomialCoefficients[i]);
			Matrix_Free(mPolynomialCoefficients[i]);
		}
		free(mPolynomialCoefficients);
		polyConvertParameters(matrix, polynomial->NbRows, polynomial->NbColumns, llPolynomialCoefficients
				      , llRows, llColumns, nbParams, param_name);
		for(i = 0; i < nbCoefficients; i++) {
			free(llPolynomialCoefficients[i]);
		}
		free(llRows);
		free(llColumns);
		free(llPolynomialCoefficients);

	} else {
		polyConvert(matrix, polynomial->NbRows, polynomial->NbColumns, nbParams, param_name);
	}
	free(matrix);

	return polynomial;
}


/* Converts a Matrix from the polylib format to long long* format */
long long *matrix2longlong(Matrix *M)
{
 	long long *matrix = calloc(sizeof(long long), M->NbRows * M->NbColumns);
	unsigned int nr,nc;
	unsigned int j,k;
	Value *p;

	p=*(M->p);
	nr = M->NbRows;
	nc = M->NbColumns;

	/* Copy matrix to long long */
	for(j = 0; j < nc; j++) {
		for(k = 0; k < nr; k++) {
			// TODO: loosing precision to long long
			long long l = value2longlong(*p++);
			matrix[nr*j+k] = l;
		}
	}
	return matrix;
}





/* Sets a value on the Matrix A, A(r, c) = v */
void Matrix_Set(int r, int c, int v, Matrix * A)
{
        Value val;
        value_init(val);
        value_set_si(val,v);
        value_assign(A->p[r][c] , val);
}

/* Converts a *longlong matrix to polylib format */
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns)
{
	Matrix *retval;
	unsigned int i, j;

	retval = Matrix_Alloc(rows, columns);
	for(i = 0; i < rows; i++) {
		for(j = 0; j < columns; j++) {
			Matrix_Set(i, j, (int) M[i*rows+j], retval);
		}
	}
	return retval;
}


/* Check  if a given set of constraints (M matrix) holds or not in the validity domain */
unsigned  checkConstraint(long long *M, unsigned int rows, unsigned int columns)
{
	Matrix *maxConstraints;
	Polyhedron *mC, *newB;

	maxConstraints = longlong2polylib(M, rows, columns);

#ifdef DEBUG
	printf("Max/Min Constraints (polylib Format): \n");
	Matrix_Print(stdout, P_VALUE_FMT, maxConstraints);
	printf("\n");

	printf("Original Validity Domain: \n");
	Print_Domain(stdout, VD, param_name);
	printf("\n");
#endif

	mC =Constraints2Polyhedron(maxConstraints, MAXRAYS);
	newB = DomainIntersection(VD, mC, MAXRAYS);

#ifdef DEBUG
	printf("Intersection betwenn Max/Min/Pos/Neg Constraints and \nOriginal Validity Domain: \n");
	Print_Domain(stdout, newB, param_name);
#endif

	if(emptyQ(newB)) {
#ifdef DEBUG
		printf("Proposed Coefficient is not the max/min/pos/neg.\n");
#endif
		return 0;
	} else {
#ifdef DEBUG
		printf("Proposed Coefficient is the max/min/pos/neg.\n");
#endif
		return 1;
	}

	/* free */
	Matrix_Free(maxConstraints);
	Domain_Free(mC);
	Domain_Free(newB);

	return 0;
}
