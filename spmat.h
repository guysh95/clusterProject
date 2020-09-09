#ifndef _SPMAT_H
#define _SPMAT_H

/*
 * spmat is a representation of a sparse matrix. In our project, we use sparse matrix as the adjacency matrix of the graph.
 * each spmat structure contains it's dimension and an array of linked lists representing only the non-zero cells of the adjacency matrix (the
 * private field). The structure also contains a few basic functions, all described bellow.
 * note that we defined another structure - linked lists. We use this structure in the private field to represent the linked list.
 */
#include "g.h"

typedef struct _spmat {
	/* Matrix size (n*n) */
	int		n;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void	(*add_row)(struct _spmat *A, const int *inputRow, int i, int verK);

	/* Frees all resources used by A */
	void	(*free)(struct _spmat *A);


	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void	(*multSpmat)(const struct _spmat *A, const void *v, double *result, int vType);

	/* returns sum of index row */
    int  (*spmatRowSum)(const struct _spmat *A, group *G, int index);

    double  (*spmat_modularityAlter)(const struct _spmat *A, int *S, int index);

	/* an implementation of the matrix as and array of linked lists */
	void	*private;
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate(int n, int nnz);
int normRowHelp_array(const struct _spmat *A, group *G, int row, int index);
void preNormRow_array(const struct _spmat *A, group *G, int row);

#endif
