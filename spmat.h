#ifndef _SPMAT_H
#define _SPMAT_H

/*
 * spmat is a representation of a sparse matrix.
 * In our project, we use sparse matrix as the adjacency matrix of the graph,
 * in order to achieve fast computing time and low the project's memory use.
 * The structure also contains a few basic functions, all described bellow.
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

    /* compute the change in B[g] modularity that depends on S[index] for mod_max */
    double  (*spmat_modularityAlter)(const struct _spmat *A, int *S, int index);

	/* an implementation of the matrix as and array of linked lists */
	void	*private;
} spmat;

/* Allocates a new array-implemented sparse matrix of size n and nnz values */
spmat* spmat_allocate(int n, int nnz);

/*
 * Assuming preNormRow_array has already been called.
 * this function accepts an index of a member in G and
 * checks if that member is a neighbour of the vertex with the index row.
 * if it is the return value will be 1, otherwise 0.
 */
int normRowHelp_array(const struct _spmat *A, group *G, int row, int index);

/*
 * This function sets a pointer to the current colind values needed
 * by the norm computation of B[g] before normRowHelp_array is being called.
 */
void preNormRow_array(const struct _spmat *A, group *G, int row);

#endif
