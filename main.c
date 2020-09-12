#include<stdio.h>
#include<stdlib.h>
#include "spmat.h"
#include "Bhat.h"
#include "g.h"
#include<time.h>
#include<math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "checkError.h"
#include "queue.h"
#define ZERO 0.0001

void algorithm3(bmat *B, group *G, int *S, queue *P, queue *O);
void writeToOutput(char *filename, queue *O);
int findPartition(bmat* B, group *G, int *S, double *eigenVec, double *randVec);
void sNoDivision(int *S, int dim, int *members, int *endMem);
void reset_eigenVec(double* eigenVec, int dim);
void unmoved_build(int *unmoved, int len);
void mod_max(int *s, bmat *B, group *G, int* unmoved, double *improve, int *indices);
off_t sizeOfFile(const char* filename);


int main(int argc, char* argv[]) {
    FILE *input;
    int dim, *inputRow, i, verK, nnz;
    int M = 0, *verDegrees;
    spmat *A;
    group *G;
    bmat *B;
    queue *P, *O;
    clock_t start = clock(), end;
    if(argc != 3){
        printf("invalid input\n");
        exit(1);
    }

    /**** Reading information from input file into A, verDegrees and M ****/
    input = fopen(argv[1], "r");
    checkFilePtr(input);
    checkFileAccess((int)fread(&dim, sizeof(int), 1, input), 1);
    inputRow = malloc(dim*sizeof(int));
    checkMalloc(inputRow);
    verDegrees = malloc(dim*sizeof(int));
    checkMalloc(verDegrees);
    nnz = sizeOfFile(argv[1])/sizeof(int) - dim - 1;
    A = spmat_allocate(dim, nnz);
    for(i = 0; i < dim; i++) {
        checkFileAccess((int) fread(&verK, sizeof(int), 1, input),1);
        checkFileAccess((int) fread(inputRow, sizeof(int), verK, input), verK);
        M += verK;
        *verDegrees = verK;
        verDegrees++;
        A->add_row(A, inputRow, verK, i);
    }
    if(M == 0){
        printf("the graph doesn't have any edges, this is undefined");
        exit(1);
    }
    verDegrees -= dim;
    fclose(input);
    /**** Allocating main G and B ****/
    G = allocate_first(dim);
    B = malloc(sizeof(bmat));
    checkMalloc(B);
    allocate_bmat(A, verDegrees, M, B);
    /*** Allocating the P and O queues ***/
    P = allocate_queue();
    O = allocate_queue();
    /*** Executing algorithm 3 to find the ideal partition
     * of the network into community structures ***/
    algorithm3(B, G, inputRow, P, O);
    free(inputRow), free_queue(P), free_bmat(B);
    /*** Writing the partition we got into the output file ***/
    writeToOutput(argv[2], O);
    free_queue(O);

    end = clock();
    printf("%d nodes: application done in %f seconds\n", dim, ((double)(end-start))/CLOCKS_PER_SEC);
    return 0;
}
/*
 * This algorithm accepts a modularity matrix B,
 * a first initialized group G, and pre allocated vector and queues S, P and O.
 * then we execute algorithm 3 as in:
 * first we insert G to P, then while P is not empty we pop the first group in P,
 * find it's ideal partition with power iteration, which is executed in findPartition, and mod_max.
 * Then we split the current group and insert the two groups we get to P
 * if they can be partitioned themselves or insert them to O otherwise.
 * At the end O holds the ideal partition into subgroups of the original G.
 */
void algorithm3(bmat *B, group *G, int *S, queue *P, queue *O){
    double *eigenVec, *randVec, *improve;
    int *unmoved, *indices, steps = 0, isEVneg;
    group *gPos, *gNeg;
    eigenVec = malloc(B->dim * sizeof(double));
    checkMalloc(eigenVec);
    randVec = malloc(B->dim * sizeof(double));
    checkMalloc(randVec);
    improve = eigenVec;
    unmoved = (int*)randVec;
    indices = unmoved + B->dim;
    insert_queue(P, G);
    while(!isEmpty(P)){
        steps++;
        G = pop_queue(P);
        isEVneg = findPartition(B, G, S, eigenVec, randVec); /* algorithm 2 */
        if(!isEVneg)
            /*
             * isEVneg == TRUE if the leading eigen value of B[g] is negative
             * or if the modularity given by B[g] and S is negative.
             */
            mod_max(S, B, G, unmoved, improve, indices); /* algorithm 4*/
        gPos = malloc(sizeof(group));
        checkMalloc(gPos);
        gNeg = malloc(sizeof(group));
        checkMalloc(gNeg);
        split(G, S, gPos, gNeg);
        if(gPos->size == 0 || gNeg->size == 0){
            if(gPos->size == 0){
                insert_queue(O, gNeg);
                free_group(gPos);
            }
            else{
                insert_queue(O, gPos);
                free_group(gNeg);
            }
        }
        else{
            if(gPos->size == 1)
                insert_queue(O, gPos);
            else
                insert_queue(P, gPos);
            if(gNeg->size == 1)
                insert_queue(O, gNeg);
            else
                insert_queue(P, gNeg);
        }
    }
    free(eigenVec), free(randVec);
}

/*
 * This function accepts an absolute path to a file and a queue O,
 * and then writes all of the information about the groups in O to that file.
 */
void writeToOutput(char *filename, queue *O){
    FILE *output;
    group *G;
    output = fopen(filename, "w");
    checkFilePtr(output);
    checkFileAccess((int)fwrite(&O->size, sizeof(int), 1, output), 1);
    while(!isEmpty(O)){
        G = pop_queue(O);
        writeGroupToFile(G, output);
        free_group(G);
    }
    fclose(output);
}

/*
 * algorithm 2
 * accepts bmat B and compute its modularity
 * returns optimal division of g, represented by +1\-1\0 vector S
 * 0 values represent indices of vertices that are not in G.
 * the function return's 1 is the eigenVal computed by the power iteration in negative,
 * or if the modularity of B[g] with the partiton S is negative.
 */
int findPartition(bmat* B, group *G, int *S, double *eigenVec, double *randVec){
    double *temp, *goal, eigenVal;
    int *members = G->members, *endMem = members + G->size;
    int dim, i;

    /* getting leading eigen val and vec of B */
    reset_eigenVec(eigenVec, B->dim);
    findEigenPair(B, G, eigenVec, &eigenVal, randVec); /* uses power iteration */
    dim = B->dim;
    if(eigenVal <= ZERO){
        /* g is indivisible */
        sNoDivision(S, dim, members, endMem);
        return 1;
    }
    /* calculating s */
    goal = eigenVec + dim; /* goal is pointer instead of using index */
    for(temp = eigenVec, i = 0; temp < goal; temp++, i++, S++){ /* calculate s with eigen value */
        if(members == endMem)
            *S = 0;
        else{
            if (i != *members)
                *S = 0;
            else{
                *S = *temp > 0 ? 1 : -1;
                members++;
            }
        }
    }
    S -= dim;
    if(computeModularity(B, G, S, randVec) <= ZERO){
        /* g is indivisible */
        sNoDivision(S, dim, members, endMem);
        return 1;
    }
    return 0;
}
/*
 * This function reset eigenVec values to 0.0 before the start of the current computation
 * this is used in findPartition
 */
void reset_eigenVec(double* eigenVec, int dim){
    double* ptr;
    for(ptr = eigenVec ; ptr < eigenVec + dim; ptr++)
        *ptr = 0.0;
}

/*
 * This function set all of S values that are relevant to G to 1.
 * The rest of S values are set to 0.
 */
void sNoDivision(int *S, int dim, int *members, int *endMem){
    int i;
    for(i = 0; i < dim; i++, S++) {
        if (members == endMem)
            *S = 0;
        else {
            if (i != *members)
                *S = 0;
            else{
                *S = 1;
                members++;
            }
        }
    }
}

/*
 * The function that optimizes an initial partition of Bhat.
 * input: the initial partition s of B, as calculated with power iteration.
 * in the end od the process s is updated to be the optimized partition of G.
 */
void mod_max(int *s, bmat *B, group *G, int* unmoved, double *improve, int *indices){
    int len;
    int *head;
    int i, j, k, max_ind, max_imp_ind, count = 0, limit = 3000 * B->dim;
    double del_q, score, max, max_imp;
    len = G->size;
    del_q = 1;
    while (del_q > ZERO) {
        if(count > limit){
            printf("aborted due to suspicion of infinite loop\n");
            exit(1);
        }
        count++;
        unmoved_build(unmoved, len); /* unmoved is a boolean array in the size of G */
        max_imp = -HUGE_VAL;
        for (k = 0; k < len; k++){  /*line 3 in the pseudo code*/
            head = unmoved;
            max = -HUGE_VAL;
            for (i = 0; i < len; i++) {
                if(*head == 0){
                    head++;
                    continue;
                }
                s[G->members[i]] = -1 * s[G->members[i]];
                /* we compute the alteration to the modularity
                 * directly instead of computing the new modularity */
                score = modularityAlteration(B, G, s, i);
                s[G->members[i]] = -1 * s[G->members[i]];
                if (score > max) {
                    max_ind = i;
                    max = score;
                }
                head++;
            }
            s[G->members[max_ind]] = -1 * s[G->members[max_ind]];
            indices[k] = max_ind;
            improve[k] = (k == 0) ? max : (improve[k - 1] + max);
            unmoved[max_ind] = 0; /* moving vertex in max_ind */
            if (improve[k] > max_imp) {
                max_imp = improve[k];
                max_imp_ind = k;
            }
        }
        /* making sure that if the max improvment to the modularity is 0
         * then the G is not divided */
        if(max_imp == 0)
            max_imp_ind = len-1;

        for (i = len - 1; i >= max_imp_ind + 1; i--) {
            j = indices[i];
            s[G->members[j]] = -1 * s[G->members[j]];
        }
        del_q = (max_imp_ind == len - 1) ? 0 : improve[max_imp_ind];
    }
}
/*
 * This function set all of the pre-allocated unmoved vector values to 1.
 */
void unmoved_build(int *unmoved, int len){
    int *ptr, *goal = unmoved + len;
    for(ptr = unmoved; ptr < goal; ptr++)
        *ptr = 1;
}

/*
 * This function accepts an absolute path to a file and returns the file's size in bytes,
 * using sys/stat.h library.
 */
off_t sizeOfFile(const char* filename){
    struct stat st;
    if(stat(filename, &st) == 0)
        return st.st_size;
    printf("can't compute file size\n");
    exit(1);
}
