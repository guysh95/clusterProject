#ifndef CHECKP_G_H
#define CHECKP_G_H
#include<stdio.h>

typedef struct newGroup{
    int *members;
    double *F;
    int size;
} group;

void free_group(group *g);
group* allocate_first(int size);
void split(group *g, int *S, group *gPos, group *gNeg);
void writeGroupToFile(group *G, FILE *output);


#endif /*CHECKP_G_H*/
