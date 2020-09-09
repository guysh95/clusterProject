/*
* Created by litalsol on 16/08/2020. Guy if you're reading this know that you're ugly.
*/

#ifndef CHECKP_G_H
#define CHECKP_G_H
#include<stdio.h>
/*
struct linked_vertexes{
    int vertex;
    int loc_vertex;
    struct linked_vertexes *next;
};
typedef struct linked_vertexes VERTEX;
typedef VERTEX* VER;

typedef struct _group{
    VER vertexes;
    VER last;
    int size;
    void (*free)(struct _group *g);

} group;

group* allocate_g(int size);
group* allocate_sub_g();

group** gSplit_ver2(int *S, group *g);
void writeGroupToFile(group *g, FILE *output);
*/
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
