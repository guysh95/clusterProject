/*
* Created by guy on 05/09/2020.
*/
#include "queue.h"
#include "g.h"
#include <stdlib.h>
#include "checkError.h"

/* allocating an empty queue and returning a pointer to it's location */
queue* allocate_queue(){
    queue *Q = malloc(sizeof(queue));
    checkMalloc(Q);
    Q->head = NULL;
    Q->last = NULL;
    Q->size = 0;
    return Q;
}

/* This function accepts a queue Q and a group G and inserts G the end of the queue */
void insert_queue(queue *Q, group *G){
    item *new = malloc(sizeof(item));
    checkMalloc(new);
    new->data = G;
    new->next = NULL;
    if(Q->head == NULL){
        Q->head = new;
        Q->last = new;
    }
    else{
        Q->last->next = new;
        Q->last = new;
    }
    Q->size++;
}
/*
 * This function pops the first group in Q and returns a pointer to that group.
 * assuming isEmpty(Q) == FALSE
 */
group * pop_queue(queue *Q){
    item *temp;
    group *G = Q->head->data;
    temp = Q->head->next;
    free(Q->head);
    Q->head = temp;
    if(Q->head == NULL)
        Q->last = NULL;
    Q->size--;
    return G;
}

/* this pointer returns TRUE(== 1) if Q doesn't have any nodes left and FALSE(== 0) otherwise. */
int isEmpty(queue *Q){
    if(Q->size == 0)
        return 1;
    return 0;
}
/*
 * This function free's all of Q and it's nodes allocated memory.
 * NOTE: if Q is not empty then because of the project's
 * specific workflow an error will be thrown
 */
void free_queue(queue *Q){
    if(isEmpty(Q)){
        free(Q);
        return;
    }
    printf("queue is not empty but has been sent to free\n");
    exit(1);
}
