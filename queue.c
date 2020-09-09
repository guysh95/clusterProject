/*
* Created by guy on 05/09/2020.
*/
#include "queue.h"
#include "g.h"
#include <stdlib.h>
#include "checkError.h"

queue* allocate_queue(){
    queue *Q = malloc(sizeof(queue));
    checkMalloc(Q);
    Q->head = NULL;
    Q->last = NULL;
    Q->size = 0;
    return Q;
}

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
 * output is the oldest bmat in Queue
 * assuming isEmpty(Q) == false
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

int isEmpty(queue *Q){
    if(Q->size == 0)
        return 1;
    return 0;
}
/* assuming queue is empty, error if not empty */
void free_queue(queue *Q){
    if(isEmpty(Q)){
        free(Q);
        return;
    }
    printf("queue is not empty but has been sent to free\n");
    exit(1);
}
