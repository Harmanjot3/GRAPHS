#include <stdio.h>
#include <stdlib.h>
#include "graphs.h"
#include <limits.h>
#define NO_PATH_FOUND 103
#define MEMORY_ALLOCATION_FAILED 110

typedef struct StackNode {
    int data;
    struct StackNode* next;
} StackNode;

typedef struct MinHeapNode {
    int vertex;
    int distance;
} MinHeapNode;

typedef struct MinHeap {
    int size;
    int capacity;
    int *pos;
    MinHeapNode **array;
} MinHeap;
////void *myMalloc(size_t size) {
    ////return malloc(size);
////}
void decreaseKey(MinHeap* minHeap, int index, int distance) {
    minHeap->array[index]->distance = distance;
    while (index != 0 && minHeap->array[index]->distance < minHeap->array[(index - 1) / 2]->distance) {
        MinHeapNode *temp = minHeap->array[index];
        minHeap->array[index] = minHeap->array[(index - 1) / 2];
        minHeap->array[(index - 1) / 2] = temp;
        minHeap->pos[temp->vertex] = (index - 1) / 2;
        minHeap->pos[minHeap->array[index]->vertex] = index;
        index = (index - 1) / 2;
    }
}

void minHeapify(MinHeap* minHeap, int index) {
    int smallest = index;
    int left = 2 * index + 1;
    int right = 2 * index + 2;

    if (left < minHeap->size && minHeap->array[left]->distance < minHeap->array[smallest]->distance) {
        smallest = left;
    }
    if (right < minHeap->size && minHeap->array[right]->distance < minHeap->array[smallest]->distance) {
        smallest = right;
    }
    if (smallest != index) {
        MinHeapNode *temp = minHeap->array[index];
        minHeap->array[index] = minHeap->array[smallest];
        minHeap->array[smallest] = temp;
        minHeap->pos[temp->vertex] = smallest;
        minHeap->pos[minHeap->array[index]->vertex] = index;
        minHeapify(minHeap, smallest);
    }
}

void push(StackNode** stack, int data) {
    StackNode* newNode = (StackNode*)malloc(sizeof(StackNode));
    newNode->data = data;
    newNode->next = *stack;
    *stack = newNode;
}

int pop(StackNode** stack) {
    if (*stack == NULL) {
        return -1;
    }
    StackNode* temp = *stack;
    int data = temp->data;
    *stack = temp->next;
    free(temp);
    return data;
}

int isEmpty(StackNode* stack) {
    return stack == NULL;
}

void initializeDijkstraTable(DijkstraTable *pTable, int startNode) {
    for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
        pTable->table[i].distance = INT_MAX;
        pTable->table[i].predecessor = -1;
    }
    pTable->table[startNode].distance = 0;
}

int allNodesVisited(bool visited[]) {
    for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
        if (!visited[i]) {
            return 0;
        }
    }
    return 1;
}

int findMinDistanceVertex(DijkstraTable *pTable, bool visited[]) {
    int minDistance = INT_MAX;
    int minVertex = -1;
    for (int v = 0; v < NUMBER_OF_VERTICES; v++) {
        if (!visited[v] && pTable->table[v].distance < minDistance) {
            minDistance = pTable->table[v].distance;
            minVertex = v;
        }
    }
    return minVertex;
}


void reverseArray(int arr[], int size) {
    int start = 0;
    int end = size - 1;
    while (start < end) {
        int temp = arr[start];
        arr[start] = arr[end];
        arr[end] = temp;
        start++;
        end--;
    }
}

MinHeap* createMinHeap(int capacity) {
    MinHeap* minHeap = (MinHeap*)malloc(sizeof(MinHeap));
    if (minHeap == NULL) {
        // Handle memory allocation error
        return NULL;
    }
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->pos = (int*)malloc(capacity * sizeof(int));
    minHeap->array = (MinHeapNode**)malloc(capacity * sizeof(MinHeapNode*));
    if (minHeap->pos == NULL || minHeap->array == NULL) {
        // Handle memory allocation error
        free(minHeap->pos);
        free(minHeap->array);
        free(minHeap);
        return NULL;
    }
    return minHeap;
}

void insertMinHeap(MinHeap* minHeap, MinHeapNode minHeapNode) {
    if (minHeap->size >= minHeap->capacity) {
        // Handle heap overflow
        return;
    }
    int i = minHeap->size;
    minHeap->array[i] = (MinHeapNode*)malloc(sizeof(MinHeapNode));
    if (minHeap->array[i] == NULL) {
        // Handle memory allocation error
        return;
    }
    minHeap->array[i]->vertex = minHeapNode.vertex;
    minHeap->array[i]->distance = minHeapNode.distance;
    minHeap->pos[minHeapNode.vertex] = i;
    minHeap->size++;
    decreaseKey(minHeap, minHeap->size - 1, minHeapNode.distance);
}

int isEmptyMinHeap(MinHeap* minHeap) {
    return minHeap->size == 0;
}

MinHeapNode extractMin(MinHeap* minHeap) {
    if (isEmptyMinHeap(minHeap)) {
        // Handle heap underflow
        MinHeapNode emptyNode = {-1, -1}; // Or any appropriate default node
        return emptyNode;
    }
    MinHeapNode* minHeapNode = minHeap->array[0];
    MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;
    minHeap->pos[minHeap->array[0]->vertex] = 0;
    minHeap->pos[minHeapNode->vertex] = minHeap->size - 1;
    minHeap->size--;
    minHeapify(minHeap, 0);
    return *minHeapNode;
}

void freeMinHeap(MinHeap* minHeap) {
    for (int i = 0; i < minHeap->size; i++) {
        free(minHeap->array[i]);
    }
    free(minHeap->array);
    free(minHeap->pos);
    free(minHeap);
}



AdjacencyMatrix* createAdjacencyMatrix(int defaultEdgeValue)
{
    AdjacencyMatrix* matrix = (AdjacencyMatrix*)myMalloc(sizeof(AdjacencyMatrix));
    if (matrix == NULL) {
        return NULL; // Memory allocation error
    }

    for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
        for (int j = 0; j < NUMBER_OF_VERTICES; j++) {
            matrix->matrix[i][j] = defaultEdgeValue;
        }
    }

    return matrix;
}

int addEdge(AdjacencyMatrix *pMatrix, int src, int dest, int weight)
{
    if (pMatrix == NULL || src < 0 || dest < 0 || src >= NUMBER_OF_VERTICES || dest >= NUMBER_OF_VERTICES) {
        return INVALID_INPUT_PARAMETER;
    }

    pMatrix->matrix[src][dest] = weight;
    return SUCCESS;
}

int addEdges(AdjacencyMatrix *pMatrix, Edge edges[], int edgeNum)
{
    if (pMatrix == NULL || edges == NULL || edgeNum <= 0) {
        return INVALID_INPUT_PARAMETER;
    }

    int successCount = 0;
    for (int i = 0; i < edgeNum; i++) {
        int result = addEdge(pMatrix, edges[i].src, edges[i].dest, edges[i].weight);
        if (result == SUCCESS) {
            successCount++;
        }
    }

    if (successCount == edgeNum) {
        return SUCCESS;
    } else if (successCount == 0) {
        return MEMORY_ALLOCATION_ERROR;
    } else {
        return PARTIAL_SUCCESS;
    }
}

int doDepthFirstTraversal(AdjacencyMatrix *pMatrix, int startingNode, int traversalOutput[])
{
    if (pMatrix == NULL || startingNode < 0 || startingNode >= NUMBER_OF_VERTICES || traversalOutput == NULL) {
        return INVALID_INPUT_PARAMETER;
    }

    // Step 2: Create visited array
    bool visited[NUMBER_OF_VERTICES];
    for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
        visited[i] = false;
    }

    // Step 3: Create stack and mark starting node as visited
    StackNode *stack = NULL;
    push(&stack, startingNode);
    visited[startingNode] = true;

    int currentIndex = 0;

    // Step 4: Depth-first traversal
    while (!isEmpty(stack)) {
        int currentNode = pop(&stack);
        traversalOutput[currentIndex++] = currentNode;

        // Explore adjacent nodes
        for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
            if (pMatrix->matrix[currentNode][i] != 0 && !visited[i]) {
                push(&stack, i);
                visited[i] = true;
            }
        }
    }

    // Step 5: Free memory if needed

    return SUCCESS;
}

int loadMatrixFromFile(AdjacencyMatrix *pMatrix, char filename[])
{
    if (pMatrix == NULL || filename == NULL) {
        return INVALID_INPUT_PARAMETER;
    }

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        return FILE_IO_ERROR;
    }

    for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
        for (int j = 0; j < NUMBER_OF_VERTICES; j++) {
            if (fscanf(file, "%d", &(pMatrix->matrix[i][j])) != 1) {
                fclose(file);
                return FILE_IO_ERROR;
            }
        }
    }

    fclose(file);
    return SUCCESS;
}

int doDijsktraAlgorithm(AdjacencyMatrix *pMatrix, DijkstraTable *pTable, int startNode) {
    if (pMatrix == NULL || pTable == NULL || startNode < 0 || startNode >= NUMBER_OF_VERTICES) {
        return INVALID_INPUT_PARAMETER;
    }
    initializeDijkstraTable(pTable, startNode);
    bool visited[NUMBER_OF_VERTICES];
    for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
        visited[i] = false;
    }
    while (!allNodesVisited(visited)) {
        int minVertex = findMinDistanceVertex(pTable, visited);
        visited[minVertex] = true;
        for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
            if (!visited[i] && pMatrix->matrix[minVertex][i] != 0) {
                int newDistance = pTable->table[minVertex].distance + pMatrix->matrix[minVertex][i];
                if (newDistance < pTable->table[i].distance) {
                    pTable->table[i].distance = newDistance;
                    pTable->table[i].predecessor = minVertex;
                }
            }
        }
    }
    return SUCCESS;
}

int findShortestPathTo(DijkstraTable *pTable, int nodeFrom, int nodeTo, int pathFound[]) {
    if (pTable == NULL || nodeFrom < 0 || nodeFrom >= NUMBER_OF_VERTICES || nodeTo < 0 || nodeTo >= NUMBER_OF_VERTICES || pathFound == NULL) {
        return INVALID_INPUT_PARAMETER;
    }

    // Check if a path exists from nodeFrom to nodeTo
    if (pTable->table[nodeTo].predecessor == -1) {
        return NO_PATH_FOUND; // No path exists
    }

    // Construct the shortest path by traversing the predecessor chain
    int currentNode = nodeTo;
    int pathIndex = 0;

    while (currentNode != nodeFrom) {
        pathFound[pathIndex++] = currentNode;
        currentNode = pTable->table[currentNode].predecessor;
    }

    // Add the starting node to the path
    pathFound[pathIndex++] = nodeFrom;

    // Reverse the pathFound array to get the correct order
    reverseArray(pathFound, pathIndex);

    return SUCCESS;
}

int addEdgeToAdjacencyList(AdjacencyList *pList, int src, int dest, int weight)
{
    if (pList == NULL || src < 0 || dest < 0 || src >= NUMBER_OF_VERTICES || dest >= NUMBER_OF_VERTICES) {
        return INVALID_INPUT_PARAMETER;
    }

    ListNode *node = (ListNode*)myMalloc(sizeof(ListNode));
    if (node == NULL) {
        return MEMORY_ALLOCATION_ERROR;
    }

    node->destNode = dest;
    node->weight = weight;
    node->next = pList->adjacencyList[src];
    pList->adjacencyList[src] = node;

    return SUCCESS;
}

int populateAdjacencyList(AdjacencyList *pList, AdjacencyMatrix *pMatrix)
{
    if (pList == NULL || pMatrix == NULL) {
        return INVALID_INPUT_PARAMETER;
    }

    for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
        for (int j = 0; j < NUMBER_OF_VERTICES; j++) {
            if (pMatrix->matrix[i][j] != 0) {
                int result = addEdgeToAdjacencyList(pList, i, j, pMatrix->matrix[i][j]);
                if (result != SUCCESS) {
                    return result;
                }
            }
        }
    }

    return SUCCESS;
}



int doDijsktraAlgorithmOnAdjacencyList(AdjacencyList *pList, DijkstraTable *pTable, int startNode) {
    if (pList == NULL || pTable == NULL || startNode < 0 || startNode >= NUMBER_OF_VERTICES) {
        return INVALID_INPUT_PARAMETER;
    }
    // Initialize Dijkstra table
    for (int i = 0; i < NUMBER_OF_VERTICES; i++) {
        pTable->table[i].distance = INT_MAX;
        pTable->table[i].predecessor = -1;
    }
    pTable->table[startNode].distance = 0;

    // Create a priority queue (min heap)
    MinHeap *priorityQueue = createMinHeap(NUMBER_OF_VERTICES);
    if (priorityQueue == NULL) {
        return MEMORY_ALLOCATION_FAILED;
    }

    // Insert startNode into the priority queue
    MinHeapNode startNodeHeapNode = {startNode, 0};
    insertMinHeap(priorityQueue, startNodeHeapNode);

    // Dijkstra's algorithm
    while (!isEmptyMinHeap(priorityQueue)) {
        MinHeapNode minNode = extractMin(priorityQueue);
        int u = minNode.vertex;

        ListNode *adjNode = pList->adjacencyList[u];
        while (adjNode != NULL) {
            int v = adjNode->destNode;
            int weight = adjNode->weight;

            if (pTable->table[u].distance != INT_MAX &&
                pTable->table[u].distance + weight < pTable->table[v].distance) {
                pTable->table[v].distance = pTable->table[u].distance + weight;
                pTable->table[v].predecessor = u;

                MinHeapNode updatedNode = {v, pTable->table[v].distance};
                insertMinHeap(priorityQueue, updatedNode);
            }
            adjNode = adjNode->next;
        }
    }

    freeMinHeap(priorityQueue);
    return SUCCESS;
}
