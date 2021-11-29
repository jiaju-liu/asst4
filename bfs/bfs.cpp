#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <vector>
#include <iostream>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
#define CHUNK_SIZE 100
using namespace std;

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
int top_down_step(
    Graph g,
    int* distances,
    int counter)
{
    int new_counter = counter+1;
    int progress = 0;
    // maybe make this dynamic
    #pragma omp parallel for //schedule(dynamic, g->num_nodes / 1000) 
    for (int i=0; i<g->num_nodes; i++) {
        if (distances[i] == counter) {
            int myProgress = 0;

            const Vertex* start = outgoing_begin(g, i);
            const Vertex* end = outgoing_end(g, i);

            for (const Vertex* node = start; node != end; node++) {
                if (distances[*node] == NOT_VISITED_MARKER) {
                    distances[*node] = new_counter;
                    myProgress++;
                }
            }
            
            if (myProgress) {
                progress += myProgress;
            }
        }
    }
    return progress;
}

int bottom_up_step(
    Graph g,
    int* distances,
    int counter)
{
    int progress = 0;
    int new_counter = counter + 1;
    #pragma omp parallel for schedule(dynamic, g->num_nodes / 200) 
    for (int i=0; i<g->num_nodes; i++) {
        if (distances[i] == NOT_VISITED_MARKER) {
            int myProgress = 0;

            const Vertex* start = incoming_begin(g, i);
            const Vertex* end = incoming_end(g, i);

            for (const Vertex* node = start; node != end; node++) {
                int incoming = *node;
                if (distances[incoming] == counter) {
                    distances[i] = new_counter;
                    myProgress++;
                    break;
                }
            }
            if (myProgress) {
                progress += myProgress;
            }
        }
    }
    return progress;
}
// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    sol->distances[ROOT_NODE_ID] = 0;
    int count = 0;

    while (top_down_step(graph, sol->distances, count++)) {
    }
}

void bfs_bottom_up(Graph graph, solution* sol)
{

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    sol->distances[ROOT_NODE_ID] = 0;
    int count = 0;

    while (bottom_up_step(graph, sol->distances, count++)) {
    }
}

void bfs_hybrid(Graph graph, solution* sol)
{

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;
    int thresh = 0.15 *graph->num_nodes;

    sol->distances[ROOT_NODE_ID] = 0;
    int count = 0;
    int progress = top_down_step(graph, sol->distances, count++);

    while (progress) {
        if (progress > thresh) {
            progress = bottom_up_step(graph, sol->distances, count++);
        } else {
            progress = top_down_step(graph, sol->distances, count++);
        }
    }
}
