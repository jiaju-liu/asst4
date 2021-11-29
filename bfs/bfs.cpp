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
bool top_down_step(
    Graph g,
    int* distances,
    int counter)
{
    int new_counter = counter+1;
    bool progress = false;
    // maybe make this dynamic
    #pragma omp parallel for //schedule(dynamic, g->num_nodes / 1000) 
    for (int i=0; i<g->num_nodes; i++) {
        if (distances[i] == counter) {
            bool myProgress = false;

            const Vertex* start = outgoing_begin(g, i);
            const Vertex* end = outgoing_end(g, i);

            for (const Vertex* node = start; node != end; node++) {
                if (distances[*node] == NOT_VISITED_MARKER) {
                    distances[*node] = new_counter;
                    myProgress = true;
                }
            }
            
            if (myProgress) {
                progress = true;
            }
        }
    }
    return progress;
}

bool bottom_up_step(
    Graph g,
    int* distances,
    int counter)
{
    int progress = 0;
    int new_counter = counter + 1;
    //vector<int> myProgress(g->num_nodes,0);
    // maybe make this dynamic
    #pragma omp parallel for schedule(dynamic, g->num_nodes / 100) 
    for (int i=0; i<g->num_nodes; i++) {
        if (distances[i] == NOT_VISITED_MARKER) {
            int myProgress = 0;

            const Vertex* start = incoming_begin(g, i);
            const Vertex* end = incoming_end(g, i);
            //int start_edge = g->incoming_starts[i];
            //int end_edge = (i == g->num_nodes - 1)
                               //? g->num_edges
                               //: g->incoming_starts[i + 1];

            // attempt to add all neighbors to the new frontier
            //for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                //int incoming = g->outgoing_edges[neighbor];
            for (const Vertex* node = start; node != end; node++) {
                int incoming = *node;
                if (distances[incoming] == counter) {
                    distances[i] = new_counter;
                    myProgress = 1;
                    break;
                }
            }
            if (myProgress) {
                progress = 1;
            }
        }
    }
    return progress != 0;
}
// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;
    int count = 0;
    bool progress = true;

    while (top_down_step(graph, sol->distances, count++)) {

// #ifdef VERBOSE
//         double start_time = CycleTimer::currentSeconds();
// #endif

//         progress = top_down_step(graph, sol->distances, count++);

// #ifdef VERBOSE
//     double end_time = CycleTimer::currentSeconds();
//     printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
// #endif
    }
}

void bfs_bottom_up(Graph graph, solution* sol)
{

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    sol->distances[ROOT_NODE_ID] = 0;
    int count = 0;
    bool progress = true;

    while (bottom_up_step(graph, sol->distances, count++)) {

// #ifdef VERBOSE
//         double start_time = CycleTimer::currentSeconds();
// #endif

//         progress = bottom_up_step(graph, sol->distances, count++);

// #ifdef VERBOSE
//     double end_time = CycleTimer::currentSeconds();
//     printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
// #endif

    }
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
}
