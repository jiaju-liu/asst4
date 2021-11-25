#include "page_rank.h"

#include <vector>
#include <cmath>
#include <omp.h>
#include <utility>
#include <iostream>
#include <chrono>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

using namespace std;


// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double* solution, double damping, double convergence)
{

    bool converged = false;
    int numNodes = num_nodes(g);
    double equal_prob = 1.0 / numNodes;
    vector<double> score_old_array(numNodes);
    vector<double> score_new_array(numNodes);
    vector<double> score_old_divide(numNodes);
    vector<double>* score_old = &score_old_array;
    vector<double>* score_new = &score_new_array;

    vector<int> no_outgoing(numNodes);

    int list_size = 0;
    #pragma omp parallel for
    for (int i = 0; i < numNodes; i++) {
        (*score_new)[i] = equal_prob;

        // No outgoing edges
        if (!outgoing_size(g, i)) {
            int my_index = 0;
            #pragma omp atomic capture
            {my_index = list_size; list_size++;}
            no_outgoing[my_index] = i;
        } 
    }
    int chunk_size = 50;
    while (!converged) {
        vector<double> *temp = score_old;
        score_old = score_new;
        score_new = temp;
        // try next time to predivide the pagerank vectors by outgoing size
        #pragma omp parallel for 
        for (int i = 0; i < numNodes; i++) {
            int outgoing = outgoing_size(g, i);
            if (outgoing) {
                score_old_divide[i] = (*score_old)[i] / outgoing;
            }
        }
        double start_val = 0;
        #pragma imp parallel for reduction(+:start_val)
        for (int j = 0; j < list_size; j++) {
          start_val += (*score_old)[no_outgoing[j]];
        }
        start_val /= numNodes;

        // compute score_new[vi] for all nodes vi:
        #pragma omp parallel for schedule(dynamic, chunk_size)
        for (int i = 0; i < numNodes; i++) {
            const Vertex* begin = incoming_begin(g, i);
            const Vertex* end = incoming_end(g, i);
            double mySum = start_val;
            for (const Vertex* v = begin; v != end; v++) {
                // outgoing can't be 0 here
                mySum += score_old_divide[*v];
            }

            mySum = (damping * mySum) + (1.0-damping) * equal_prob;
            (*score_new)[i] = mySum;
        }
        double global_diff = 0;
        #pragma omp parallel for reduction(+:global_diff)
        for (int i = 0; i < numNodes; i++) {
            global_diff += fabs((*score_new)[i] - (*score_old)[i]);
        }
        converged = (global_diff < convergence);
    }

    #pragma omp parallel for 
    for (int i = 0; i < numNodes; i++) {
        solution[i] = (*score_new)[i];
    }
}
