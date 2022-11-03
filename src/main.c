#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <ctype.h>
#include <unistd.h>
#include "graph.c"
#include "align.c"

int main () {
    
    int numNodes, numEdges;
    
    struct SeqNode *nodes = generateGraphFromGFA("test.gfa", &numNodes, &numEdges);
    
    bitalign_aligner(nodes, "ACGTCATGCAGTCGTAACGTAGTCGTCACAGTCAGTCGTAGCTAT", 0, 0, 3, 41, 0, 5, 200, 0, 1, 1, 1, 1);
    bitalign_aligner(nodes, "ACGTCATGCAGTCGTAACGTAGTCGTCACAGTCAGTCGTAGCTAGTA", 0, 0, 3, 41, 1, 5, 200, 0, 1, 1, 1, 1);
    
    deleteGraph(nodes, numNodes);
     
    return(0);
    
    
}
