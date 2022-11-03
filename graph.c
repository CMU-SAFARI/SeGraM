#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>

struct SeqNode
{
    int nodeID;
    short numOut;
    //short numIn;
    char *seq;
    int *outNodes;
    //int *inNodes;
};

void removeAll(char *str, char toRemove)
{
    int i, j;
    int len = (int)strlen(str);

    for(i=0; i<len; i++)
    {
        if(str[i] == toRemove)
        {
            for(j=i; j<len; j++)
            {
                str[j] = str[j+1];
            }

            len--;

            // If a character is removed then make sure i doesn't increments
            i--;
        }
    }
}

struct SeqNode* generateGraphFromGFA(char *filename, int *numNodes, int *numEdges)
{
    *numNodes = 0;
    *numEdges = 0;
    
    char *line = NULL;
    
    size_t len_t = 0;

    FILE *fp;
    fp = fopen(filename, "r");
    
    while ((getline(&line, &len_t, fp)) != -1)
    {
        line[strlen(line) - 1] = '\0';
        
        if (line[0] == 'S')
        {
            (*numNodes)++;
        }
    }
    
    struct SeqNode *nodes = (struct SeqNode *) malloc(*numNodes * sizeof(struct SeqNode));
    
    fseek(fp, 0, SEEK_SET);
    
    while ((getline(&line, &len_t, fp)) != -1)
    {
        line[strlen(line) - 1] = '\0';
        
        if (line[0] == 'S')
        {
            char *ptr = strtok(line, "\t");
            ptr = strtok(NULL, "\t");
            
            int ind = atoi(ptr)-1;
            
            nodes[ind].nodeID = ind;
            //nodes[ind].numIn = 0;
            nodes[ind].numOut = 0;
            nodes[ind].outNodes = NULL;
            //nodes[ind].inNodes = NULL;
            
            ptr = strtok(NULL, "\t");
            //removeAll(ptr, 'N');
            
            nodes[ind].seq = (char *)malloc((strlen(ptr)+1) * sizeof(char));
            strcpy(nodes[ind].seq, ptr);
        }
    }
    
    fseek(fp, 0, SEEK_SET);
    
    while ((getline(&line, &len_t, fp)) != -1)
    {
        line[strlen(line) - 1] = '\0';
        
        if (line[0] == 'L')
        {
            char *ptr = strtok(line, "\t");
            ptr = strtok(NULL, "\t");
            
            int from = atoi(ptr)-1;
            
            ptr = strtok(NULL, "\t");
            ptr = strtok(NULL, "\t");
            
            int to = atoi(ptr)-1;
            
            nodes[from].numOut++;
            //nodes[to].numIn++;
            
            if (nodes[from].numOut == 1)
            {
                nodes[from].outNodes = (int *) malloc(sizeof(int));
            }
            else
            {
                nodes[from].outNodes = (int *) realloc(nodes[from].outNodes, nodes[from].numOut * sizeof(int));
            }
            
            (nodes[from].outNodes)[(nodes[from].numOut)-1] = to;
            
            /*
            if (nodes[to].numIn == 1)
            {
                nodes[to].inNodes = (int *) malloc(sizeof(int));
            }
            else
            {
                nodes[to].inNodes = (int *) realloc(nodes[to].inNodes, nodes[to].numIn * sizeof(int));
            }
            
            (nodes[to].inNodes)[(nodes[to].numIn)-1] = from;
            */
            
            (*numEdges)++;
        }
    }
    
    fclose(fp);
    
    if (line)
        free(line);
    
    return nodes;
}

void getGraphStats(struct SeqNode *nodes, int numNodes, char *filename)
{
    FILE *fp;
    fp = fopen(filename, "w");
    
    int maxHopLengthGlobal = 0;
    int maxHopSeqLengthGlobal = 0;
    
    for (int i = 0; i < numNodes; i++)
    {
        int nodeID = nodes[i].nodeID;
        int seqLen = (int)strlen(nodes[i].seq);
        
        //fprintf(stderr, "%d - %d - %s\n", nodeID, seqLen, nodes[i].seq);
        
        int maxHopLength = 0;
        for (int j = 0; j < nodes[i].numOut; j++)
        {
            if (((nodes[i].outNodes)[j] - nodeID) > maxHopLength)
            {
                maxHopLength = ((nodes[i].outNodes)[j] - nodeID);
            }
        }
        
        int maxHopSeqLength = 0;
        for (int j = 0; j < nodes[i].numOut; j++)
        {
            int curHopSeqLength = 0;
            for (int k = nodeID+1; k < (nodes[i].outNodes)[j] ; k++)
            {
                curHopSeqLength += strlen(nodes[k].seq);
            }
            
            if (curHopSeqLength > maxHopSeqLength)
            {
                maxHopSeqLength = curHopSeqLength;
            }
        }
        
        fprintf(fp, "id: %d\t#in: %d\t#out: %d\tlen: %d\tmaxHop: %d\tmaxHopSeq: %d\n", i, nodes[i].numOut, nodes[i].numOut, seqLen, maxHopLength, maxHopSeqLength);
        
        if (maxHopLength > maxHopLengthGlobal)
            maxHopLengthGlobal = maxHopLength;
        if (maxHopSeqLength > maxHopSeqLengthGlobal)
            maxHopSeqLengthGlobal = maxHopSeqLength;
    }
    
    fprintf(stderr, "%s\nmaxHop: %d\tmaxHopSeq: %d\n", filename, maxHopLengthGlobal, maxHopSeqLengthGlobal);
    
    fclose(fp);
}

void deleteGraph(struct SeqNode *nodes, int numNodes)
{
    for (int i = 0; i < numNodes; i++)
    {
        /*
        if (nodes[i].inNodes)
        {
            free(nodes[i].inNodes);
        }
        */
        
        if (nodes[i].outNodes)
        {
            free(nodes[i].outNodes);
        }
        free(nodes[i].seq);
    }
    free(nodes);
}

