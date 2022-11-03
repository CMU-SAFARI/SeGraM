#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <ctype.h>

unsigned long long *generatePatternBitmasksACGT(char* pattern, int m)
{
    int count = ceil(m/64.0);
        
    int len = 4*count; // A,C,G,T
        
    unsigned long long *patternBitmasks = (unsigned long long *) malloc(len * sizeof(unsigned long long));
        
    unsigned long long max = ULLONG_MAX;

    // Initialize the pattern bitmasks
    for (int i=0; i < len; i++)
    {
        patternBitmasks[i] = max;
    }

    // Update the pattern bitmasks
    int index;
    for (int i=0; i < m; i++)
    {
        index = count - ((m-i-1) / 64) - 1;
        if ((pattern[i] == 'A') || (pattern[i] == 'a'))
        {
            patternBitmasks[0*count + index] &= ~(1ULL << ((m-i-1) % 64));
        }
        else if ((pattern[i] == 'C') || (pattern[i] == 'c'))
        {
            patternBitmasks[1*count + index] &= ~(1ULL << ((m-i-1) % 64));
        }
        else if ((pattern[i] == 'G') || (pattern[i] == 'g'))
        {
            patternBitmasks[2*count + index] &= ~(1ULL << ((m-i-1) % 64));
        }
        else if ((pattern[i] == 'T') || (pattern[i] == 't'))
        {
            patternBitmasks[3*count + index] &= ~(1ULL << ((m-i-1) % 64));
        }
    }
    
    return patternBitmasks;
}

void bitalignTB(unsigned long long ***tracebackMatrix, int count, int W, int O, int m, int n, int minError, int *ed, int *textConsumed, int *patternConsumed, unsigned long long mask, char *lastChar, char *lastChar2, char *lastChar3, int *charCount, int *charCount2, int *charCount3,  char *CIGARstr, char *CIGARstr2, char *MD, char *text, char *subText, unsigned long long *patternBitmasks, struct SeqNode *nodes, int startNode, int endNode, int nodeStartEnd[endNode-startNode+1][3], int subOriginalNode[n], bool subIsNodeEnd[n], int textCur, int textEnd, int hopLimit, int *countM, int *countS, int *countD, int *countI, int *countOpen, int *countExtend, bool *isFirst, int scoreS, int scoreOpen, int scoreExtend)
{
    int ind;
    
    int curPattern = m-1;
    int curText = 0;
    int curError = minError;
    
    *textConsumed = 0;
    *patternConsumed = 0;
    
    while ((*textConsumed < (W-O)) && (*patternConsumed < (W-O)) && (curPattern >= 0) && (curError >= 0))
    {
        char c = subText[curText];
        
        unsigned long long curBitmask[count];
        
        if ((c == 'A') || (c == 'a'))
        {
            for (int a=0; a<count; a++)
            {
                curBitmask[a] = patternBitmasks[0*count + a];
            }
        }
        else if ((c == 'C') || (c == 'c'))
        {
            for (int a=0; a<count; a++)
            {
                curBitmask[a] = patternBitmasks[1*count + a];
            }
        }
        else if ((c == 'G') || (c == 'g'))
        {
            for (int a=0; a<count; a++)
            {
                curBitmask[a] = patternBitmasks[2*count + a];
            }
        }
        else if ((c == 'T') || (c == 't'))
        {
            for (int a=0; a<count; a++)
            {
                curBitmask[a] = patternBitmasks[3*count + a];
            }
        }
        
        ind = count - (curPattern / 64) - 1;
        
        if (subIsNodeEnd[curText])
        {
            int numOut = nodes[subOriginalNode[curText]].numOut;
            
            int relOffset[numOut];
            
            for (int g=0; g < numOut; g++)
            {
                int outNodeID = nodes[subOriginalNode[curText]].outNodes[g];
                int relStartOffset = -1;
                
                for (int i = 0; i < (endNode-startNode+1); i++)
                {
                    if (nodeStartEnd[i][0] == outNodeID)
                    {
                        relStartOffset = nodeStartEnd[i][1];
                        break;
                    }
                }
                
                if ((relStartOffset >= textCur) & (relStartOffset < textEnd))
                {
                    if (relStartOffset-textCur - curText <= hopLimit)
                    {
                        relOffset[g] = relStartOffset;
                    }
                    else
                    {
                        relOffset[g] = -1;
                        continue;
                    }
                }
                else
                {
                    relOffset[g] = -1;
                    continue;
                }
            }
            
            unsigned long long deletion[numOut][count], substitution[numOut][count], match[numOut][count], insertion[numOut][count];
            
            if (curError < 1)
            {
                for (int g = 0; g < numOut; g++)
                {
                    if (relOffset[g] != -1)
                    {
                        match[g][0] = tracebackMatrix[relOffset[g]][curError][0] << 1;
                        for (int a=1; a<count; a++)
                        {
                            match[g][a-1] |= (tracebackMatrix[relOffset[g]][curError][a] >> 63);
                            match[g][a] = tracebackMatrix[relOffset[g]][curError][a] << 1;
                        }
                        for (int a=0; a<count; a++)
                        {
                            match[g][a] |= curBitmask[a];
                        }
                        
                        for (int a=0; a<count; a++)
                        {
                            substitution[g][a] = ULLONG_MAX;
                            deletion[g][a] = ULLONG_MAX;
                            insertion[g][a] = ULLONG_MAX;
                        }
                    }
                }
            }
            else
            {
                for (int g = 0; g < numOut; g++)
                {
                    if (relOffset[g] != -1)
                    {
                        for (int a=0; a<count; a++)
                        {
                            deletion[g][a] = tracebackMatrix[relOffset[g]][curError-1][a];
                        }
                        
                        substitution[g][0] = deletion[g][0] << 1;
                        for (int a=1; a<count; a++)
                        {
                            substitution[g][a-1] |= (deletion[g][a] >> 63);
                            substitution[g][a] = deletion[g][a] << 1;
                        }
                        
                        match[g][0] = tracebackMatrix[relOffset[g]][curError][0] << 1;
                        for (int a=1; a<count; a++)
                        {
                            match[g][a-1] |= (tracebackMatrix[relOffset[g]][curError][a] >> 63);
                            match[g][a] = tracebackMatrix[relOffset[g]][curError][a] << 1;
                        }
                        
                        for (int a=0; a<count; a++)
                        {
                            match[g][a] |= curBitmask[a];
                        }
                    }
                    else
                    {
                        for (int a=0; a<count; a++)
                        {
                            substitution[g][a] = ULLONG_MAX;
                            deletion[g][a] = ULLONG_MAX;
                            match[g][a] = ULLONG_MAX;
                        }
                    }
                    
                    insertion[g][0] = tracebackMatrix[curText][curError-1][0] << 1;
                    for (int a=1; a<count; a++)
                    {
                        insertion[g][a-1] |= (tracebackMatrix[curText][curError-1][a] >> 63);
                        insertion[g][a] = tracebackMatrix[curText][curError-1][a] << 1;
                    }
                }
            }
            
            for (int g = 0; g < numOut; g++)
            {
                // affine-deletion
                if (*lastChar=='D' && ((deletion[g][ind] & mask) == 0))
                {
                    if (*lastChar == 'D')
                    {
                        *charCount += 1;
                        *countExtend += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                            sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                        }
                        *charCount = 1;
                        *lastChar = 'D';
                        *countOpen += 1;
                    }
                    if (*lastChar2 == 'D')
                    {
                        *charCount2 += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                            sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                        }
                        *charCount2 = 1;
                        *lastChar2 = 'D';
                    }
                    if (*lastChar3 == 'M')
                    {
                        sprintf(MD, "%s%d^%c", MD, *charCount3, text[curText]);
                        *lastChar3 = 'D';
                        *charCount3 = 0;
                    }
                    else if (*lastChar3 == 'D')
                    {
                        sprintf(MD, "%s%c", MD, text[curText]);
                        *lastChar3 = 'D';
                        *charCount3 = 0;
                    }
                    else
                    {
                        sprintf(MD, "%s^%c", MD, text[curText]);
                        *lastChar3 = 'D';
                        *charCount3 = 0;
                    }
                    *textConsumed += relOffset[g]-curText;
                    curText = relOffset[g];
                    curError -= 1;
                    *countD += 1;
                    *ed += 1;
                    break;
                }
                // affine-insertion
                else if (*lastChar=='I' && ((insertion[g][ind] & mask) == 0))
                {
                    if (*lastChar == 'I')
                    {
                        *charCount += 1;
                        *countExtend += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                            sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                        }
                        *charCount = 1;
                        *lastChar = 'I';
                        *countOpen += 1;
                    }
                    if (*lastChar2 == 'I')
                    {
                        *charCount2 += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                            sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                        }
                        *charCount2 = 1;
                        *lastChar2 = 'I';
                    }
                    curPattern -= 1;
                    curError -= 1;
                    mask = 1ULL << (curPattern % 64);
                    *countI += 1;
                    *ed += 1;
                    *patternConsumed += 1;
                    break;
                }
                // match
                else if ((match[g][ind] & mask) == 0)
                {
                    if (*lastChar == 'M')
                    {
                        *charCount += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                            sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                        }
                        *charCount = 1;
                        *lastChar = 'M';
                    }
                    if (*lastChar2 == 'M')
                    {
                        *charCount2 += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                            sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                        }
                        *charCount2 = 1;
                        *lastChar2 = 'M';
                    }
                    if (*lastChar3 == 'M')
                    {
                        *charCount3 += 1;
                    }
                    else
                    {
                        *charCount3 = 1;
                        *lastChar3 = 'M';
                    }
                    *textConsumed += relOffset[g]-curText;
                    curText = relOffset[g];
                    *patternConsumed += 1;
                    curPattern -= 1;
                    mask = 1ULL << (curPattern % 64);
                    *countM += 1;
                    break;
                }
                // substitution
                else if ((substitution[g][ind] & mask) == 0)
                {
                    if (*lastChar == 'S')
                    {
                        *charCount += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                            sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                        }
                        *charCount = 1;
                        *lastChar = 'S';
                    }
                    if (*lastChar2 == 'M')
                    {
                        *charCount2 += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                            sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                        }
                        *charCount2 = 1;
                        *lastChar2 = 'M';
                    }
                    if (*lastChar3 == 'M')
                    {
                        sprintf(MD, "%s%d%c", MD, *charCount3, text[curText]);
                        *lastChar3 = 'S';
                        *charCount3 = 0;
                    }
                    else
                    {
                        sprintf(MD, "%s%c", MD, text[curText]);
                        *lastChar3 = 'S';
                        *charCount3 = 0;
                    }
                    *textConsumed += relOffset[g]-curText;
                    curText = relOffset[g];
                    curPattern -= 1;
                    curError -= 1;
                    mask = 1ULL << (curPattern % 64);
                    *countS += 1;
                    *ed += 1;
                    *patternConsumed += 1;
                    break;
                }
                // deletion
                else if ((deletion[g][ind] & mask) == 0)
                {
                    if (*lastChar == 'D')
                    {
                        *charCount += 1;
                        *countExtend += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                            sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                        }
                        *charCount = 1;
                        *lastChar = 'D';
                        *countOpen += 1;
                    }
                    if (*lastChar2 == 'D')
                    {
                        *charCount2 += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                            sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                        }
                        *charCount2 = 1;
                        *lastChar2 = 'D';
                    }
                    if (*lastChar3 == 'M')
                    {
                        sprintf(MD, "%s%d^%c", MD, *charCount3, text[curText]);
                        *lastChar3 = 'D';
                        *charCount3 = 0;
                    }
                    else if (*lastChar3 == 'D')
                    {
                        sprintf(MD, "%s%c", MD, text[curText]);
                        *lastChar3 = 'D';
                        *charCount3 = 0;
                    }
                    else
                    {
                        sprintf(MD, "%s^%c", MD, text[curText]);
                        *lastChar3 = 'D';
                        *charCount3 = 0;
                    }
                    *textConsumed += relOffset[g]-curText;
                    curText = relOffset[g];
                    curError -= 1;
                    *countD += 1;
                    *ed += 1;
                    break;
                }
                // insertion
                else if ((insertion[g][ind] & mask) == 0)
                {
                    if (*lastChar == 'I')
                    {
                        *charCount += 1;
                        *countExtend += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                            sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                        }
                        *charCount = 1;
                        *lastChar = 'I';
                        *countOpen += 1;
                    }
                    if (*lastChar2 == 'I')
                    {
                        *charCount2 += 1;
                    }
                    else
                    {
                        if (!*isFirst)
                        {
                            sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                            sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                        }
                        *charCount2 = 1;
                        *lastChar2 = 'I';
                    }
                    curPattern -= 1;
                    curError -= 1;
                    mask = 1ULL << (curPattern % 64);
                    *countI += 1;
                    *ed += 1;
                    *patternConsumed += 1;
                    break;
                }
            }
        }
        else
        {
            unsigned long long deletion[count], substitution[count], match[count], insertion[count];
            
            if (curError < 1)
            {
                match[0] = tracebackMatrix[curText+1][curError][0] << 1;
                for (int a=1; a<count; a++)
                {
                    match[a-1] |= (tracebackMatrix[curText+1][curError][a] >> 63);
                    match[a] = tracebackMatrix[curText+1][curError][a] << 1;
                }
                for (int a=0; a<count; a++)
                {
                    match[a] |= curBitmask[a];
                }
                
                for (int a=0; a<count; a++)
                {
                    substitution[a] = ULLONG_MAX;
                    deletion[a] = ULLONG_MAX;
                    insertion[a] = ULLONG_MAX;
                }
            }
            else
            {
                for (int a=0; a<count; a++)
                {
                    deletion[a] = tracebackMatrix[curText+1][curError-1][a];
                }
                
                substitution[0] = deletion[0] << 1;
                for (int a=1; a<count; a++)
                {
                    substitution[a-1] |= (deletion[a] >> 63);
                    substitution[a] = deletion[a] << 1;
                }
                
                insertion[0] = tracebackMatrix[curText][curError-1][0] << 1;
                for (int a=1; a<count; a++)
                {
                    insertion[a-1] |= (tracebackMatrix[curText][curError-1][a] >> 63);
                    insertion[a] = tracebackMatrix[curText][curError-1][a] << 1;
                }
                
                match[0] = tracebackMatrix[curText+1][curError][0] << 1;
                for (int a=1; a<count; a++)
                {
                    match[a-1] |= (tracebackMatrix[curText+1][curError][a] >> 63);
                    match[a] = tracebackMatrix[curText+1][curError][a] << 1;
                }
                
                for (int a=0; a<count; a++)
                {
                    match[a] |= curBitmask[a];
                }
            }
            
            // affine-deletion
            if (*lastChar=='D' && ((deletion[ind] & mask) == 0))
            {
                curText += 1;
                curError -= 1;
                if (*lastChar == 'D')
                {
                    *charCount += 1;
                    *countExtend += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                        sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                    }
                    *charCount = 1;
                    *lastChar = 'D';
                    *countOpen += 1;
                }
                if (*lastChar2 == 'D')
                {
                    *charCount2 += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                        sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                    }
                    *charCount2 = 1;
                    *lastChar2 = 'D';
                }
                if (*lastChar3 == 'M')
                {
                    sprintf(MD, "%s%d^%c", MD, *charCount3, text[curText-1]);
                    *lastChar3 = 'D';
                    *charCount3 = 0;
                }
                else if (*lastChar3 == 'D')
                {
                    sprintf(MD, "%s%c", MD, text[curText-1]);
                    *lastChar3 = 'D';
                    *charCount3 = 0;
                }
                else
                {
                    sprintf(MD, "%s^%c", MD, text[curText-1]);
                    *lastChar3 = 'D';
                    *charCount3 = 0;
                }
                *countD += 1;
                *ed += 1;
                *textConsumed += 1;
            }
            // affine-insertion
            else if (*lastChar=='I' && ((insertion[ind] & mask) == 0))
            {
                curPattern -= 1;
                curError -= 1;
                mask = 1ULL << (curPattern % 64);
                if (*lastChar == 'I')
                {
                    *charCount += 1;
                    *countExtend += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                        sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                    }
                    *charCount = 1;
                    *lastChar = 'I';
                    *countOpen += 1;
                }
                if (*lastChar2 == 'I')
                {
                    *charCount2 += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                        sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                    }
                    *charCount2 = 1;
                    *lastChar2 = 'I';
                }
                *countI += 1;
                *ed += 1;
                *patternConsumed += 1;
            }
            // match
            else if ((match[ind] & mask) == 0)
            {
                curText += 1;
                curPattern -= 1;
                mask = 1ULL << (curPattern % 64);
                if (*lastChar == 'M')
                {
                    *charCount += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                        sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                    }
                    *charCount = 1;
                    *lastChar = 'M';
                }
                if (*lastChar2 == 'M')
                {
                    *charCount2 += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                        sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                    }
                    *charCount2 = 1;
                    *lastChar2 = 'M';
                }
                if (*lastChar3 == 'M')
                {
                    *charCount3 += 1;
                }
                else
                {
                    *charCount3 = 1;
                    *lastChar3 = 'M';
                }
                *countM += 1;
                *textConsumed += 1;
                *patternConsumed += 1;
            }
            // substitution
            else if ((substitution[ind] & mask) == 0)
            {
                curText += 1;
                curPattern -= 1;
                curError -= 1;
                mask = 1ULL << (curPattern % 64);
                if (*lastChar == 'S')
                {
                    *charCount += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                        sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                    }
                    *charCount = 1;
                    *lastChar = 'S';
                }
                if (*lastChar2 == 'M')
                {
                    *charCount2 += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                        sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                    }
                    *charCount2 = 1;
                    *lastChar2 = 'M';
                }
                if (*lastChar3 == 'M')
                {
                    sprintf(MD, "%s%d%c", MD, *charCount3, text[curText-1]);
                    *lastChar3 = 'S';
                    *charCount3 = 0;
                }
                else
                {
                    sprintf(MD, "%s%c", MD, text[curText-1]);
                    *lastChar3 = 'S';
                    *charCount3 = 0;
                }
                *countS += 1;
                *ed += 1;
                *textConsumed += 1;
                *patternConsumed += 1;
            }
            // deletion
            else if ((deletion[ind] & mask) == 0)
            {
                curText += 1;
                curError -= 1;
                if (*lastChar == 'D')
                {
                    *charCount += 1;
                    *countExtend += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                        sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                    }
                    *charCount = 1;
                    *lastChar = 'D';
                    *countOpen += 1;
                }
                if (*lastChar2 == 'D')
                {
                    *charCount2 += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                        sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                    }
                    *charCount2 = 1;
                    *lastChar2 = 'D';
                }
                if (*lastChar3 == 'M')
                {
                    sprintf(MD, "%s%d^%c", MD, *charCount3, text[curText-1]);
                    *lastChar3 = 'D';
                    *charCount3 = 0;
                }
                else if (*lastChar3 == 'D')
                {
                    sprintf(MD, "%s%c", MD, text[curText-1]);
                    *lastChar3 = 'D';
                    *charCount3 = 0;
                }
                else
                {
                    sprintf(MD, "%s^%c", MD, text[curText-1]);
                    *lastChar3 = 'D';
                    *charCount3 = 0;
                }
                *countD += 1;
                *ed += 1;
                *textConsumed += 1;
            }
            // insertion
            else if ((insertion[ind] & mask) == 0)
            {
                curPattern -= 1;
                curError -= 1;
                mask = 1ULL << (curPattern % 64);
                if (*lastChar == 'I')
                {
                    *charCount += 1;
                    *countExtend += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr, "%s%d", CIGARstr, *charCount);
                        sprintf(CIGARstr, "%s%c", CIGARstr, *lastChar);
                    }
                    *charCount = 1;
                    *lastChar = 'I';
                    *countOpen += 1;
                }
                if (*lastChar2 == 'I')
                {
                    *charCount2 += 1;
                }
                else
                {
                    if (!*isFirst)
                    {
                        sprintf(CIGARstr2, "%s%d", CIGARstr2, *charCount2);
                        sprintf(CIGARstr2, "%s%c", CIGARstr2, *lastChar2);
                    }
                    *charCount2 = 1;
                    *lastChar2 = 'I';
                }
                *countI += 1;
                *ed += 1;
                *patternConsumed += 1;
            }
        }
        
        *isFirst = false;
    }
}

void bitalignDC(struct SeqNode *nodes, char *pattern, int startNode, int startOffset, int endNode, int endOffset, int k, int hopLimit, int W, int O, int scoreM, int scoreS, int scoreOpen, int scoreExtend)
{
    int n = 0;
    
    for (int i = startNode; i <= endNode; i++)
    {
        if (i == startNode)
            n += strlen(nodes[i].seq)-startOffset;
        else if (i == endNode)
            n += endOffset+1;
        else
            n += strlen(nodes[i].seq);
    }
    
    int m = (int)strlen(pattern);
    
    char *text = (char *)malloc((n+1) * sizeof(char));
    strcpy(text, "");
    
    bool isNodeStart[n];
    bool isNodeEnd[n];
    int  originalNode[n];
    
    int nodeStartEnd[endNode-startNode+1][3];       // <nodeID, startTextOffset, endTextOffset>
    
    int ind = 0;
    int start = 0;
    int end = 0;
    
    int nodeInd = 0;
    
    for (int i = startNode; i <= endNode; i++)
    {
        if (i == startNode)
        {
            start = startOffset;
            end = (int)strlen(nodes[i].seq);
        }
        else if (i == endNode)
        {
            start = 0;
            end = endOffset+1;
        }
        else
        {
            start = 0;
            end = (int)strlen(nodes[i].seq);
        }
        
        //fprintf(stderr, "%d - %d - %s\n", start, end, nodes[i].seq);
        
        nodeStartEnd[nodeInd][0] = i;   //nodeID
        
        for (int j = start; j < end; j++)
        {
            strncat(text, &(nodes[i].seq)[j], 1);
            isNodeStart[ind] = false;
            isNodeEnd[ind] = false;
            if (j == start)
            {
                isNodeStart[ind] = true;
                nodeStartEnd[nodeInd][1] = ind;   //startTextOffset
            }
            if (j == end-1)
            {
                isNodeEnd[ind] = true;
                nodeStartEnd[nodeInd][2] = ind;   //endTextOffset
            }
            originalNode[ind] = i;
            ind++;
        }
        nodeInd++;
    }
    text[n] = '\0';
    
    //fprintf(stderr, "%s - %d\n", text, n);
    
    /*
    for (int i = 0; i < ind; i++)
    {
        fprintf(stderr, "%d (%d): %s - %s\n", i, originalNode[i], isNodeStart[i]?"true":"false", isNodeEnd[i]?"true":"false");
    }
    
    for (int i = 0; i < nodeInd; i++)
    {
        fprintf(stderr, "%d (%d) - %d - %d\n", nodeStartEnd[i][0], i, nodeStartEnd[i][1], nodeStartEnd[i][2]);
    }
    */
    
    unsigned long long max = ULLONG_MAX;
    
    int ed = 0;
    char CIGARstr[m+k];
    strcpy(CIGARstr, "");
    char CIGARstr2[m+k];
    strcpy(CIGARstr2, "");
    char MD[m+k];
    strcpy(MD, "");
    int charCount = 0;
    int charCount2 = 0;
    int charCount3 = 0;
    char lastChar = '0';
    char lastChar2 = '0';
    char lastChar3 = '0';
    
    bool isFirst = true;
    
    int countM = 0;
    int countS = 0;
    int countI = 0;
    int countD = 0;
    int countOpen = 0;
    int countExtend = 0;
    
    int textCur = 0;
    int patternCur = 0;
    int textEnd = 0;
    int patternEnd = 0;
    
    int newError = 0;
    if (W > k)
    {
        newError = k;
    }
    else
    {
        newError = W*0.3;
    }
    
    unsigned long long ***tracebackMatrix;
    tracebackMatrix = (unsigned long long ***) malloc(W * sizeof(unsigned long long **));
    for (int x=0; x < W; x++)
    {
        tracebackMatrix[x] = (unsigned long long **) malloc((newError + 1) * sizeof(unsigned long long *));
        for (int y=0; y < (newError + 1); y++)
        {
            tracebackMatrix[x][y] = (unsigned long long *) malloc(ceil(W/64.0) * sizeof(unsigned long long));
        }
    }
    
    while ((patternCur < m) && (textCur < n))
    {
        textEnd = textCur + W;
        if (textEnd > n)
        {
            textEnd = n;
        }
        patternEnd = patternCur + W;
        if (patternEnd > m)
        {
            patternEnd = m;
        }
        
        int n_sub = textEnd - textCur;
        int m_sub = patternEnd - patternCur;
        
        char subText[n_sub+1];
        char subPattern[m_sub+1];
        
        bool subIsNodeStart[n_sub];
        bool subIsNodeEnd[n_sub];
        int subOriginalNode[n_sub];
        
        int y = 0;
        for (int x=textCur; x < textEnd; x++)
        {
            subText[y] = text[x];
            subIsNodeStart[y] = isNodeStart[x];
            subIsNodeEnd[y] = isNodeEnd[x];
            subOriginalNode[y++] = originalNode[x];
        }
        subText[y] = '\0';
        
        y = 0;
        for (int x=patternCur; x < patternEnd; x++)
        {
            subPattern[y++] = pattern[x];
        }
        subPattern[y] = '\0';
        
        unsigned long long *patternBitmasks = generatePatternBitmasksACGT(subPattern, m_sub);
        
        int count = ceil(m_sub/64.0);
        int rem = m_sub % 64;
        
        unsigned long long max1;
        
        if (rem == 0)
        {
            //max1 = MSB_MASK;
            max1 = 1ULL << 63;
        }
        else
        {
            max1 = 1ULL << (rem-1);
        }
        
        // Initialize the bit arrays R
        int len1 = (newError+1) * count;
        unsigned long long R[len1];
        
        for (int i=0; i < len1; i++)
        {
            R[i] = max;
        }
        
        unsigned long long oldR[len1];
        
        unsigned long long substitution[count], insertion[count], match[count], deletion[count], curBitmask[count];
        
        // now traverse the text in opposite direction (i.e., forward), generate partial tracebacks at each checkpoint
        for (int i=n_sub-1; i >= 0; i--)
        {
            char c = subText[i];
            
            if ((c == 'A') || (c == 'a') || (c == 'C') || (c == 'c') || (c == 'G') || (c == 'g') || (c == 'T') || (c == 't'))
            {
                // copy the content of R to oldR
                for (int itR=0; itR<len1; itR++)
                {
                    oldR[itR] = R[itR];
                }
                
                if ((c == 'A') || (c == 'a'))
                {
                    for (int a=0; a<count; a++)
                    {
                        curBitmask[a] = patternBitmasks[0*count + a];
                    }
                }
                else if ((c == 'C') || (c == 'c'))
                {
                    for (int a=0; a<count; a++)
                    {
                        curBitmask[a] = patternBitmasks[1*count + a];
                    }
                }
                else if ((c == 'G') || (c == 'g'))
                {
                    for (int a=0; a<count; a++)
                    {
                        curBitmask[a] = patternBitmasks[2*count + a];
                    }
                }
                else if ((c == 'T') || (c == 't'))
                {
                    for (int a=0; a<count; a++)
                    {
                        curBitmask[a] = patternBitmasks[3*count + a];
                    }
                }
                
                if (subIsNodeEnd[i])
                {
                    int numOut = nodes[subOriginalNode[i]].numOut;
                    
                    int relOffset[numOut];
                    
                    for (int g=0; g < numOut; g++)
                    {
                        int outNodeID = nodes[subOriginalNode[i]].outNodes[g];
                        int relStartOffset = -1;
                        
                        for (int i = 0; i < nodeInd; i++)
                        {
                            if (nodeStartEnd[i][0] == outNodeID)
                            {
                                relStartOffset = nodeStartEnd[i][1];
                                break;
                            }
                        }
                        
                        if ((relStartOffset >= textCur) & (relStartOffset < textEnd))
                        {
                            if (relStartOffset-textCur - i <= hopLimit)
                            {
                                relOffset[g] = relStartOffset;
                            }
                            else
                            {
                                relOffset[g] = -1;
                                continue;
                            }
                        }
                        else
                        {
                            relOffset[g] = -1;
                            continue;
                        }
                    }
                    
                    for (int itR=0; itR<len1; itR++)
                    {
                        R[itR] = max;
                    }
                    
                    unsigned long long hopMatch[count], hopSubstitution[count], hopDeletion[count];
                    
                    // calculating R[0] by ANDing match BVs for all hops
                    for (int g = 0; g < numOut; g++)
                    {
                        if (relOffset[g] != -1)
                        {
                            hopMatch[0] = tracebackMatrix[relOffset[g]][0][0] << 1;
                            for (int a=1; a<count; a++)
                            {
                                hopMatch[a-1] |= (tracebackMatrix[relOffset[g]][0][a] >> 63);
                                hopMatch[a] = tracebackMatrix[relOffset[g]][0][a] << 1;
                            }
                            for (int a=0; a<count; a++)
                            {
                                hopMatch[a] |= curBitmask[a];
                                R[a] &= hopMatch[a];
                            }
                        }
                    }
                    for (int a=0; a<count; a++)
                    {
                        tracebackMatrix[i][0][a] = R[a];
                    }
        
                    // calculating R[d] for all hops
                    for (int d=1; d <= newError; d++)
                    {
                        for (int g = 0; g < numOut; g++)
                        {
                            if (relOffset[g] != -1)
                            {
                                for (int a=0; a<count; a++)
                                {
                                    hopDeletion[a] = tracebackMatrix[relOffset[g]][d-1][a];
                                }
                                
                                hopSubstitution[0] = hopDeletion[0] << 1;
                                for (int a=1; a<count; a++)
                                {
                                    hopSubstitution[a-1] |= (hopDeletion[a] >> 63);
                                    hopSubstitution[a] = hopDeletion[a] << 1;
                                }
                                
                                hopMatch[0] = tracebackMatrix[relOffset[g]][d][0] << 1;
                                for (int a=1; a<count; a++)
                                {
                                    hopMatch[a-1] |= (tracebackMatrix[relOffset[g]][d][a] >> 63);
                                    hopMatch[a] = tracebackMatrix[relOffset[g]][d][a] << 1;
                                }
                                
                                for (int a=0; a<count; a++)
                                {
                                    hopMatch[a] |= curBitmask[a];
                                }
                                
                                for (int a=0; a<count; a++)
                                {
                                    R[d*count+a] &= hopDeletion[a] & hopSubstitution[a] & hopMatch[a];
                                }
                            }
                        }
                        
                        int index = (d-1) * count;
                        
                        insertion[0] = R[index] << 1;
                        for (int a=1; a<count; a++)
                        {
                            insertion[a-1] |= (R[index+a] >> 63);
                            insertion[a] = R[index+a] << 1;
                            R[d*count+a] &= insertion[a];
                        }
                        
                        for (int a=0; a<count; a++)
                        {
                            tracebackMatrix[i][d][a] = R[d*count+a];
                        }
                    }
                }
                else
                {
                    // update R[0] by first shifting oldR[0] and then ORing with pattern bitmask
                    R[0] = oldR[0] << 1;
                    for (int a=1; a<count; a++)
                    {
                        R[a-1] |= (oldR[a] >> 63);
                        R[a] = oldR[a] << 1;
                    }
                    for (int a=0; a<count; a++)
                    {
                        R[a] |= curBitmask[a];
                    }
                    
                    for (int a=0; a<count; a++)
                    {
                        tracebackMatrix[i][0][a] = R[a];
                    }
                    
                    for (int d=1; d <= newError; d++)
                    {
                        int index = (d-1) * count;
                        
                        for (int a=0; a<count; a++)
                        {
                            deletion[a] = oldR[index + a];
                        }
                        
                        substitution[0] = deletion[0] << 1;
                        for (int a=1; a<count; a++)
                        {
                            substitution[a-1] |= (deletion[a] >> 63);
                            substitution[a] = deletion[a] << 1;
                        }
                        
                        insertion[0] = R[index] << 1;
                        for (int a=1; a<count; a++)
                        {
                            insertion[a-1] |= (R[index+a] >> 63);
                            insertion[a] = R[index+a] << 1;
                        }
                        
                        index += count;
                        
                        match[0] = oldR[index] << 1;
                        for (int a=1; a<count; a++)
                        {
                            match[a-1] |= (oldR[index+a] >> 63);
                            match[a] = oldR[index+a] << 1;
                        }
                        
                        for (int a=0; a<count; a++)
                        {
                            match[a] |= curBitmask[a];
                        }
                        
                        for (int a=0; a<count; a++)
                        {
                            R[index+a] = deletion[a] & substitution[a] & insertion[a] & match[a];
                            tracebackMatrix[i][d][a] = R[index+a];
                        }
                    }
                }
            }
        }
        
        int minError = -1;
        unsigned long long mask = max1;
        
        for (int t=0; t <= newError; t++)
        {
            if ((R[t*count] & mask) == 0)
            {
                minError = t;
                break;
            }
        }
        
        //fprintf(stderr, "%d", minError);
        
        if (minError == -1)
        {
            printf("No alignment found!\t");
            
            free(patternBitmasks);
            
            for (int x=0; x < W; x++)
            {
                for (int y=0; y < (newError + 1); y++)
                {
                    free(tracebackMatrix[x][y]);
                }
                free(tracebackMatrix[x]);
            }
            free(tracebackMatrix);
            
            return;
        }
        
        int textConsumed;
        int patternConsumed;
        
        bitalignTB(tracebackMatrix, count, W, O, m_sub, n_sub, minError, &ed, &textConsumed, &patternConsumed, mask, &lastChar, &lastChar2, &lastChar3, &charCount, &charCount2, &charCount3,  CIGARstr, CIGARstr2, MD, text, subText, patternBitmasks, nodes, startNode, endNode, nodeStartEnd, subOriginalNode, subIsNodeEnd, textCur, textEnd, hopLimit, &countM, &countS, &countD, &countI, &countOpen, &countExtend, &isFirst, scoreS, scoreOpen, scoreExtend);
        
        textCur += textConsumed;
        patternCur += patternConsumed;
        
        free(patternBitmasks);
    }
    
    for (int x=0; x < W; x++)
    {
        for (int y=0; y < (newError + 1); y++)
        {
            free(tracebackMatrix[x][y]);
        }
        free(tracebackMatrix[x]);
    }
    free(tracebackMatrix);
    
    sprintf(CIGARstr, "%s%d", CIGARstr, charCount);
    sprintf(CIGARstr, "%s%c", CIGARstr, lastChar);
    
    sprintf(CIGARstr2, "%s%d", CIGARstr2, charCount2);
    sprintf(CIGARstr2, "%s%c", CIGARstr2, lastChar2);
    
    if (lastChar3 == 'M')
    {
        sprintf(MD, "%s%d", MD, charCount3);
    }
    
    int score = countM*scoreM-countS*scoreS-countOpen*(scoreOpen+scoreExtend)-countExtend*scoreExtend;
    
    printf("ED:%d\tAS:%d\t%s\t%s\t%s\n", ed, score, CIGARstr, CIGARstr2, MD);
    printf("#M:%d\t#S:%d\t#O:%d\t#E:%d\n", countM, countS, countOpen, countExtend);
    
}

void bitalign_aligner(struct SeqNode *nodes, char *pattern, int startNode, int startOffset, int endNode, int endOffset, int k, int hopLimit, int W, int O, int scoreM, int scoreS, int scoreOpen, int scoreExtend)
{
    bitalignDC(nodes, pattern, startNode, startOffset, endNode, endOffset, k, hopLimit, W, O, scoreM, scoreS, scoreOpen, scoreExtend);
}
