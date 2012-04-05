#ifndef MUMMERGPU_GOLD_H
#define MUMMERGPU_GOLD_H

#include "common.h"

extern "C"
void computeGold(MatchResults* results,
                 char* refstr,
                 char* queries,
                 int* queryAddrs,
                 int* queryLengths,
                 PixelOfNode* nodeTexture,
                 PixelOfChildren* childrenTexture,
                 int numQueries,
                 int match_length,
                 int rc);

#endif
