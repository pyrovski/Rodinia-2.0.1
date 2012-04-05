#include <stdio.h>
#include <string.h>


#define ulong4 uint32_t
#define int2 int32_t
#define uint4 uint32_t

#include "omp.h"
#include "mummergpu.h"
// Matches are reported as a node in the suffix tree,
// plus a distance up the node's parent link for partial
// matches on the patch from the root to the node

//#define OMP
#define N_THREADS 4


static const int maxdim = 4096;

#define __VERBOSE___no

#ifdef __VERBOSE__
#define XPRINTF(...)  printf(__VA_ARGS__)
#else
#define XPRINTF(...)  do{}while(0)
#endif

#define WARP_SIZE 16

#if REORDER_TREE
#define fNID "%d,%d"
#define NID(addr) (addr & 0x0000FFFF), ((addr & 0xFFFF0000)>>16)
#define GOROOT(addr) addr = 0x00010000
//#define GOROOT(addr) addr.x = 0; addr.y = 1
#else
#define fNID "%d"
#define NID(addr) addr
#define GOROOT(addr) addr = 1
#endif


#define GETQCHAR(qrypos) queries[qrypos]
#define RESULT_SPAN 1
#define MATCH_BASE(match_coords, qryid) (MatchCoord*)match_coords + qryAddr - (qryid * (min_match_len + 1))

#define GETRCHAR(refpos) getRef(refpos, ref)



#if MERGETEX

#if TREE_ACCESS_HISTOGRAM
#define GETNODE(addr, two_level)         getMerged(nodes, childrenarr, addr,  0, NULL, NULL).node
#define GETNODEHIST(addr, two_level)     getMerged(nodes, childrenarr, addr,  0, node_hist, child_hist).node
#define GETCHILDREN(addr, two_level)     getMerged(nodes, childrenarr, addr,  1, NULL, NULL).children
#define GETCHILDRENHIST(addr, two_level) getMerged(nodes, childrenarr, addr,  1, node_hist, child_hist).children
#else
#define GETNODE(addr, two_level)         getMerged(nodes, childrenarr, addr,  0).node
#define GETNODEHIST(addr, two_level)     getMerged(nodes, childrenarr, addr,  0).node
#define GETCHILDREN(addr, two_level)     getMerged(nodes, childrenarr, addr,  1).children
#define GETCHILDRENHIST(addr, two_level) getMerged(nodes, childrenarr, addr,  1).children
#endif

#else

#if TREE_ACCESS_HISTOGRAM
#define GETNODEHIST(addr, two_level)    getNode(addr,  nodes, node_hist)
#define GETNODE(addr, two_level)        getNode(addr,  nodes, NULL)
#else
#define GETNODEHIST(addr, two_level)    getNode(addr,  nodes)
#define GETNODE(addr, two_level)        getNode(addr,  nodes)
#endif


#if TREE_ACCESS_HISTOGRAM
#define GETCHILDRENHIST(addr, two_level)    getChildren(addr,  childrenarr, child_hist)
#define GETCHILDREN(addr, two_level)        getChildren(addr,  childrenarr, NULL)
#else
#define GETCHILDRENHIST(addr, two_level)    getChildren(addr,  childrenarr)
#define GETCHILDREN(addr, two_level)        getChildren(addr,  childrenarr)
#endif

#endif

#define SHIFT_QUERIES(queries, qryAddr) queries += qryAddr
#define SET_RESULT(c, r, e, q, m, rc) set_result(c, r, e, q, m, rc)
//////////////////////////////////
/// getRef
//////////////////////////////////

static inline char getRef(int refpos, const char* ref) {
  return ref[refpos];
}


union SingleNode {
  PixelOfNode node;
  PixelOfChildren children;
};


//////////////////////////////////
/// getNode
//////////////////////////////////

inline const PixelOfNode getNode(unsigned int cur,  const PixelOfNode* nodes
#if TREE_ACCESS_HISTOGRAM
                    , int* node_hist
#endif
                   ) {
#if TREE_ACCESS_HISTOGRAM
  int id = addr2id(cur);
  if (node_hist) {
    node_hist[id]++;
  }
#endif

#if REORDER_TREE
  return *(nodes + (cur & 0x0000FFFF) + (((cur & 0xFFFF0000)>>16) * MAX_TEXTURE_DIMENSION));
#else
  return *(nodes + cur);
#endif
}

//////////////////////////////////
/// getChildren
//////////////////////////////////

inline const PixelOfChildren getChildren(unsigned int cur, const PixelOfChildren* childrenarr
#if TREE_ACCESS_HISTOGRAM
                            , int* child_hist
#endif
                           ) {
#if TREE_ACCESS_HISTOGRAM
  int id = addr2id(cur);
  if (child_hist) {
    child_hist[id]++;
  }
#endif

#if REORDER_TREE
  return *(childrenarr +  (cur & 0x0000FFFF) + (((cur & 0xFFFF0000)>>16) * MAX_TEXTURE_DIMENSION));
#else
  return *(childrenarr + cur);
#endif
}

#if MERGETEX

//////////////////////////////////
/// getMerged
//////////////////////////////////

inline SingleNode getMerged(PixelOfNode * nodes,
                     PixelOfChildren * childrenarr,
                     unsigned int cur,
                     int   getChildrenData
#if TREE_ACCESS_HISTOGRAM
                     , int* node_hist
                     , int* child_hist
#endif
                    ) {
  SingleNode n;
//	TextureAddress cur = _cur;
#if !REORDER_TREE
  //cur.x *= 2;
  unsigned int x = cur * 2;
  int useChildrenForData = 0;

  if (x >= MAX_TEXTURE_DIMENSION*MAX_TEXTURE_DIMENSION) {
    x -= MAX_TEXTURE_DIMENSION*MAX_TEXTURE_DIMENSION;
    useChildrenForData = 1;
  }

#else
  unsigned short x = cur & 0x0000FFFF;
  unsigned short y = (cur & 0xFFFF0000) >> 16;
  int useChildrenForData = 0;

  // WARNING INSANE HACK TO WORK AROUND NVCC BUG

  goto TEST;
MASK:

  x &= 0x7FF;
  x *= 2;

  goto INC;
TEST:

  if (x >= 2048) {
    useChildrenForData = 1;
  }

  goto MASK;
INC:

#endif

  x += getChildrenData;

#if !REORDER_TREE
  cur = x;
#else
  cur = (y << 16) | x;
#endif

  if (useChildrenForData) {
    n.children = getChildren(cur, childrenarr
#if TREE_ACCESS_HISTOGRAM
                             , child_hist
#endif
                            );
  } else {
    n.node =  getNode(cur, nodes
#if TREE_ACCESS_HISTOGRAM
                      , node_hist
#endif
                     );
  }
  return n;
}

#endif


//////////////////////////////////
/// set_result
//////////////////////////////////

inline void set_result(unsigned int cur,
                MatchCoord* result,
                int edge_match_length,
                int qry_match_len,
                int min_match_len,
                int rc
               ) {
  if (qry_match_len > min_match_len) {
    edge_match_length |= rc;
    result->node.data = cur;
    result->edge_match_length = edge_match_length;
  } else {
    XPRINTF("  match too short (%d < %d)\n", qry_match_len, min_match_len);
  }
}

inline unsigned int arrayToAddress(const unsigned char arr[3]) {
#if REORDER_TREE
  return (arr[0] | ((arr[2] & 0xF) << 8)) | ((arr[1] | ((arr[2] & 0xF0) << 4)) << 16);
#else
  return MK3(arr);
#endif
}

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>=(b)?(a):(b))

//! @todo OpenACC
/*! @todo what are the sizes of arrays to pass to the device?
  How long should they reside on the device?
*/
int kernel_gold(MatchResults* results,
		const char* queries,
		const PixelOfNode* nodes,
		const PixelOfChildren* childrenarr,
		const char* ref,
                const int* queryAddrs,
                const int* queryLengths,
                const int numQueries,
                const int min_match_len
#if TREE_ACCESS_HISTOGRAM
                ,int* node_hist,
                int* child_hist
#endif
		) {
  int qryid;
  //#pragma acc data copy() local()
  {
    #pragma acc kernels loop independent gang
    for (qryid = 0; qryid < numQueries; ++qryid) {
      XPRINTF("> qryid: %d\n", qryid);
      
      int qlen = queryLengths[qryid];
      int qryAddr = queryAddrs[qryid];
      
      unsigned int cur = 0;
      
      int mustmatch = 0;
      int qry_match_len = 0;
      MatchCoord* match_coords = results->h_match_coords;
      MatchCoord * result = MATCH_BASE(match_coords, qryid);
      
      SHIFT_QUERIES(queries, qryAddr);
      
      int last = qlen - min_match_len;
      int qrystart;
      //! @todo are these loop iterations independent?
      for (qrystart = 0; qrystart <= last; qrystart++, result += RESULT_SPAN){
	//PixelOfNode node;
	unsigned int node_start;
	unsigned int prev;
	
	if ((cur == 0) || (qry_match_len < 1)) {
	  // start at root of tree
	  GOROOT(cur);
	  qry_match_len = 1;
	  mustmatch = 0;
	}

	char c = GETQCHAR(qrystart + qry_match_len);

	XPRINTF("In node ("fNID"): starting with %c [%d] =>  \n",
		NID(cur), c, qry_match_len);

	unsigned int refpos = 0;
	while ((c != '\0')) {
	  XPRINTF("Next edge to follow: %c (%d)\n", c, qry_match_len);

	  const PixelOfChildren children = GETCHILDRENHIST(cur, false);
	  prev = cur;

	  switch (c) {
	  case 'A':
	    cur = arrayToAddress(children.a);
	    break;
	  case 'C':
	    cur = arrayToAddress(children.c);
	    break;
	  case 'G':
	    cur = arrayToAddress(children.g);
	    break;
	  case 'T':
	    cur = arrayToAddress(children.t);
	    break;
	  default:
	    cur = 0;
	    break;
	  };

	  //arrayToAddress(next, cur);

	  XPRINTF(" In node: ("fNID")\n", NID(cur));

	  // No edge to follow out of the node
	  if (cur == 0) {
	    XPRINTF(" no edge\n");
	    SET_RESULT(prev, result, 0, qry_match_len, min_match_len, FORWARD);

	    qry_match_len -= 1;
	    mustmatch = 0;

	    goto NEXT_SUBSTRING;
	  }

	  PixelOfNode node;
	  node = GETNODEHIST(cur, true);
	  node_start = MK3(node.start);
	  unsigned int node_end = MK3(node.end);

	  XPRINTF(" Edge coordinates: %d - %d\n", node_start, node_end);
	  {
	    int edgelen = node_end - node_start + 1;
	    int edge_matchlen = node_start + mustmatch;
	    int past_node_end = node_end + 1;
	    int dist_to_edge_end = mustmatch - edgelen;
	    if (mustmatch) {
	      refpos = min(edge_matchlen, past_node_end);
	      qry_match_len += min(edgelen, mustmatch);
	      mustmatch = max(dist_to_edge_end, 0);
	    } else {
	      // Try to walk the edge, the first char definitely matches
	      qry_match_len++;
	      refpos = node_start + 1;
	    }
	  }

	  c = GETQCHAR(qrystart + qry_match_len);

	  while (refpos <= node_end && c != '\0') {
	    char r = GETRCHAR(refpos);

	    XPRINTF(" Edge cmp ref: %d %c, qry: %d %c\n", refpos, r, qry_match_len, c);

	    if (r != c) {
	      // mismatch on edge
	      XPRINTF("mismatch on edge: %d, edge_pos: %d\n", qry_match_len, refpos - (node_start));
	      goto RECORD_RESULT;
	    }

	    qry_match_len++;
	    refpos++;

	    c = GETQCHAR(qrystart + qry_match_len);
	  }
	} // while(c)

	XPRINTF("end of string\n");

      RECORD_RESULT: {
	  //PixelOfNode node;
	  //node.data = getnodehist(cur, false);
	  SET_RESULT(cur, result, refpos - node_start, qry_match_len,
		     min_match_len, FORWARD);

	  mustmatch = refpos - node_start;
	  qry_match_len -= mustmatch + 1;
	}
      NEXT_SUBSTRING: {
	  PixelOfNode node;
	  node = GETNODEHIST(prev, false);
	  cur = arrayToAddress(node.suffix);
	}
	//XPRINTF(" following suffix link. mustmatch:%d qry_match_len:%d sl:("fNID")\n",
	//       mustmatch, qry_match_len, NID(cur));
      } //for(int qrystart = 0; qrystart <= last; qrystart++, result += RESULT_SPAN)
    } //for (qryid = 0; qryid < numQueries; ++qryid)
  }

  return 0;
}




//! @todo this is the target function for OpenACC
void computeGold(MatchResults* results,
                 char* refstr,
                 char* queries,
                 int* queryAddrs,
                 int* queryLengths,
                 PixelOfNode* nodeTexture,
                 PixelOfChildren* childrenTexture,
                 int numQueries,
                 int match_length,
                 int rc) {

  if (rc == REVERSE) {

  } else {

      kernel_gold(results,
                  queries,
                  nodeTexture,
                  childrenTexture,
                  refstr,
                  queryAddrs,
                  queryLengths,
                  numQueries,
                  match_length
#if TREE_ACCESS_HISTOGRAM
                  ,int* node_hist,
                  int* child_hist
#endif
                 );
  }
}

