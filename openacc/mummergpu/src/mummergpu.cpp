// Includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <sys/time.h>

#include "vector_types.h"

// includes, kernels
#include "common.h"
#include "suffix-tree.h"

#include "mummergpu.h"

int USE_PRINT_KERNEL = 1;

#define BREATHING_ROOM (16 * 1024 * 1024)
#define BASES_PER_TREE_PAGE 8388608
//#define BASES_PER_TREE_PAGE 7000000
#define BLOCKSIZE 256
unsigned int cuda_calls = 0;
void trap_dbg() {
  fprintf(stderr, "Trapped\n");
}

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void runTest( int argc, char** argv);

extern "C"
void computeGold(MatchResults* results,
                 char* refstr,
                 char* queries,
                 int* queryAddrs,
                 int* queryLengths,
                 PixelOfNode* nodeTexture,
                 PixelOfChildren* childrenTexture,
                 int numQueries,
                 int mismatch_length,
                 int rc);

extern "C"
void getReferenceString(const char * filename, char** refstr, size_t* reflen);

extern "C"
void createTreeTexture(const char * filename,
                       PixelOfNode** nodeTexture,
                       PixelOfChildren** childrenTexture,
                       unsigned int* width,
                       unsigned int* node_height,
                       unsigned int* children_height,
                       AuxiliaryNodeData** aux_data,
                       int* num_match_coords,
                       int min_match_len,
                       Statistics* statistics,
                       const char * dotfilename,
                       const char * texfilename);

extern "C"
void getQueriesTexture(int qfile,
                       char** queryTexture,
                       size_t* queryLength,
                       int** queryAddrs,
                       char*** queryNames,
                       int** queryLengths,
                       unsigned int* numQueries,
                       unsigned int* num_match_coords,
                       unsigned int device_memory_avail,
                       int min_match_length,
                       bool rc);

/*
extern "C"
int lookupNumLeaves(ReferencePage * page, TextureAddress addr);
*/

void printAlignments(ReferencePage* page,
                     Alignment* alignments,
                     char* query,
                     int qrylen,
                     TextureAddress nodeid,
                     int qrypos,
                     int edge_depth,
                     int min_match,
                     bool rc,
                     bool forwardcoordinates);

int  countLeafNodes(int nodeid);

extern "C"
void mapQueriesEndToEnd(MatchContext* ctx,
                        ReferencePage* page,
                        MatchInfo* h_matches,
                        unsigned int numMatches,
                        Alignment* h_alignments,
                        unsigned int numAligments);

char *  createTimer() {
  unsigned int * ptr = (unsigned int *) malloc(sizeof(struct Timer_t));
  memset(ptr, 0, sizeof(struct Timer_t));
  return (char *) ptr;
}

void startTimer(char * ptr) {
  gettimeofday(&(((struct Timer_t *)ptr)->start_m), NULL);
}

void stopTimer(char * ptr) {
  gettimeofday(&(((struct Timer_t *)ptr)->end_m), NULL);
}

float getTimerValue(char * ptr) {
  Timer_t * timer = (Timer_t*) ptr;

  if (timer == NULL) {
    fprintf(stderr, "Uninitialized timer!!!\n");
    return 0.0;
  }

  if (timer->end_m.tv_sec == 0) {
    stopTimer(ptr);
  }

  return  (float) (1000.0 * (timer->end_m.tv_sec - timer->start_m.tv_sec)
                   + (0.001 *  (timer->end_m.tv_usec - timer->start_m.tv_usec)));
}

void deleteTimer(char * ptr) {
  free((Timer_t *)ptr);
}

extern "C"
int createReference(const char* fromFile, Reference* ref) {
  if (!fromFile || !ref)
    return -1;

  char * loadreftimer = createTimer();
  startTimer(loadreftimer);

  getReferenceString(fromFile, &(ref->str), &(ref->len));

  stopTimer(loadreftimer);
  ref->t_load_from_disk += getTimerValue(loadreftimer);
  deleteTimer(loadreftimer);

  return 0;
}

extern "C"
int destroyReference(Reference* ref) {
  free(ref->h_node_tex_array);
  free(ref->h_children_tex_array);
  free(ref->str);
#if REORDER_REF
  free(ref->h_ref_array);
#endif

  free(ref->aux_data);
#if TREE_ACCESS_HISTOGRAM
  free(ref->h_node_hist);
  free(ref->h_child_hist);
#endif
  ref->str = NULL;
  ref->len = 0;

  return 0;
}

extern "C"
int createQuerySet(const char* fromFile, QuerySet* queries) {

  fprintf(stderr, "Opening %s...\n", fromFile);
  int qfile = open(fromFile, O_RDONLY);

  if (qfile == -1) {
    fprintf(stderr, "Can't open %s: %d\n", fromFile, errno);
    exit (1);
  }

  queries->qfile = qfile;

  return 0;
}

extern "C"
int destroyQuerySet(QuerySet* queries) {

  if (queries->qfile)
    close(queries->qfile);

  return 0;
}

extern "C"
void printStringForError(int err) {

}

extern "C"
int createMatchContext(Reference* ref,
                       QuerySet* queries,
                       MatchResults* matches,
                       bool on_cpu,
                       int min_match_length,
                       char* stats_file,
                       bool reverse,
                       bool forwardreverse,
                       bool forwardcoordinates,
                       bool showQueryLength,
                       char* dotfilename,
                       char* texfilename,
                       MatchContext* ctx) {

  ctx->queries = queries;
  ctx->ref = ref;
  ctx->full_ref = ref->str;
  ctx->full_ref_len = ref->len;

  ctx->on_cpu = on_cpu;
  ctx->min_match_length = min_match_length;
  ctx->stats_file = stats_file;
  ctx->reverse = reverse;
  ctx->forwardreverse = forwardreverse;
  ctx->forwardcoordinates = forwardcoordinates;
  ctx->show_query_length = showQueryLength;
  ctx->dotfilename = dotfilename;
  ctx->texfilename = texfilename;
  return 0;
}


extern "C"
int destroyMatchContext(MatchContext* ctx) {
  free(ctx->full_ref);
  //destroyReference(ctx->ref);
  destroyQuerySet(ctx->queries);
  return 0;
}

void buildReferenceTexture(Reference* ref,
                           char* full_ref,
                           size_t begin,
                           size_t end,
                           int min_match_len,
                           char* dotfilename,
                           char* texfilename,
                           Statistics* statistics) {
  fprintf(stderr, "Building reference texture...\n");

  PixelOfNode* nodeTexture = NULL;
  PixelOfChildren * childrenTexture = NULL;

  unsigned int width = 0;
  unsigned int node_height = 0;
  unsigned int children_height = 0;

  AuxiliaryNodeData* aux_data = NULL;
  int num_nodes;

  char * loadreftimer = createTimer();
  startTimer(loadreftimer);

  ref->len = end - begin + 3;
  ref->str = (char*)malloc(ref->len);
  ref->str[0] = 's';
  strncpy(ref->str + 1, full_ref + begin, ref->len - 3);
  strcpy(ref->str + ref->len - 2, "$");

  stopTimer(loadreftimer);
  statistics->t_ref_from_disk += getTimerValue(loadreftimer) + ref->t_load_from_disk;
  deleteTimer(loadreftimer);

  createTreeTexture(ref->str,
                    &nodeTexture,
                    &childrenTexture,
                    &width,
                    &node_height,
                    &children_height,
                    &aux_data,
                    &num_nodes,
                    min_match_len,
                    statistics,
                    dotfilename,
                    texfilename);

  ref->h_node_tex_array = nodeTexture;
  ref->h_children_tex_array = childrenTexture;
  ref->tex_width = width;
  ref->tex_node_height = node_height;
  ref->tex_children_height = children_height;

#if TREE_ACCESS_HISTOGRAM
  ref->h_node_hist = (int*)calloc(width * node_height, sizeof(int));
  ref->h_child_hist = (int*)calloc(width * children_height, sizeof(int));
#endif

  ref->aux_data = aux_data;
  ref->num_nodes = num_nodes;

  ref->bytes_on_board = (width * node_height * sizeof(PixelOfNode)) +
                        (width * children_height * sizeof(PixelOfChildren));
  fprintf(stderr, "This tree will need %ld bytes on the board\n", ref->bytes_on_board);

#if REORDER_REF
  char * reordertimer = createTimer();
  startTimer(reordertimer);

  unsigned int refpitch = ref->pitch = 65536;
  int numrows = ceil(ref->len / ((float)refpitch));
  int blocksize = 4;
  numrows += blocksize;

  int refstrsize = numrows * refpitch;
  ref->h_ref_array = (char *) malloc(refstrsize);
  ref->bytes_on_board += refstrsize;

  fprintf(stderr, "The refstr (reordered) requires %d bytes\n", refstrsize);

  int z_max = numrows * refpitch;
  for (int z = 0; z < z_max; z++) {
    ref->h_ref_array[z] = 'Z';
  }

  int x, y;
  int maxx = 0, maxy = 0;

  size_t reflen = ref->len;
  char* refstr = ref->str;


  int block_dim = refpitch * blocksize;
  for (int i = 0; i < reflen; i++) {
    int bigx = i % (block_dim); // ref string reorder
    int bigy = i / (block_dim);

    y = bigy * blocksize + bigx % blocksize;
    x = bigx / blocksize;

    //   printf("%d: (%d,%d)=%c\n", i, x, y, refstr[i]);

    assert(x < refpitch);
    assert(y < numrows);

    ref->h_ref_array[y*refpitch+x] = refstr[i];

    if (x > maxx) {
      maxx = x;
    }
    if (y > maxy) {
      maxy = y;
    }
  }

  if ((maxx >= refpitch) || (maxy >= numrows)) {
    fprintf(stderr, "ERROR: maxx: %d refpitch: %d, maxy: %d numrows: %d\n",
            maxx,    refpitch,     maxy,    numrows);

    exit(1);
  }
  stopTimer(reordertimer);
  if (statistics)
    statistics->t_reorder_ref_str += getTimerValue(reordertimer);
  deleteTimer(reordertimer);
#else
  fprintf(stderr, "The refstr requires %ld bytes\n", ref->len);
  ref->bytes_on_board += ref->len;
#endif


}

void loadReferenceTexture(MatchContext* ctx) {
  Reference* ref = ctx->ref;
  int numrows = ceil(ref->len / ((float)ref->pitch));
  int blocksize = 4;
  numrows += blocksize;

  cudaChannelFormatDesc refTextureDesc =
    cudaCreateChannelDesc(8, 0, 0, 0, cudaChannelFormatKindSigned);

  if (!ctx->on_cpu) {

  } else {
    ref->d_ref_array = NULL;
  }
}


void unloadReferenceString(Reference* ref) {
  ref->d_ref_array = NULL;
}


//loads a tree and text for [begin, end) in the reference
void loadReference(MatchContext* ctx) {

  Reference* ref = ctx->ref;

  ref->bytes_on_board = 0;

  loadReferenceTexture(ctx);

  if (!ctx->on_cpu) {

  } else {
    ref->d_node_tex_array = NULL;
    ref->d_children_tex_array = NULL;
  }
}



void dumpQueryBlockInfo(QuerySet* queries) {
  fprintf(stderr, "\tProcessing queries %s to %s\n",
          queries->h_names[0],
          queries->h_names[queries->count-1]);
}

void loadQueries(MatchContext* ctx) {
  QuerySet* queries = ctx->queries;
  queries->bytes_on_board = 0;

  unsigned int numQueries = queries->count;

  if (!ctx->on_cpu) {

  } else {
    queries->d_addrs_tex_array = NULL;
    queries->d_tex_array = NULL;
    queries->d_lengths_array = NULL;
    fprintf(stderr, " allocated %ld bytes\n", 2 * numQueries*sizeof(int) + queries->texlen);
  }


}


// Computes the location of the first MatchCoord for a given query.  NOTE:
// Do NOT use this function if COALESCED_QUERIES == 1
inline int match_coord_addrs(int qryid, int qry_addrs, int match_length) {
  return qry_addrs - qryid * (match_length + 1);
}

// Construct the offset table for a set of queries.  This table will be used
// by the printing functions, and if COALESCED_QUERIES == 1, by the matching
// kernel.
void buildCoordOffsetArray(MatchContext* ctx,
                           int** h_coord_offset_array,
                           unsigned int* num_coords) {
  int numCoords = 0;
  int match_length = ctx->min_match_length;
  int numQueries = ctx->queries->count;
  int* lengths = ctx->queries->h_lengths_array;

  int* coord_offsets = (int*)calloc(numQueries, sizeof(int));

#if COALESCED_QUERIES

  for (unsigned int i = 0; i < numQueries; i += WARP_SIZE) {
    // Every query in this warp will need at least this many coords
    int max_num_coords = 0;
    for (unsigned int j = 0; j < WARP_SIZE && (i + j) < numQueries; ++j) {
      int num_coords = lengths[i + j] - match_length + 1;
      if ( max_num_coords <  num_coords)
        max_num_coords = num_coords;
    }

    unsigned int block_size = max_num_coords * WARP_SIZE;

    for (unsigned int j = 0; j < WARP_SIZE && (i + j) < numQueries; ++j) {
      ctx->results.h_coord_tex_array[i + j] = numCoords + j;
    }
    numCoords += block_size;
  }
#else
  for (unsigned int i = 0; i < numQueries; ++i) {
    int qryoffset = ctx->queries->h_addrs_tex_array[i];
    coord_offsets[i] = match_coord_addrs(i, qryoffset, match_length);
  }
  if (numQueries > 0) {
    unsigned int last_qry = numQueries - 1;
    unsigned int last_qry_len = lengths[last_qry] - match_length + 1;
    numCoords = coord_offsets[last_qry] + last_qry_len;
    fprintf(stderr, "Need %d match coords for this result array\n",
            numCoords);
  }
#endif
  *num_coords = numCoords;
  *h_coord_offset_array = coord_offsets;
}


void loadResultBuffer(MatchContext* ctx) {
  unsigned int numQueries = ctx->queries->count;

  assert (numQueries);

  char* offsettimer = createTimer();
  startTimer(offsettimer);

  buildCoordOffsetArray(ctx,
                        &(ctx->results.h_coord_tex_array),
                        &(ctx->results.numCoords));

  stopTimer(offsettimer);
  ctx->statistics.t_build_coord_offsets += getTimerValue(offsettimer);
  deleteTimer(offsettimer);

  unsigned int numCoords = ctx->results.numCoords;
  fprintf(stderr, "Allocating result array for %d queries (%ld bytes) ...",
          numQueries, numCoords*sizeof(MatchCoord) );

  size_t boardFreeMemory = 0;
  size_t total_mem = 0;

  fprintf(stderr,"board free memory: %ld total memory: %ld\n",
          boardFreeMemory, total_mem);

  ctx->results.h_match_coords = (MatchCoord*)calloc( numCoords, sizeof(MatchCoord));
  /* pinned memory did not improve runtime
  cudaMallocHost(&ctx->results.h_match_coords, numCoords * sizeof(MatchCoord));
  memset(ctx->results.h_match_coords, 0, numCoords * sizeof(MatchCoord));
  */
  if (ctx->results.h_match_coords == NULL) {
    trap_dbg();
    exit(EXIT_FAILURE);
  }

  if (!ctx->on_cpu) {

  } else {
    ctx->results.d_match_coords = NULL;
  }

  fprintf(stderr, "done\n");
}




int flushOutput();
int addToBuffer(char* string);

char numbuffer[32];

MatchCoord* coordForQueryChar(MatchContext* ctx,
                              unsigned int qryid,
                              unsigned int qrychar) {
  MatchResults* results = &(ctx->results);
  MatchCoord* coords = results->h_match_coords;
#if COALESCED_QUERIES
  return coords + results->h_coord_tex_array[qryid] + qrychar * WARP_SIZE;
#else
  return coords + results->h_coord_tex_array[qryid] + qrychar;
#endif
}

void coordsToPrintBuffers(MatchContext* ctx,
                          ReferencePage* page,
                          MatchInfo** matches,
                          Alignment** alignments,
                          unsigned int mem_avail,
                          unsigned int* coord_idx,
                          unsigned int* match_idx,
                          unsigned int* align_idx,
                          unsigned int* nextqry,
                          unsigned int* nextqrychar) {
  unsigned int numQueries = ctx->queries->count;
  int match_length = ctx->min_match_length;
  unsigned int cidx = *coord_idx;
  unsigned int midx = 0;

  unsigned int numCoords = ctx->results.numCoords;

  unsigned int numMatches = 0;
  unsigned int numAlignments = 0;

  int DEBUG = 0;
  if (DEBUG && cidx == 0) {
    for (int j = 0; j < numCoords; ++j) {
      MatchCoord * coord = ctx->results.h_match_coords+j;
      if (coord->node.data > 0 && !(coord->edge_match_length & FRMASK)) {
        //fprintf(stdout, "node: %d\n",
        //        coord->node);
        fprintf(stdout, "node: %d leaves:%d\n",
                coord->node.data, lookupNumLeaves(page, coord->node));
      }
    }
    exit(0);
  }


  // How much can we fit into mem_avail?
  for (int j = cidx; j < numCoords; ++j) {
    MatchCoord* coord = ctx->results.h_match_coords + j;

    int queryAlignments = 0;
    int queryMatches = 0;

    if (coord->node.data > 0 && !(coord->edge_match_length & FRMASK)) {
      int numLeaves = lookupNumLeaves(page, coord->node);
      queryAlignments += numLeaves;
      queryMatches++;
    }
    int allMatches    = numMatches    + queryMatches;
    int allAlignments = numAlignments + queryAlignments;

    int neededSize = allMatches * sizeof(MatchInfo) + allAlignments * sizeof(Alignment);

    if (neededSize > mem_avail || (allMatches/BLOCKSIZE) >= MAX_GRID_DIMENSION) {
      // adding this match won't fit on the board
      break;
    }

    ++cidx;
    numMatches    = allMatches;
    numAlignments = allAlignments;
  }

  MatchInfo* M = (MatchInfo*)calloc(numMatches, sizeof(MatchInfo));
  unsigned int alignmentOffset = 0;

  int qry = *nextqry;
  int qrychar = *nextqrychar;
  bool set_full = false;
  while (qry < numQueries) {
    // h_lengths_array doesn't count the 'q' at the beginning of each query
    int qlen = ctx->queries->h_lengths_array[qry] + 1 - match_length;

    while (qrychar < qlen) {
      if (midx >= numMatches) {
        set_full = true;
        break;
      }

      MatchCoord* coord = coordForQueryChar(ctx, qry, qrychar);

      if (coord->node.data > 0 && !(coord->edge_match_length & FRMASK)) {
        MatchInfo m;
        m.resultsoffset = alignmentOffset;
        m.qrystartpos = qrychar;
        m.matchnode = coord->node;
        m.edgematch = coord->edge_match_length;
        m.numLeaves = lookupNumLeaves(page, m.matchnode);
        m.queryid = qry;

        alignmentOffset += m.numLeaves;
        M[midx++] = m;
      }

      ++qrychar;
    }

    if (set_full)
      break;

    ++qry;
    qrychar = 0;
  }

  *coord_idx = cidx;
  *match_idx = midx;
  *align_idx = alignmentOffset;
  *matches = M;
  *nextqry = qry;
  *nextqrychar = qrychar;
  fprintf(stderr, "Allocing %ld bytes of host memory for %d alignments\n",
          alignmentOffset * sizeof(Alignment), numAlignments);
  *alignments = (struct Alignment *) calloc(alignmentOffset, sizeof(Alignment));
  //cudaMallocHost((void**)alignments, numAlignments * sizeof(Alignment));
}


// TODO: need reverse-complement printing support
void runPrintOnCPU(MatchContext* ctx, ReferencePage* page,
                   MatchInfo* h_matches,
                   unsigned int numMatches,
                   Alignment* alignments,
                   unsigned int numAlignments) {
  unsigned int min_match_length = ctx->min_match_length;

  int* addrs = ctx->queries->h_addrs_tex_array;
  int* lengths = ctx->queries->h_lengths_array;
  char* qrychars = ctx->queries->h_tex_array;

  if (!numMatches)
    return;

  int qry = -1;
  unsigned int qrylen;

  for (int i = 0; i < numMatches; ++i) {
    MatchInfo& match = h_matches[i];
    if (match.queryid != qry) {
      qry = match.queryid;
      qrylen = lengths[qry];
    }
    if (!(match.edgematch & FRMASK)) {
      printAlignments(page,
                      alignments + match.resultsoffset,
#if COALESCED_QUERIES
                      qrychars + sizeof(int) * addrs[qry],
#else
                      qrychars + addrs[qry],
#endif
                      qrylen,
                      match.matchnode,
                      match.qrystartpos,
                      match.edgematch,
                      min_match_length,
                      0,
                      ctx->forwardcoordinates);
    }
  }
}

int addMatchToBuffer(int left_in_ref, int qrypos, int matchlen);

void getExactAlignments(MatchContext * ctx, ReferencePage * page, bool on_cpu) {
  assert(!ctx->reverse && !ctx->forwardreverse);

  size_t boardFreeMemory;

  if (!on_cpu) {
  } else {
    boardFreeMemory = 256 * 1024 * 1024;
  }

#ifdef __DEVICE_EMULATION__
  boardFreeMemory = 512 * 1024 * 1024;
#endif

  boardFreeMemory -= BREATHING_ROOM;
  fprintf(stderr, "board free memory: %lu\n", boardFreeMemory);

  int rTotalMatches = 0;
  int rTotalAlignments = 0;
  int totalRounds = 0;
  unsigned int last_coord = ctx->results.numCoords;
  unsigned int next_coord = 0;
  unsigned int nextqry = 0;
  unsigned int nextqrychar = 0;
  int lastqry = -1;
  while (next_coord < last_coord) {
    // see how many queries will fit on the board
    totalRounds++;

    unsigned int numMatches = 0;
    unsigned int numAlignments = 0;
    MatchInfo* h_matches = NULL;
    Alignment* h_alignments = NULL;
    int coord_left = next_coord;
    char* btimer = createTimer();
    startTimer(btimer);
    coordsToPrintBuffers(ctx, page, &h_matches, &h_alignments, boardFreeMemory,
                         &next_coord, &numMatches, &numAlignments, &nextqry, &nextqrychar);
    stopTimer(btimer);

    float btime = getTimerValue(btimer);
    ctx->statistics.t_coords_to_buffers += btime;
    fprintf(stderr, "buffer prep time= %f\n", btime);
    deleteTimer(btimer);

    fprintf(stderr, "Round %d: Printing results for match coords [%d-%d) of %d using %d matches and %d alignments\n",
            totalRounds, coord_left, next_coord, last_coord, numMatches, numAlignments);

    if (numMatches == 0)
      continue;

    char buf[256];
    //assert(qryend > qrystart);

    rTotalAlignments += numAlignments;
    rTotalMatches += numMatches;

    char* ktimer = createTimer();
    startTimer(ktimer);
    if (on_cpu) {
      runPrintOnCPU(ctx, page, h_matches, numMatches,
                    h_alignments, numAlignments);
    } else {
    }
    stopTimer(ktimer);

    float ktime = getTimerValue(ktimer);
    ctx->statistics.t_print_kernel += ktime;
    fprintf(stderr, "print kernel time= %f\n", ktime);
    deleteTimer(ktimer);

    // char* stimer = createTimer();
    // startTimer(stimer);
    // mapQueriesEndToEnd(ctx,
    //                    page,
    //                    h_matches,
    //                    numMatches,
    //                    h_alignments,
    // 				   numAlignments);
    //
    // stopTimer(stimer);
    //
    // float stime = getTimerValue(stimer);
    // fprintf(stderr, "postprocess time= %f\n", stime);
    // deleteTimer(stimer);

    //flushOutput();

    //Process the alignments
    char* otimer = createTimer();
    startTimer(otimer);

    for (int m = 0; m < numMatches; m++) {
      int base = h_matches[m].resultsoffset;
      for (int i = 0; i < h_matches[m].numLeaves; i++) {
        // See if there are any more left maximal alignments for this match
        if (h_alignments[base+i].left_in_ref == 0) {
          break;
        }

        if (h_matches[m].queryid != lastqry) {
          lastqry = h_matches[m].queryid;
          addToBuffer("> ");
          addToBuffer(*(ctx->queries->h_names + lastqry));
          addToBuffer("\n");
        }

        sprintf(buf, "%d\t%d\t%d\n",
                h_alignments[base+i].left_in_ref,
                h_matches[m].qrystartpos + 1,
                h_alignments[base+i].matchlen);
        addToBuffer(buf);

        // addMatchToBuffer(h_alignments[base+i].left_in_ref,
        // 								 h_matches[m].qrystartpos + 1,
        // 								h_alignments[base+i].matchlen);

      }
    }


    flushOutput();

    stopTimer(otimer);
    ctx->statistics.t_results_to_disk += getTimerValue(otimer);
    deleteTimer(otimer);

    free(h_matches);
    free(h_alignments);
    //cudaFreeHost((void*)h_alignments);

  }
  free(ctx->results.h_coord_tex_array);
  //cudaFreeHost(ctx->results.h_match_coords);
  free(ctx->results.h_match_coords);
  ctx->results.h_coord_tex_array = NULL;
  ctx->results.h_match_coords = NULL;

  fprintf(stderr, "Finished processing %d matches and %d potential alignments in %d rounds\n",
          rTotalMatches, rTotalAlignments, totalRounds);
}

int getQueryBlock(MatchContext* ctx, size_t device_mem_avail) {
  QuerySet* queries = ctx->queries;
  char * queryTex = NULL;
  int* queryAddrs = NULL;
  int* queryLengths = NULL;
  unsigned int numQueries;
  unsigned int num_match_coords;
  size_t queryLen;
  char** names;

  fprintf(stderr, "Loading query block... ");

  char* queryreadtimer = createTimer();
  startTimer(queryreadtimer);

  getQueriesTexture(queries->qfile,
                    &queryTex,
                    &queryLen,
                    &queryAddrs,
                    &names,
                    &queryLengths,
                    &numQueries,
                    &num_match_coords,
                    device_mem_avail,
                    ctx->min_match_length,
                    ctx->reverse || ctx->forwardreverse);

  stopTimer(queryreadtimer);
  ctx->statistics.t_queries_from_disk += getTimerValue(queryreadtimer);
  deleteTimer(queryreadtimer);

  queries->h_tex_array = queryTex;
  queries->count = numQueries;
  queries->h_addrs_tex_array = queryAddrs;
  queries->texlen = queryLen;
  queries->h_names = names;
  queries->h_lengths_array = queryLengths;

  ctx->results.numCoords = num_match_coords;

  fprintf(stderr, "done.\n");

  return numQueries;
}

void destroyQueryBlock(QuerySet* queries) {
  free(queries->h_tex_array);
  queries->h_tex_array = NULL;

  for (int i = 0; i < queries->count; ++i)
    free(queries->h_names[i]);

  free(queries->h_names);

  queries->count = 0;
  queries->texlen = 0;

  free(queries->h_addrs_tex_array);
  queries->h_addrs_tex_array = NULL;

  free(queries->h_lengths_array);
  queries->h_lengths_array = NULL;
}

void resetStats(Statistics* stats) {
  stats->t_end_to_end = 0.0;
  stats->t_match_kernel = 0.0;
  stats->t_print_kernel = 0.0;
  stats->t_queries_to_board = 0.0;
  stats->t_match_coords_to_board = 0.0;
  stats->t_match_coords_from_board = 0.0;
  stats->t_tree_to_board = 0.0;
  stats->t_ref_str_to_board = 0.0;
  stats->t_queries_from_disk = 0.0;
  stats->t_ref_from_disk = 0.0;
  stats->t_results_to_disk = 0.0;
  stats->t_tree_construction = 0.0;
  stats->t_tree_reorder = 0.0;
  stats->t_tree_flatten = 0.0;
  stats->t_reorder_ref_str = 0.0;
  stats->t_build_coord_offsets = 0.0;
  stats->t_coords_to_buffers = 0.0;
  stats->bp_avg_query_length = 0.0;

#if TREE_ACCESS_HISTOGRAM
  if (stats->node_hist_size) {
    free(stats->node_hist);
    stats->node_hist = NULL;
    stats->node_hist_size = 0;
  }

  if (stats->child_hist_size) {
    free(stats->child_hist);
    stats->child_hist = NULL;
    stats->child_hist_size = 0;
  }
#endif
}

void writeStatisticsFile(Statistics* stats,
                         char* stats_filename,
                         char* node_hist_filename = NULL,
                         char* child_hist_filename = NULL) {
  if (stats_filename) {
    FILE* f = fopen(stats_filename, "w");

    if (!f) {
      fprintf(stderr, "WARNING: could not open %s for writing\n", stats_filename);
    } else {
      fprintf(f, "Q");
      fprintf(f, ",R");
      fprintf(f, ",T");
      fprintf(f, ",m");
      fprintf(f, ",r");
      fprintf(f, ",t");
      fprintf(f, ",n");
      fprintf(f, ",Total");
      fprintf(f, ",Match kernel");
      fprintf(f, ",Print Kernel");
      fprintf(f, ",Queries to board");
      fprintf(f, ",Match coords to board");
      fprintf(f, ",Match coords from board");
      fprintf(f, ",Tree to board");
      fprintf(f, ",Ref str to board");
      fprintf(f, ",Queries from disk");
      fprintf(f, ",Ref from disk");
      fprintf(f, ",Output to disk");
      fprintf(f, ",Tree construction");
      fprintf(f, ",Tree reorder");
      fprintf(f, ",Tree flatten");
      fprintf(f, ",Ref reorder");
      fprintf(f, ",Build coord table");
      fprintf(f, ",Coords to buffers");
      fprintf(f, ",Avg qry length");
      fprintf(f, "\n");

      fprintf(f, "%d", QRYTEX);
      fprintf(f, ",%d", REFTEX);
      fprintf(f, ",%d", TREETEX);
      fprintf(f, ",%d", MERGETEX);
      fprintf(f, ",%d", REORDER_REF);
      fprintf(f, ",%d", REORDER_TREE);
      fprintf(f, ",%d", RENUMBER_TREE);
      fprintf(f, ",%f", stats->t_end_to_end);
      fprintf(f, ",%f", stats->t_match_kernel);
      fprintf(f, ",%f", stats->t_print_kernel);
      fprintf(f, ",%f", stats->t_queries_to_board);
      fprintf(f, ",%f", stats->t_match_coords_to_board);
      fprintf(f, ",%f", stats->t_match_coords_from_board);
      fprintf(f, ",%f", stats->t_tree_to_board);
      fprintf(f, ",%f", stats->t_ref_str_to_board);
      fprintf(f, ",%f", stats->t_queries_from_disk);
      fprintf(f, ",%f", stats->t_ref_from_disk);
      fprintf(f, ",%f", stats->t_results_to_disk);
      fprintf(f, ",%f", stats->t_tree_construction);
      fprintf(f, ",%f", stats->t_tree_reorder);
      fprintf(f, ",%f", stats->t_tree_flatten);
      fprintf(f, ",%f", stats->t_reorder_ref_str);
      fprintf(f, ",%f", stats->t_build_coord_offsets);
      fprintf(f, ",%f", stats->t_coords_to_buffers);
      fprintf(f, ",%f", stats->bp_avg_query_length);
      fprintf(f,"\n");

      fclose(f);
    }
  }
#if TREE_ACCESS_HISTOGRAM
  if (node_hist_filename) {
    FILE* f = fopen(node_hist_filename, "w");
    if (!f) {
      fprintf(stderr, "WARNING: could not open %s for writing\n", node_hist_filename);
    } else {
      for (unsigned int i = 0; i < ctx->statistics.node_hist_size; ++i)
        fprintf(f, "%d\t%d\n", i, ctx->statistics.node_hist[i]);
    }

  }

  if (child_hist_filename) {
    FILE* f = fopen(child_hist_filename, "w");
    if (!f) {
      fprintf(stderr, "WARNING: could not open %s for writing\n", child_hist_filename);
    } else {
      for (unsigned int i = 0; i < ctx->statistics.child_hist_size; ++i)
        fprintf(f, "%d\t%d\n", i, ctx->statistics.child_hist[i]);
    }

  }

  float total_node_hits = 0;
  float tree_top_node_hits = 0;

  float total_child_hits = 0;
  float tree_top_child_hits = 0;

  for (unsigned int i = 0; i < ctx->statistics.node_hist_size; ++i) {
    total_node_hits +=ctx->statistics.node_hist[i];
    if (i < 256) {
      tree_top_node_hits += ctx->statistics.node_hist[i];
    }
  }

  for (unsigned int i = 0; i < ctx->statistics.child_hist_size; ++i) {
    total_child_hits +=ctx->statistics.child_hist[i];
    if (i < 256) {
      tree_top_child_hits += ctx->statistics.child_hist[i];
    }
  }

  fprintf(stderr, "Tree top node  hits (%d/%d) = %f percent\n",(int)tree_top_node_hits, (int)total_node_hits, tree_top_node_hits /total_node_hits);
  fprintf(stderr, "Tree top child hits (%d/%d) = %f percent\n",(int)tree_top_child_hits, (int)total_child_hits, tree_top_child_hits /total_child_hits);
#endif
}

void matchOnCPU(MatchContext* ctx, bool doRC) {
  //TODO: CPU is matching is disabled.
  if (doRC) {
    // Match the reverse complement of the queries to the ref
    computeGold(&ctx->results,
                ctx->ref->str,
                ctx->queries->h_tex_array,
                ctx->queries->h_addrs_tex_array,
                ctx->queries->h_lengths_array,
                (PixelOfNode*)(ctx->ref->h_node_tex_array),
                (PixelOfChildren*)(ctx->ref->h_children_tex_array),
                ctx->queries->count,
                ctx->min_match_length,
                REVERSE);
  } else {
    computeGold(&ctx->results,
                ctx->ref->str,
                ctx->queries->h_tex_array,
                ctx->queries->h_addrs_tex_array,
                ctx->queries->h_lengths_array,
                (PixelOfNode*)(ctx->ref->h_node_tex_array),
                (PixelOfChildren*)(ctx->ref->h_children_tex_array),
                ctx->queries->count,
                ctx->min_match_length,
                FORWARD);
  }
}


void matchQueryBlockToReferencePage(MatchContext* ctx,
                                    ReferencePage* page,
                                    bool reverse_complement) {
  char*  ktimer = createTimer();

  fprintf(stderr, "Memory footprint is:\n\tqueries: %ld\n\tref: %ld\n\tresults: %ld\n",
          ctx->queries->bytes_on_board,
          ctx->ref->bytes_on_board,
          ctx->results.bytes_on_board);

  startTimer(ktimer);
  if (ctx->on_cpu) {
    matchOnCPU(ctx, reverse_complement);
  } else {

  }
  stopTimer(ktimer);

  float ktime = getTimerValue(ktimer);
  ctx->statistics.t_match_kernel += ktime;
  fprintf(stderr, "match kernel time= %f\n", ktime);
  deleteTimer(ktimer);

}


int matchSubset(MatchContext* ctx,
                ReferencePage* page) {

  loadQueries(ctx);

  fprintf(stderr,
          "Matching queries %s - %s against ref coords %d - %d\n",
          ctx->queries->h_names[0],
          ctx->queries->h_names[ctx->queries->count - 1],
          page->begin,
          page->end);

  loadResultBuffer(ctx);

  // TODO: renable RC support by calling this twice /w reverse/fwdreverse
  // idiom.
  matchQueryBlockToReferencePage(ctx, page, false);

  if (USE_PRINT_KERNEL && !ctx->on_cpu) {
  } else {
    getExactAlignments(ctx, page, true);
  }

  flushOutput();
  return 0;
}

int getFreeDeviceMemory(bool on_cpu) {
  size_t free_mem = 0;

  // We have to 'prime' CUDA by making an allocation here.  cuMemGetInfo
  // will return zeroes until we do a malloc.
  if (!on_cpu) {

  } else {
    free_mem = 804585472; // pretend we are on a 8800 GTX
  }

  return free_mem;
}

int matchQueriesToReferencePage(MatchContext* ctx, ReferencePage* page) {
  fprintf(stderr, "Beginning reference page %p\n", page);

  int free_mem = getFreeDeviceMemory(ctx->on_cpu);

  int available_mem = free_mem - page->ref.bytes_on_board - BREATHING_ROOM;
  ctx->ref = &(page->ref);
  loadReference(ctx);

  while (getQueryBlock(ctx, available_mem)) {
    matchSubset(ctx, page);
    ctx->statistics.bp_avg_query_length =
      ctx->queries->texlen / (float)(ctx->queries->count) - 2;
    destroyQueryBlock(ctx->queries);
  }

  unloadReferenceString(ctx->ref);
  lseek(ctx->queries->qfile, 0, SEEK_SET);
  return 0;
}



void initReferencePages( MatchContext* ctx , int* num_pages, ReferencePage** pages_out) {
  unsigned int bases_in_ref = ctx->full_ref_len - 3;
  unsigned int page_size = BASES_PER_TREE_PAGE < bases_in_ref ?
                           BASES_PER_TREE_PAGE : bases_in_ref;
  unsigned int num_reference_pages = ceil((bases_in_ref + 0.0) / page_size);
  fprintf(stderr, "Stream will use %d pages for %d bases, page size = %d\n",
          num_reference_pages, bases_in_ref, page_size);

  unsigned int page_overlap = MAX_QUERY_LEN + 1;
  ReferencePage* pages = (ReferencePage*) calloc(num_reference_pages,
                         sizeof(ReferencePage));

  pages[0].begin = 1;
  pages[0].end = pages[0].begin +
                 page_size  +
                 ceil(page_overlap / 2.0) + 1; //the 1 is for the 's' at the beginning
  pages[0].shadow_left = -1;
  pages[0].id = 0;

  for (int i = 1; i < num_reference_pages - 1; ++i) {
    pages[i].begin = pages[i - 1].end - page_overlap;
    pages[i].end = pages[i].begin + page_size +  page_overlap;
    pages[i - 1].shadow_right = pages[i].begin;
    pages[i].shadow_left = pages[i-1].end;
    pages[i].id = i;
  }

  if (num_reference_pages > 1) {
    int last_page = num_reference_pages - 1;
    pages[last_page].begin = pages[last_page - 1].end - page_overlap;
    pages[last_page].end = ctx->full_ref_len - 1;
    pages[last_page - 1].shadow_right = pages[last_page].begin;
    pages[last_page].shadow_right = -1;
    pages[last_page].shadow_left = pages[last_page - 1].end;
    pages[last_page].id = last_page;
  }

  *pages_out = pages;
  *num_pages = num_reference_pages;
}

int streamReferenceAgainstQueries(MatchContext* ctx) {
  int num_reference_pages = 0;
  ReferencePage* pages = NULL;
  initReferencePages(ctx, &num_reference_pages, &pages);


  buildReferenceTexture(&(pages[0].ref),
                        ctx->full_ref,
                        pages[0].begin,
                        pages[0].end,
                        ctx->min_match_length,
                        ctx->dotfilename,
                        ctx->texfilename,
                        &(ctx->statistics));


  matchQueriesToReferencePage(ctx, &pages[0]);
  destroyReference(&(pages[0].ref));

  for (int i = 1; i < num_reference_pages - 1; ++i) {

    buildReferenceTexture(&(pages[i].ref),
                          ctx->full_ref,
                          pages[i].begin,
                          pages[i].end,
                          ctx->min_match_length,
                          NULL,
                          NULL,
                          &(ctx->statistics));

    matchQueriesToReferencePage(ctx, &pages[i]);
    destroyReference(&(pages[i].ref));
  }

  if (num_reference_pages > 1) {
    int last_page = num_reference_pages - 1;
    buildReferenceTexture(&(pages[last_page].ref),
                          ctx->full_ref,
                          pages[last_page].begin,
                          pages[last_page].end,
                          ctx->min_match_length,
                          NULL,
                          NULL,
                          &(ctx->statistics));

    matchQueriesToReferencePage(ctx, &pages[last_page]);
    destroyReference(&(pages[last_page].ref));
  }
  free(pages);
  return 0;
}


extern "C"
int matchQueries(MatchContext* ctx) {
  assert(sizeof(struct PixelOfNode) == sizeof(uint4));
  assert(sizeof(struct PixelOfChildren) == sizeof(uint4));

#if TREE_ACCESS_HISTOGRAM
  ctx->statistics.node_hist_size = 0;
  ctx->statistics.child_hist_size = 0;
#endif

  resetStats(&(ctx->statistics));

  char* ttimer = createTimer();
  startTimer(ttimer);

  int ret;

  fprintf(stderr, "Streaming reference pages against all queries\n");
  ret = streamReferenceAgainstQueries(ctx);

  stopTimer(ttimer);
  ctx->statistics.t_end_to_end += getTimerValue(ttimer);
  deleteTimer(ttimer);

  writeStatisticsFile(&(ctx->statistics), ctx->stats_file, "node_hist.out", "child_hist.out");

  return ret;
}



