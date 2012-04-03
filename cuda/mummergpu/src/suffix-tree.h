#ifndef SUFFIX_TREE_H
#define SUFFIX_TREE_H

#include "mummergpu.h"
#include "suffix-tree.hh"
#include "cpucommon.h"

extern "C"
inline int lookupNumLeaves(ReferencePage * page, TextureAddress addr) {
  unsigned int nodeid = addr2id(addr);
  TextureAddress printParent = NODE_PRINTPARENT(nodeid);
  nodeid = addr2id(printParent);
  return NODE_NUMLEAVES(nodeid);
}

#endif
