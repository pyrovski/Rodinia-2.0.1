#ifndef SUFFIX_TREE_HH
#define SUFFIX_TREE_HH
#include "common.cu"


inline int addr2id(TextureAddress addr) {
#if MERGETEX && REORDER_TREE
  // shift x'a 12th bit as y's 13th
  addr.y |= (addr.x & 0x800) << 1;
  addr.x &= 0x7FF;

  int blocky = addr.y & 0x1F;
  int bigy = addr.y >> 5;
  int bigx = (addr.x << 5) + blocky;
  return bigx + (bigy << 16);

#elif REORDER_TREE
  int blocky = addr.y & 0x1F;
  int bigy = addr.y >> 5;
  int bigx = (addr.x << 5) + blocky;
  return bigx + (bigy << 17);

#elif MERGETEX
  return addr.x;

#else
  return addr.x;

#endif
}

#endif
