#ifndef CPUCOMMON_H
#define CPUCOMMON_H

#define NODE_LENGTH(x)      (page->ref.aux_data[x].length)
#define NODE_PRINTPARENT(x) (page->ref.aux_data[x].printParent)
#define NODE_NUMLEAVES(x)   (page->ref.aux_data[x].numleaves)

#if REORDER_TREE
#define GETNODE(node_addr)     (((PixelOfNode*)    (page->ref.h_node_tex_array))     + (node_addr.x) + (node_addr.y * MAX_TEXTURE_DIMENSION))
#define GETCHILDREN(node_addr) (((PixelOfChildren*)(page->ref.h_children_tex_array)) + (node_addr.x) + (node_addr.y * MAX_TEXTURE_DIMENSION))
#define PADDR(node_addr)        node_addr.x << "," << node_addr.y
#else
#define GETCHILDREN(node_addr) (((PixelOfChildren*)(page->ref.h_children_tex_array)) + (node_addr.x))
#define GETNODE(node_addr)     (((PixelOfNode*)    (page->ref.h_node_tex_array))     + (node_addr.x))
#define PADDR(node_addr)       node_addr.x
#endif

#endif
