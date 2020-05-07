#ifndef FILE_CHUNKGROUP_SEEN
#define FILE_CHUNKGROUP_SEEN

#include <list>
#include "Chunk.h"


class ChunkGroup {
 public:
  ChunkGroup();
  
  
  ChunkGroup(std::list<Chunk>& chunks, double rat);
  
  ~ChunkGroup() {}
  
  double get_ratio() { return ratio; }
  float get_base_freq() { return members.front().get_freq_center(); }
  double get_rel_base_freq() { return members.front().get_rel_freq_center(); }
  int get_size() { return members.size(); }

  Chunk get_super_chunk() { return super_chunk; }
  std::list<Chunk> get_chunks() { return members; }
 private:
  std::list<Chunk> members;
  double ratio;
  Chunk super_chunk;
  //stkFrames* filtered_group
  
};

#endif
