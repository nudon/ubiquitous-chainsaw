#ifndef FILE_CHUNKGROUPMATCH_SEEN
#define FILE_CHUNKGROUPMATCH_SEEN

//will also have a lot more funcs in here to use during output phase
#include <list>
#include "ChunkGroup.h"
#include "ChunkGroupCompare.h"

class ChunkGroupMatch {
 public :
  ChunkGroupMatch() {}
  ChunkGroupMatch(ChunkGroup& to_match);
  ~ChunkGroupMatch() {}
  
  double get_score() { return score; }

  ChunkGroup get_orig() { return orig; }
  ChunkGroup get_match() { return match; }

  static std::list<ChunkGroupMatch> best_match(std::list<ChunkGroup>& rec,
					       std::list<ChunkGroup>& rep,
					       ChunkGroupCompare& comp);

  void best_group_match(std::list<ChunkGroup> &many_chunks, ChunkGroupCompare &comp);

  int get_active_start();
  int get_active_end();
  //header for sorting things based on member fields if needed
  //static bool comp_fields(ChunkGroupMatch &a, ChunkGroupMatch &b);
  
 private:
  ChunkGroup orig;
  ChunkGroup match;
  double score;


  
  void set_match(ChunkGroup&  m)  { match = m; }
  void set_score(double s) { score = s; }
};

#endif
