#ifndef FILE_SOUNDPATCH_SEEN
#define FILE_SOUNDPATCH_SEEN

#include <vector>
#include <Stk.h>
#include "ChunkMatch.h"
#include "ChunkGroupMatch.h"
#include "Song.h"

class SoundPatch {
 public :
  SoundPatch() {}
  SoundPatch(ChunkMatch& match, Song& src);
  SoundPatch(ChunkGroupMatch& match, Song& src);
  ~SoundPatch() {};

  int get_length() { return end_time - start_time + 1; }
  
  //time_i is time index in song output, 
  void tick_out(int time_i, stk::StkFrames& out);

  static bool comp_active_start(ChunkMatch &a, ChunkMatch &b);
  static bool comp_active_end(ChunkMatch &a, ChunkMatch &b);
  
 private :
  int start_time;
  int end_time;
  std::vector<stk::StkFrames> sav;

  //i is a direct index to the stkframes
  void get_frame(int i, stk::StkFrames& out);
  
};

#endif
