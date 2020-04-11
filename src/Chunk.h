#ifndef FILE_CHUNK_SEEN
#define FILE_CHUNK_SEEN
#include <LentPitShift.h>
#include "ChunkFilter.h"

//#include <math.h>
class Chunk {
 public:
  Chunk() {
    min_freq = -1;
    max_freq = -1;
    min_time = -1;
    max_time = -1;
    bin_size = -1;
    chunk_id = -1;
    
  }
  Chunk(int minf, int maxf, int mint, int maxt, int bins, int binl, int cid) {
    min_freq = minf;
    max_freq = maxf;
    min_time = mint;
    max_time = maxt;
    bin_size = bins;
    chunk_id = cid;
  }
  
  //~Chunk() { if (chunk_filter != NULL) delete(chunk_filter);}
  ~Chunk() {}

  
  float get_freq_center() { return (max_freq + min_freq) / 2.0; }
  float get_time_center() { return (max_time + min_time) / 2.0; }

  float get_freq_margin() { return (max_freq - min_freq) / 2.0; }
  float get_time_margin() { return (max_time - min_time) / 2.0; }

  int get_time_length() { return max_time - min_time + 1; }
  int get_freq_length() { return max_freq - min_freq + 1; }
  
  int get_time_start() { return std::round(get_time_center() - get_time_margin()); }
  int get_time_end() { return std::round(get_time_center() + get_time_margin()); }

  double get_rel_freq_center() { return get_freq_center() / bin_size; }
  double get_rel_freq_margin() { return get_freq_margin() / bin_size; }

  int get_bin_size() { return bin_size; }
  int get_chunk_id() const { return chunk_id; }
  int get_chunk_size() { return 4 *  get_time_margin() * get_freq_margin(); }

  void set_chunk_id(int nid) { chunk_id = nid; }

  static bool comp_time_center(Chunk& a, Chunk& b) { return a.get_time_center() < b.get_time_center(); }

  static bool  comp_chunk_id(Chunk& a, Chunk& b) { return a.get_chunk_id() < b.get_chunk_id(); }

  static bool  comp_chunk_freq(Chunk& a, Chunk& b) { return a.get_freq_center() < b.get_freq_center(); }

  void make_chunk_filter();
  
  ChunkFilter get_filter() { return chunk_filter; }

  //void make_pit_shift() { pit = stk::LentPitShift(1, bin_size); }
  //stk::LentPitShift &getPitShift() { return pit; }

  bool operator == (Chunk other) const;
  
 private:
  int chunk_id;
  int min_freq;
  int max_freq;
  int min_time;
  int max_time;
  int bin_size;
  ChunkFilter chunk_filter;
  //stk::LentPitShift pit;
  
};



#endif
