#ifndef FILE_CHUNK_SEEN
#define FILE_CHUNK_SEEN
#include "ChunkFilter.h"

class Chunk {
 public:
  Chunk() {
    min_freq = -1;
    max_freq = -1;
    min_time = -1;
    max_time = -1;
    bin_size = -1;
    bin_len = -1;
    chunk_id = -1;
    
  }
  Chunk(int minf, int maxf, int mint, int maxt, int bins, int binl, int cid) {
    min_freq = minf;
    max_freq = maxf;
    min_time = mint;
    max_time = maxt;
    bin_size = bins;
    bin_len = binl;
    chunk_id = cid;
  }
  
  //~Chunk() { if (chunk_filter != NULL) delete(chunk_filter);}
  ~Chunk() {}

  
  float get_freq_center() { return (max_freq + min_freq) / 2.0; }
  float get_time_center() { return (max_time + min_time) / 2.0; }

  //wait my margin is just diff in min and max? though it was dist from center to edge
  //my active time might be wrong then, or wherever I start mixing things in matchChunks. would be using twice what it should
  int get_freq_margin();
  int get_time_margin();

  double get_rel_freq_center() { return get_freq_center() / bin_size; }
  double get_rel_freq_margin() { return get_freq_margin() / bin_size; }

  int get_bin_size() { return bin_size; }
  int get_bin_len() { return bin_len; }
  int get_chunk_id() { return chunk_id; }
  int get_chunk_size() { return 4 *  get_time_margin() * get_freq_margin(); }

  

  void make_chunk_filter();
  
  ChunkFilter get_filter();

  bool operator == (Chunk other);
  
 private:
  int chunk_id;
  int min_freq;
  int max_freq;
  int min_time;
  int max_time;
  int bin_size;
  int bin_len;
  ChunkFilter chunk_filter;
  
};



#endif
