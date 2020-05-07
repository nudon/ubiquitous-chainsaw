#include "ChunkGroup.h"

using std::min;
using std::max;
using std::list;

static void get_group_stats(list<Chunk>& pool, int* start_time, int* end_time, int* start_freq, int* end_freq);

ChunkGroup::ChunkGroup(list<Chunk>& chunks, double rat) {
  int min_time = -1, max_time = -1, min_freq = -1, max_freq = -1;
  members = std::list<Chunk>(chunks);
  ratio = rat;
  get_group_stats(chunks, &min_time, &max_time, &min_freq, &max_freq);
  super_chunk = Chunk(min_freq, max_freq, min_time, max_time, -1,-1,-1); 
}

ChunkGroup::ChunkGroup() {
  ratio = -1;
}

void get_group_stats(list<Chunk>& pool, int* start_time, int* end_time, int* start_freq, int* end_freq) {
  int start_t = -1, end_t = -1, start_f = -1, end_f = -1;
  bool start = true;
  for (Chunk& c : pool) {
    if (start) {
      start_t = c.get_time_start();
      start_f = c.get_freq_start();
      end_t = c.get_time_end();  
      end_f = c.get_freq_end();
      start = false;
    }
    else {
      start_t = min(start_t, c.get_time_start());
      start_f = min(start_f, c.get_freq_start());
      end_t = max(end_t, c.get_time_end());
      end_f = max(end_f, c.get_freq_end());
    }
  }
  if (start_time != NULL) {
    *start_time = start_t;
  }
  if (start_freq != NULL) {
    *start_freq = start_f;
  }
  if (end_time != NULL) {
    *end_time = end_t;
  }
  if (end_freq != NULL) {
    *end_freq = end_f;
  }
}
