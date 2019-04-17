#include "Chunk.h"
#include <iostream>

void Chunk::make_chunk_filter() {
  chunk_filter = ChunkFilter(get_freq_center(), get_freq_margin(),get_bin_size());
}

int Chunk::get_freq_margin() {
  int diff = max_freq - min_freq;
  if (diff % 2 != 0) {
    diff -= 1;
  }
  return diff / 2;
}

int Chunk::get_time_margin() {
  int diff = max_time - min_time;
  if (diff % 2 != 0) {
    diff -= 1;
  }
  return diff / 2;
}


ChunkFilter Chunk::get_filter() {
  return chunk_filter;
}

bool Chunk:: operator == (Chunk other) {
  return get_chunk_id() == other.get_chunk_id();
}
