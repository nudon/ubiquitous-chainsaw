#include "Chunk.h"
#include <iostream>

void Chunk::make_chunk_filter() {
  chunk_filter = ChunkFilter(get_freq_center(), get_freq_margin(),get_bin_size());
}


ChunkFilter Chunk::get_filter() {
  return chunk_filter;
}

bool Chunk:: operator == (Chunk other) {
  return get_chunk_id() == other.get_chunk_id();
}
