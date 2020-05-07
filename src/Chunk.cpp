#include "Chunk.h"

void Chunk::make_chunk_filter() {
  chunk_filter = ChunkFilter(get_freq_center(), get_freq_margin(),get_bin_size());
}

bool Chunk:: operator == (Chunk other) const{
  return get_chunk_id() == other.get_chunk_id();
}
