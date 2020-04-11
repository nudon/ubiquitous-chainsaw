#include "list"
#include "Song.h"
#include "SongEmbedder.h"
using std::string;
using std::list;


int main(int argc, char* args[]) {
  //parse args
  std::string song1 = "./sound/in1.wav";
  std::string song2 = "./sound/in2.wav";
  if (argc == 3) {
    song1 = args[1];
    song2 = args[2];
  }
  Song a (song1);
  Song b (song2);
  int samples_per_slice = 1024 * 2;
  int edge_snazr = 2;
  int chunk_snazr = 1;
  int match_snazr = 1;

  a.test_chunk_grouping(samples_per_slice, edge_snazr, chunk_snazr);
  //a.make_chunks(samples_per_slice, edge_snazr, chunk_snazr);

  return 0;
  
  SongEmbedder embeder (a, b, samples_per_slice, edge_snazr, chunk_snazr, match_snazr);
  embeder.funk();

  
  
  return 0;
}
