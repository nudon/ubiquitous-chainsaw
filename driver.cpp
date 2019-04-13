#include "list"
#include "Song.h"
#include "SongEmbedder.h"
using std::string;
using std::list;


int main() {
  //parse args
  std::string song1 = "./sound/in1.wav";
  std::string song2 = "./sound/in2.wav";
  Song a (song1);
  Song b (song2);
  SongEmbedder embeder (a, b);
  embeder.funk();
  return 0;
}
