#include "list"
#include "Song.h"
#include "SongEmbedder.h"
using std::string;
using std::list;


int main() {
  //parse args
  std::string song1 = "./sound/tsu.wav";
  std::string song2 = "./sound/tsu.wav";
  Song a (song1);
  Song b (song2);
  SongEmbedder embeder (a, b);
  embeder.funk();
}
