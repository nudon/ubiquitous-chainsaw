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
  SongEmbedder embeder (a, b);
  embeder.funk();
  return 0;
}
