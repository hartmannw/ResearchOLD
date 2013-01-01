#include "SpeechFeatures.h"

int main(int argc, char *argv[])
{
  fileutilities::SpeechFeatures sf;
  sf.ReadCepFile(argv[1]);
  sf.WritePfileAscii(argv[2]);
  return 0;
}
