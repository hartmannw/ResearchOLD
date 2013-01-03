#include "MultiBestPath.h"
#include "SpeechFeatures.h"
#include<string>

int main()
{
  std::string fname;
  fileutilities::SpeechFeatures sf 
  fname = "../../babel/features/BP_101/scripted/training/mog_pgram/BABEL_BP_101_11422_20111019_143344_D2_scripted.pgram";

  sf.ReadHtkFile(fname);
  return 0;
}
