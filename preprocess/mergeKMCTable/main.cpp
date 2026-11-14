#include <string>
#include <deque>
#include <iostream>
#include <fstream>
#include "help.hpp"

void helper(void)
{
  std::cout << "This is the function to merge huge KMC files into KMC table!\n"
            << "call by ~ <file_list> (<outname>)\n"
            << "     <file_list> is a pure text where each line is one filename from KMC\n"
            << "     <outname> is set to be final.tsv by default\n";
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 3) { helper(); return 0; }
  std::deque<std::string> que;
  std::ifstream fin(argv[1]);
  if(!fin.is_open()){ helper(); std::cout << argv[1] << " not exits!\n"; return 0; }
  std::string line;
  while(getline(fin, line)) que.push_back(line);
  fin.close();
  size_t tmp_id = 0;
  while(que.size() > 1)
  {
    std::string fname1, fname2;
    fname1 = que.front(); que.pop_front();
    fname2 = que.front(); que.pop_front();
    std::string fname3 = "tmp_" + std::to_string(tmp_id); tmp_id++;
    que.push_back(fname3);
    if(fname1.find("txt") != std::string::npos && fname2.find("txt") != std::string::npos)
    {
      // two KMC table
      std::string id1 = fname1.substr(0, fname1.find('.'));
      std::string id2 = fname2.substr(0, fname2.find('.'));
      StringVector kmer1, kmer2, table1, table2;
      loadKMC(fname1, id1, kmer1, table1);
      loadKMC(fname2, id2, kmer2, table2);
      merge(kmer1, table1, kmer2, table2, fname3);
    }
    else if(fname1.find("txt") != std::string::npos || fname2.find("txt") != std::string::npos)
    {
      // KMC table with sparseTable
      if(fname2.find("txt") != std::string::npos) std::swap(fname1, fname2); // fname2 is KMC table
      std::string id1 = fname1.substr(0, fname1.find('.'));
      StringVector kmer1, table1;
      loadKMC(fname1, id1, kmer1, table1);
      merge(fname2, kmer1, table1, fname3);
    }
    else
    {
      // two sparseTable
      merge(fname1, fname2, fname3);
    }
  }
  line = que.front();
  if(argc == 2) dump(line);
  else          dump(line, argv[2]);
  return 0;
}
