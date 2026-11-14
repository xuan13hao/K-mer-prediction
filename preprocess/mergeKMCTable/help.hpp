#ifndef _MERGE_
#define _MERGE_
#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;
typedef vector<string> StringVector;
int MIN_COUNT = 10;

#ifdef FLAG_VERBOSE
  bool FF_VEBOSE = true;
#else
  bool FF_VEBOSE = false;
#endif

void loadSparseLine(ifstream& fin, string& kmer, string& row)
{
  kmer.clear(); row.clear();
  size_t size;
  fin.read((char*)&size, sizeof(size_t));
  kmer.resize(size);
  fin.read((char*)&kmer[0], size);

  fin.read((char*)&size, sizeof(size_t));
  row.resize(size);
  fin.read((char*)&row[0], size);
}

void loadKMC(const string& fname, const string& id, StringVector& kmers, StringVector& table)
{
  clock_t t;
  ifstream fin(fname);
  string line;
  while(getline(fin, line))
  {
    StringVector cont;
    boost::split( cont, line, boost::is_any_of( " \t" ) );
    if (stoi(cont[1]) >=  MIN_COUNT)
    {
      kmers.push_back(cont[0]);
      table.push_back(id+"\t"+cont[1]);
    }
  }
  fin.close();
  if(FF_VEBOSE) cout << "loadKMC " << fname << " takes " << ((double)clock() - t)/CLOCKS_PER_SEC << " s\n";
}

void loadSparse(const string& fname, StringVector& kmers, StringVector& table)
{
  clock_t t;
  ifstream fin(fname, ios::binary);
  int cnt;
  while(!fin.eof())
  {
    cnt++;
    size_t size; string tmp;
    fin.read((char*)&size, sizeof(size_t));
    tmp.resize(size);
    fin.read((char*)&tmp[0], size);
    kmers.push_back(tmp);
    tmp.clear();

    fin.read((char*)&size, sizeof(size_t));
    tmp.resize(size);
    fin.read((char*)&tmp[0], size);
    table.push_back(tmp);
    if (fin.peek() == EOF ) break;
  }
  fin.close();
  if(FF_VEBOSE) cout << "loadSparse " << fname << " takes " << ((double)clock() - t)/CLOCKS_PER_SEC << " s\n";
}

void dump(const string& fname, const string& outname="k30.all.tsv")
{
  clock_t t;
  ifstream fin(fname, ios::binary);
  ofstream fout(outname);
  string l1, l2;
  while(!fin.eof())
  {
    if (fin.peek() == EOF ) break; loadSparseLine(fin, l1, l2);
    fout << l1 << "\t" << l2 << "\n";
  }
  fin.close();
  fout.close();
  if(FF_VEBOSE) cout << "dump takes " << ((double)clock() - t)/CLOCKS_PER_SEC << " s\n";
}

void merge(
  const StringVector& kmer1,
  const StringVector& table1,
  const StringVector& kmer2,
  const StringVector& table2,
  const string& outname)
{
  clock_t t;
  ofstream fout;
  fout.open(outname, ios::binary);
  int p=0, q=0;
  while(p < kmer1.size() && q < kmer2.size())
  {
    if(kmer1[p] == kmer2[q])
    {
      size_t size = kmer1[p].size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer1[p][0], size);
      string tmp = table1[p] + '\t' + table2[q];
      size = tmp.size();
      fout.write((char*)&size, sizeof(size_t));
      fout.write((char*)&tmp[0],  size);
      p++; q++;
    }
    else if (kmer1[p] < kmer2[q])
    {
      size_t size = kmer1[p].size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer1[p][0], size);
      size = table1[p].size();
      fout.write((char*)&size,         sizeof(size_t));
      fout.write((char*)&table1[p][0], size);
      p++;
    }
    else
    {
      size_t size = kmer2[q].size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer2[q][0], size);
      size = table2[q].size();
      fout.write((char*)&size,         sizeof(size_t));
      fout.write((char*)&table2[q][0], size);
      q++;
    }
  }
  while(p < kmer1.size())
  {
    size_t size = kmer1[p].size();
    fout.write((char*)&size,        sizeof(size_t));
    fout.write((char*)&kmer1[p][0], size);
    size = table1[p].size();
    fout.write((char*)&size,         sizeof(size_t));
    fout.write((char*)&table1[p][0], size);
    p++;
  }

  while(q < kmer2.size())
  {
    size_t size = kmer2[q].size();
    fout.write((char*)&size,        sizeof(size_t));
    fout.write((char*)&kmer2[q][0], size);
    size = table2[q].size();
    fout.write((char*)&size,         sizeof(size_t));
    fout.write((char*)&table2[q][0], size);
    q++;
  }
  fout.close();
  if(FF_VEBOSE) cout << "merge two table takes " << ((double)clock() - t)/CLOCKS_PER_SEC << " s\n";
}

void merge(
  const string& fname1,
  const StringVector& kmer2,
  const StringVector& table2,
  const string& outname)
{
  clock_t t;
  ifstream fin1(fname1, ios::binary);
  ofstream fout(outname, ios::binary);
  int q=0;
  string kmer1, row1;
  loadSparseLine(fin1, kmer1, row1);
  while(!fin1.eof() && q < kmer2.size())
  {
    if(kmer1 == kmer2[q])
    {
      size_t size = kmer1.size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer1[0], size);
      string tmp = row1 + "\t" + table2[q];
      size = tmp.size();
      fout.write((char*)&size, sizeof(size_t));
      fout.write((char*)&tmp[0],  size);
      if (fin1.peek() == EOF) break; loadSparseLine(fin1, kmer1, row1);
      q++;
    }
    else if (kmer1 < kmer2[q])
    {
      size_t size = kmer1.size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer1[0], size);
      size = row1.size();
      fout.write((char*)&size,         sizeof(size_t));
      fout.write((char*)&row1[0], size);
      if (fin1.peek() == EOF) break; loadSparseLine(fin1, kmer1, row1);
    }
    else
    {
      size_t size = kmer2[q].size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer2[q][0], size);
      size = table2[q].size();
      fout.write((char*)&size,         sizeof(size_t));
      fout.write((char*)&table2[q][0], size);
      q++;
    }
  }

  while(!fin1.eof())
  {
    size_t size = kmer1.size();
    fout.write((char*)&size,        sizeof(size_t));
    fout.write((char*)&kmer1[0], size);
    size = row1.size();
    fout.write((char*)&size,         sizeof(size_t));
    fout.write((char*)&row1[0], size);
    if (fin1.peek() == EOF) break; loadSparseLine(fin1, kmer1, row1);
  }
  while(q < kmer2.size())
  {
    size_t size = kmer2[q].size();
    fout.write((char*)&size,        sizeof(size_t));
    fout.write((char*)&kmer2[q][0], size);
    size = table2[q].size();
    fout.write((char*)&size,         sizeof(size_t));
    fout.write((char*)&table2[q][0], size);
    q++;
  }
  fout.close();
  if(FF_VEBOSE) cout << "merge two table takes " << ((double)clock() - t)/CLOCKS_PER_SEC << " s\n";

}

void merge(string& fname1, string& fname2, string& outname)
{
  ifstream fin1, fin2;
  ofstream fout;
  fin1.open(fname1, ios::binary);
  fin2.open(fname2, ios::binary);
  fout.open(outname, ios::binary);
  string kmer1, row1, kmer2, row2;
  loadSparseLine(fin1, kmer1, row1);
  loadSparseLine(fin2, kmer2, row2);
  while(!fin1.eof() && !fin2.eof())
  {
    if(kmer1 == kmer2)
    {
      size_t size = kmer1.size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer1[0], size);
      string tmp = row1 + "\t" + row2;
      size = tmp.size();
      fout.write((char*)&size, sizeof(size_t));
      fout.write((char*)&tmp[0],  size);
      if (fin1.peek() == EOF) break; loadSparseLine(fin1, kmer1, row1);
      if (fin2.peek() == EOF) break; loadSparseLine(fin2, kmer2, row2);
    }
    else if (kmer1 < kmer2)
    {
      size_t size = kmer1.size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer1[0], size);
      size = row1.size();
      fout.write((char*)&size,         sizeof(size_t));
      fout.write((char*)&row1[0], size);
      if (fin1.peek() == EOF) break; loadSparseLine(fin1, kmer1, row1);
    }
    else
    {
      size_t size = kmer2.size();
      fout.write((char*)&size,        sizeof(size_t));
      fout.write((char*)&kmer2[0], size);
      size = row2.size();
      fout.write((char*)&size,         sizeof(size_t));
      fout.write((char*)&row2[0], size);
      if (fin2.peek() == EOF) break; loadSparseLine(fin2, kmer2, row2);
    }
  }
  while(!fin1.eof())
  {
    size_t size = kmer1.size();
    fout.write((char*)&size,        sizeof(size_t));
    fout.write((char*)&kmer1[0], size);
    size = row1.size();
    fout.write((char*)&size,         sizeof(size_t));
    fout.write((char*)&row1[0], size);
    if (fin1.peek() == EOF) break; loadSparseLine(fin1, kmer1, row1);
  }
  while(!fin2.eof())
  {
    size_t size = kmer2.size();
    fout.write((char*)&size,        sizeof(size_t));
    fout.write((char*)&kmer2[0], size);
    size = row2.size();
    fout.write((char*)&size,         sizeof(size_t));
    fout.write((char*)&row2[0], size);
    if (fin2.peek() == EOF) break; loadSparseLine(fin2, kmer2, row2);
  }
  fin1.close();
  fin2.close();
  fout.close();
}


#endif
