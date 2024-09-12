/*

MIT License 
Copyright (c) 2024 Junyan Dai and Erin Molloy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
#ifndef READ_SEARCH_SPACE
#define READ_SEARCH_SPACE
#include<fstream>
#include<cstdio>
#include<iostream>
#include<regex>
#include<string>
#include<set>
#include<vector>
#include<dirent.h>
#include "whereami++.hpp"
#include "whereami.hpp"
#include "utility.hpp"
#include <boost/unordered_set.hpp>
#include<unordered_set>
#include<boost/algorithm/string.hpp>

#endif
std::string findAstralJar() {
  std::string edir = whereami::executable_dir();
  std::cout << "Executable in " << edir << std::endl;
  DIR* d = opendir((edir + std::string("/Astral")).c_str());
  if (d == NULL)
    std::cout << "Directory " << edir + std::string("/Astral") << " could not be accessed" << std::endl;
  std::cout << strerror(errno) << std::endl;
  dirent* entry;
  while(d && (entry = readdir(d))) {
    std::string fname = entry->d_name;
    std::cout << "testing " << fname << std::endl;
    if (fname.substr(0, 7) == "astral." && fname.substr(fname.size() - 4) == ".jar") {
      std::cout << "Found" << std::endl;
      return edir + std::string("/Astral/") + fname;
    }
  }
  std::cout << "Couldn't find ASTRAL jar" << std::endl;
  return "";
}


void get_clade(std::unordered_set<Bipartition> &clades, boost::dynamic_bitset<> &Os,std::string line, boost::unordered_map<std::string, unsigned int> &label2index, std::vector<std::string> &labels, std::vector<std::string> &outgroup) {
  
  line = line.substr(1,line.size()-2);
  
  std::vector<std::string> taxa;
  boost::split(taxa, line, boost::is_any_of("{,}\t "),boost::token_compress_on);
  
  boost::dynamic_bitset<> bs(labels.size());
  
  for (auto t: taxa) {
    bs.set(label2index[t]);
  }
  
  Bipartition bp(bs);
  unsigned int cond = 0;
  for (auto t: outgroup) {
    if (bs[label2index[t]] == 1) cond++;
  }

  if (cond == 0 || cond == bs.count()) {
    clades.insert(bp);
    
  } else if (cond == outgroup.size()) {
    bp.flip();
    clades.insert(bp);
    
  } 
}


std::unordered_set<Bipartition> read_search_space(std::string const& filename, boost::unordered_map<std::string, unsigned int> &label2index, std::vector<std::string> &labels, std::vector<std::string> &outgroup, std::string memory) {
  std::stringstream clade_stream;
  std::string astralpath = findAstralJar();
#ifdef _WIN32
  std::string s = "java " + memory + " -jar " + astralpath + " -i " + filename + " -k searchspace_norun -o NUL";
  //if (guided_tree != "") {
//	s = "java -jar " + astralpath + " -i " + filename + " -j " + guided_tree + " -k searchspace_norun -o NUL";
  //}
#else
  std::string s = "java " + memory + " -jar " + astralpath + " -i " + filename + " -k searchspace_norun -o /dev/null";
  //if (guided_tree != "") {
  //	s = "java -jar " + astralpath + " -i " + filename + " -j " + guided_tree + " -k searchspace_norun -o /dev/null";
  //}
#endif
  char buffer[128];

  FILE* stream = popen(s.c_str(), "r");
  std::stringstream result;
  while(!feof(stream)) {
    if(fgets(buffer,128,stream) != NULL) {
      result << buffer;
    }
  }

  
  std::stringstream cladestream_mapped(result.str());
  std::string line;

  std::stringstream unmapped;
  std::unordered_set<Bipartition> clades;
  boost::dynamic_bitset<> Outs(labels.size());

  for (auto e : outgroup) {
    Outs.set(label2index[e]);
  }
  
  while (!cladestream_mapped.eof()) {
    getline(cladestream_mapped, line);
    if (line[0] != '{') {
      continue;
    }
    get_clade(clades, Outs, line, label2index, labels, outgroup);
  }

  boost::dynamic_bitset<> tbs(labels.size());
  tbs.flip();
  Bipartition triv(tbs);
  clades.insert(triv);
  return clades;
          
}
