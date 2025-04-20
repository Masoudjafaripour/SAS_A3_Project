/*
Copyright (c) 2023 Grid-based Path Planning Competition and Contributors <https://gppc.search-conference.org/>

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

#ifndef GPPC_ENTRY_H
#define GPPC_ENTRY_H

#include <vector>
#include <string>

// include common used class in GPPC
#include "GPPC.h"

typedef GPPC::xyLoc xyLoc;

void PreprocessMap(const std::vector<bool> &bits, int width, int height, const std::string &filename);
void *PrepareForSearch(const std::vector<bool> &bits, int width, int height, const std::string &filename);

void summarize_expansions();
extern std::vector<int> expansion_log;
extern long long heuristic_micro_time_global;


/*
return true if the pathfinding is completed (even if not path exist), 
usually this function always return true;

return false if the pathfinding is not completed and requires further function calls, 
e.g.:
  the shortest path from s to g is <s,v1,v2,g>
  GetPath(data, s, g, path);  // get the prefix <s, v1>
  GetPath(data, v1, g, path); // get the prefix <s,v1,v2>
  GetPath(data, v2, g, path); // get the entire <s,v1,v2,g>
  
*/
bool GetPath(void *data, xyLoc s, xyLoc g, std::vector<xyLoc> &path);

std::string GetName();

#endif // GPPC_ENTRY_H
