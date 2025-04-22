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

#include <cstdio>
#include <ios>
#include <numeric>
#include <algorithm>
#include <string>
#include <unistd.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "ScenarioLoader.h"
#include "Timer.h"
#include "Entry.h"
#include "validator/ValidatePath.hpp"
#include <filesystem>  // Add at top


std::string datafile, mapfile, scenfile, flag;
const std::string index_dir = "index_data";
constexpr double PATH_FIRST_STEP_LENGTH = 20.0;
std::vector<bool> mapData;
int width, height;
bool pre   = false;
bool run   = false;
bool check = false;

namespace fs = std::filesystem;


void LoadMap(const char *fname, std::vector<bool> &map, int &width, int &height)
{
  FILE *f;
  f = std::fopen(fname, "r");
  if (f)
  {
    std::fscanf(f, "type octile\nheight %d\nwidth %d\nmap\n", &height, &width);
    map.resize(height*width);
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        char c;
        do {
          std::fscanf(f, "%c", &c);
        } while (std::isspace(c));
        map[y*width+x] = (c == '.' || c == 'G' || c == 'S');
      }
    }
    std::fclose(f);
  }
}

double euclidean_dist(const xyLoc& a, const xyLoc& b) {
  int dx = std::abs(b.x - a.x);
  int dy = std::abs(b.y - a.y);
  double res = std::sqrt(dx * dx + dy * dy);
  return res;
}

double GetPathLength(const std::vector<xyLoc>& path)
{
  double len = 0;
  for (int x = 0; x < (int)path.size()-1; x++)
    len += euclidean_dist(path[x], path[x+1]);
  return len;
}

// returns -1 if valid path, otherwise id of segment where invalidness was detetcted
int ValidatePath(const std::vector<xyLoc>& thePath)
{
  return inx::ValidatePath(mapData, width, height, thePath);
}

int total_num_experiments = 0;

int RunExperiment(void* data) {
  Timer t;
  ScenarioLoader scen(scenfile.c_str());
  std::vector<xyLoc> thePath;

  // std::string resultfile = "result.csv";
  // std::ofstream fout(resultfile);
  // const std::string header = "map,scen,experiment_id,path_size,path_length,ref_length,time_cost,20steps_cost,max_step_time";
  // fout << header << std::endl;

  for (int x = 0; x < scen.GetNumExperiments(); x++)
  {
    xyLoc s, g;
    s.x = scen.GetNthExperiment(x).GetStartX();
    s.y = scen.GetNthExperiment(x).GetStartY();
    g.x = scen.GetNthExperiment(x).GetGoalX();
    g.y = scen.GetNthExperiment(x).GetGoalY();

    thePath.clear();
    typedef Timer::duration dur;
    dur max_step = dur::zero(), tcost = dur::zero(), tcost_first = dur::zero();
    bool done = false, done_first = false;
    do {
      t.StartTimer();
      done = GetPath(data, s, g, thePath);
      t.EndTimer();
      max_step = std::max(max_step, t.GetElapsedTime());
      tcost += t.GetElapsedTime();
      if (!done_first) {
        tcost_first += t.GetElapsedTime();
        done_first = GetPathLength(thePath) >= PATH_FIRST_STEP_LENGTH - 1e-6;
      }
    } while (!done);
    double plen = done?GetPathLength(thePath): 0;
    double ref_len = scen.GetNthExperiment(x).GetDistance();


    // fout << std::setprecision(9) << std::fixed;
    // fout << mapfile  << "," << scenfile       << ","
    //      << x        << "," << thePath.size() << ","
    //      << plen     << "," << ref_len        << ","
    //      << tcost.count() << "," << tcost_first.count() << ","
    //      << max_step.count() << std::endl;

    if (check) {
      std::printf("%d %d %d %d", s.x, s.y, g.x, g.y);
      int validness = ValidatePath(thePath);
      if (validness < 0) {
        std::printf(" valid");
      } else {
        std::printf(" invalid-%d", validness);
      }
      std::printf(" %d", static_cast<int>(thePath.size()));
      for (const auto& it: thePath) {
        std::printf(" %d %d", it.x, it.y);
      }
      std::printf(" %.5f\n", plen);
    }

  }
  return scen.GetNumExperiments();
}

void print_help(char **argv) {
  std::printf("Invalid Arguments\nUsage %s <flag> <map> <scenario>\n", argv[0]);
  std::printf("Flags:\n");
  std::printf("\t-full : Preprocess map and run scenario\n");
  std::printf("\t-pre : Preprocess map\n");
  std::printf("\t-run : Run scenario without preprocessing\n");
  std::printf("\t-check: Run for validation\n");
}

bool parse_argv_1(int argc, char **argv) {
  if (argc < 2) return false;
  flag = std::string(argv[1]);
  if (flag == "-full") pre = run = true;
  else if (flag == "-pre") pre = true;
  else if (flag == "-run") run = true;
  else if (flag == "-check") run = check = true;
  else return false;

  if (argc < 3) return false;
  mapfile = std::string(argv[2]);

  if (run) {
    if (argc < 4) return false;
    scenfile = std::string(argv[3]);
  }
  return true;
}


bool parse_argv(int argc, char **argv) {
  if (argc < 2) return false;
  flag = std::string(argv[1]);
  if (flag == "-full") pre = run = true;
  else if (flag == "-pre") pre = true;
  else if (flag == "-run") run = true;
  else if (flag == "-check") run = check = true;
  else return false;

  if (argc < 3) return false;
  mapfile = std::string(argv[2]);

  // ✅ Allow -run with only a folder path
  if (run && !fs::is_directory(mapfile)) {
    if (argc < 4) return false;
    scenfile = std::string(argv[3]);
  }

  return true;
}

std::string basename(const std::string& path) {
  std::size_t l = path.find_last_of('/');
  if (l == std::string::npos) l = 0;
  else l += 1;
  std::size_t r = path.find_last_of('.');
  if (r == std::string::npos) r = path.size()-1;
  return path.substr(l, r-l);
}

std::vector<int> all_expansions;

std::string parent_folder_name(const std::string& path) {
  std::size_t last_slash = path.find_last_of('/');
  if (last_slash == std::string::npos || last_slash == 0) return "";
  
  std::size_t second_last_slash = path.find_last_of('/', last_slash - 1);
  if (second_last_slash == std::string::npos)
    return path.substr(0, last_slash);
  else
    return path.substr(second_last_slash + 1, last_slash - second_last_slash - 1);
}


void summarize_and_export_stats(const std::vector<int>& log,
                                 const std::string& csv_filename,
                                 int total_problems,
                                 long long heuristic_time_ms)
{
    if (log.empty()) return;

    int n = log.size();
    double mean = accumulate(log.begin(), log.end(), 0.0) / n;

    std::vector<int> sorted = log;
    std::sort(sorted.begin(), sorted.end());
    double median = (n % 2 == 0)
        ? (sorted[n/2 - 1] + sorted[n/2]) / 2.0
        : sorted[n/2];

    double sum_sq = 0.0;
    for (int x : log) sum_sq += (x - mean) * (x - mean);
    double stdev = sqrt(sum_sq / (n - 1));
    double ci_95 = 1.96 * stdev / sqrt(n);

    // Save to CSV
    std::ofstream fout(csv_filename);
    fout << "total_problems,mean,median,ci_95,heuristic_time_ms\n";
    fout << total_problems << ","
         << mean << ","
         << median << ","
         << ci_95 << ","
         << heuristic_time_ms << "\n";
    fout.close();
}


int main(int argc, char **argv)
{
  if (!parse_argv(argc, argv)) {
    print_help(argv);
    return 1;
  }

  // redirect stdout/stderr
  std::freopen("run.stdout", "w", stdout);
  std::freopen("run.stderr", "w", stderr);

  long long total_heuristic_time = 0;

  if (flag == "-run" && fs::is_directory(mapfile)) {
    // Loop through all .map files in the directory
    for (const auto& entry : fs::directory_iterator(mapfile)) {
      if (entry.path().extension() == ".map") {
        std::string map_path = entry.path().string();
        std::string scen_path = map_path + ".scen";

        if (!fs::exists(scen_path)) {
          std::cerr << "Skipping " << map_path << ": no .scen file found.\n";
          continue;
        }

        mapfile = map_path;
        scenfile = scen_path;

        // Load and run as usual
        mapData.clear();
        LoadMap(mapfile.c_str(), mapData, width, height);
        datafile = index_dir + "/" + GetName() + "-" + basename(mapfile);
        void* reference = PrepareForSearch(mapData, width, height, datafile);
        // RunExperiment(reference);

        int num_exps = RunExperiment(reference);
        total_num_experiments += num_exps;

        //  Summary for this map
        summarize_expansions();
        all_expansions.insert(all_expansions.end(), expansion_log.begin(), expansion_log.end());
        expansion_log.clear();  // Reset for next map

        total_heuristic_time += heuristic_micro_time_global;

        if (total_num_experiments > 20000) { // more than actual 20k 
          std::cerr << "total_num_experiments is" << total_num_experiments << std::endl;
          break;
        }
      }
    }


    std::cout << "Total number of experiments: " << total_num_experiments << std::endl;
    std::string parent_folder = parent_folder_name(mapfile);  // mapfile is "benchmarks/temp"
    std::string output_csv = parent_folder + "_summary.csv";
    // std::cout << "\n=== Global Summary Across All Maps ===\n";
    summarize_and_export_stats(all_expansions, output_csv, all_expansions.size(), total_heuristic_time);


    return 0;
  }

  // Default single-run behavior
  LoadMap(mapfile.c_str(), mapData, width, height);
  datafile = index_dir + "/" + GetName() + "-" + basename(mapfile);

  if (pre)
    PreprocessMap(mapData, width, height, datafile);

  if (!run)
    return 0;

  void *reference = PrepareForSearch(mapData, width, height, datafile);

  char argument[256];
  std::sprintf(argument, "pmap -x %d | tail -n 1 > run.info", getpid());
  std::system(argument);
  RunExperiment(reference);
  std::sprintf(argument, "pmap -x %d | tail -n 1 >> run.info", getpid());
  std::system(argument);
  return 0;
}

