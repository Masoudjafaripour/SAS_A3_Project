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

#include <fstream>
using std::ifstream;
using std::ofstream;

#include "ScenarioLoader.h"

/** 
 * Loads the experiments from the scenario file. 
 */
ScenarioLoader::ScenarioLoader(const char* fname)
{
	strncpy(scenName, fname, 1000);
  ifstream sfile(fname,std::ios::in);
  
  float ver;
  std::string first;
  sfile>>first;

  // Check if a version number is given
  if(first != "version"){
    ver = 0.0;
    sfile.seekg(0,std::ios::beg);
  }
  else{
    sfile>>ver;
  }

  int sizeX = 0, sizeY = 0; 
  int bucket;
  std::string map;  
  int xs, ys, xg, yg;
  double dist;

  // Read in & store experiments
  if (ver==0.0){
    while(sfile>>bucket>>map>>xs>>ys>>xg>>yg>>dist) {
      Experiment exp(xs,ys,xg,yg,bucket,dist,map);
      experiments.push_back(exp);
    }
  }
  else if(ver==1.0){
    while(sfile>>bucket>>map>>sizeX>>sizeY>>xs>>ys>>xg>>yg>>dist){
      Experiment exp(xs,ys,xg,yg,sizeX,sizeY,bucket,dist,map);
      experiments.push_back(exp);
    }
  }
  else{
    printf("Invalid version number.\n");
  }
}

void ScenarioLoader::Save(const char *fname)
{
//	strncpy(scenName, fname, 1024);
	ofstream ofile(fname);
	
	float ver = 1.0;
	ofile<<"version "<<ver<<std::endl;
	
	
	for (unsigned int x = 0; x < experiments.size(); x++)
	{
		ofile<<experiments[x].bucket<<"\t"<<experiments[x].map<<"\t"<<experiments[x].scaleX<<"\t";
		ofile<<experiments[x].scaleY<<"\t"<<experiments[x].startx<<"\t"<<experiments[x].starty<<"\t";
		ofile<<experiments[x].goalx<<"\t"<<experiments[x].goaly<<"\t"<<experiments[x].distance<<std::endl;
	}
}

void ScenarioLoader::AddExperiment(Experiment which)
{
	experiments.push_back(which);
}

