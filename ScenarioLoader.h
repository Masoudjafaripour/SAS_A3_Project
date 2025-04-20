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

#ifndef GPPC_SCENARIOLOADER_H
#define GPPC_SCENARIOLOADER_H

#include <vector>
#include <cstring>
#include <string>

static const int kNoScaling = -1;

/** 
 * Experiments stored by the ScenarioLoader class. 
 */
class ScenarioLoader;

class Experiment {
public:
	Experiment(int sx,int sy,int gx,int gy,int b, double d, std::string m)
    :startx(sx),starty(sy),goalx(gx),goaly(gy),scaleX(kNoScaling),scaleY(kNoScaling),bucket(b),distance(d),map(m){}
	Experiment(int sx,int sy,int gx,int gy,int sizeX, int sizeY,int b, double d, std::string m)
    :startx(sx),starty(sy),goalx(gx),goaly(gy),scaleX(sizeX),scaleY(sizeY),bucket(b),distance(d),map(m){}
	int GetStartX() const {return startx;}
	int GetStartY() const {return starty;}
	int GetGoalX() const {return goalx;}
	int GetGoalY() const {return goaly;}
	int GetBucket() const {return bucket;}
	double GetDistance() const {return distance;}
	void GetMapName(char* mymap) const {strcpy(mymap,map.c_str());}
	const char *GetMapName() const { return map.c_str(); }
	int GetXScale() const {return scaleX;}
	int GetYScale() const {return scaleY;}
	
private:
	friend class ScenarioLoader;
	int startx, starty, goalx, goaly;
	int scaleX;
	int scaleY;
	int bucket;
	double distance;
	std::string map;
};

/** A class which loads and stores scenarios from files.  
 * Versions currently handled: 0.0 and 1.0 (includes scale). 
 */

class ScenarioLoader{
public:
	ScenarioLoader() { scenName[0] = 0; }
	ScenarioLoader(const char *);
	void Save(const char *);
	int GetNumExperiments(){return experiments.size();}
	const char *GetScenarioName() { return scenName; }
	Experiment GetNthExperiment(int which)
	{return experiments[which];}
	void AddExperiment(Experiment which);
private:
	char scenName[1024];
	std::vector<Experiment> experiments;
};

#endif // GPPC_SCENARIOLOADER_H
