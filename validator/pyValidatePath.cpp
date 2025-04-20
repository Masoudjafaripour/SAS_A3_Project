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

#include <pybind11/pybind11.h>
#include "ValidatePath.hpp"

namespace py = pybind11;

struct xyLoc {
    int16_t x;
    int16_t y;
};

struct Checker{
    std::vector<bool> map;
    int width;
    int height;
    Checker(py::list& theMap, int width, int height):width(width),height(height)
    {
        map.resize(py::len(theMap));
        for(int i = 0 ; i<py::len(theMap);i++){
            map[i] = theMap[i].cast<bool>();
        }

    }
    int validatePath(py::list thePath){
        std::vector<xyLoc> path;
        path.resize(py::len(thePath));
        for(int i=0;i<py::len(thePath);i++){
            xyLoc loc;
            loc.x = thePath[i].attr("x").cast<int>();
            loc.y = thePath[i].attr("y").cast<int>();
            path[i] = loc;
        }

        return inx::ValidatePath(map, width, height, path);
    };
};




PYBIND11_MODULE(Grid_Path_Checker, m) {
    py::class_<Checker>(m, "Grid_Path_Checker")
        .def(py::init<py::list&, int, int>())
        .def("validatePath", &Checker::validatePath);
}


