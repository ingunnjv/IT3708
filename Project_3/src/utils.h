#ifndef PROJECT_3_UTILS_H
#define PROJECT_3_UTILS_H

#include <utility>
#include <vector>

std::pair<int, int> id2xy(int id);
bool isElementInVector(int element, std::vector<int> vec);
void callPythonGanttChartPlotter(string solutionFileName);

#endif //PROJECT_3_UTILS_H

