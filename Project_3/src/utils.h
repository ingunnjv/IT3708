#ifndef PROJECT_3_UTILS_H
#define PROJECT_3_UTILS_H

#include <utility>
#include <vector>
#include <string>

std::pair<int, int> id2xy(int id);
bool isElementInVector(int element, std::vector<int> vec);
void callPythonGanttChartPlotter(std::string solutionFileName);
void printScoreToScreen(double best_makespan, double optimal_makespan);

#endif //PROJECT_3_UTILS_H

