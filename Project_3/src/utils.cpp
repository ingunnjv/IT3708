#include "utils.h"
using namespace std;

std::pair<int, int> id2xy(int id, int num_cols){
    int x = id / num_cols;
    int y = id % num_cols;
    pair<int,int> xy = make_pair(x,y);
    return xy;
};

bool isElementInVector(int element, std::vector<int> vec){
    for (auto &el: vec){
        if (el == element){
            return true;
        }
    }
    return false;
}

void callPythonGanttChartPlotter(string solutionFileName){
    /* Create gantt chart from python script */
    printf("Print Gantt chart of best solution..\n");
    string command = "python \"..\\src\\run_gantt.py\"";
    string args = " " + solutionFileName;
    command += args;
    system(command.c_str());
}

void printScoreToScreen(double best_makespan, double optimal_makespan){
    printf("\nOptimizer algorithm terminated successfuly\n");
    printf("Shortest makespan obtained: %.1f\n", best_makespan);
    printf("Acceptable makespan value: %.1f\n", optimal_makespan);
    printf("The best solution obtained is within %.2f %% of the optimal solution\n", 100.0*(best_makespan/optimal_makespan - 1.0));
}