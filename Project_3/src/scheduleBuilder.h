#ifndef PROJECT_3_SCHEDULEBUILDER_H
#define PROJECT_3_SCHEDULEBUILDER_H

#include "jssp.h"

struct schedule{
    double makespan;
    std::string filename;
    std::vector<std::vector<schedule_block>> machine_schedules;
};


void buildSchedule(schedule &schedule, const std::vector<std::pair<task *, task *>> &path, int k, JSSP *jssp);
std::pair<int, int> findJobInMachineSchedules(int task_id, const std::vector<std::vector<schedule_block>> &machine_schedules);
void saveScheduleAsCSV(schedule &schedule, JSSP *jssp);



#endif //PROJECT_3_SCHEDULEBUILDER_H
