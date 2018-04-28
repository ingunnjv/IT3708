#ifndef PROJECT_3_SCHEDULE_BUILDER_H
#define PROJECT_3_SCHEDULE_BUILDER_H

#include "jssp.h"


struct schedule_block{
    task* task;
    double start_time;
};

struct schedule{
    double makespan;
    std::vector<std::vector<schedule_block>> machine_schedules;
};

void buildSchedule(schedule &schedule, const std::vector<std::pair<task *, task *>> &path, JSSP *jssp);
std::pair<int, int> findJobInMachineSchedules(int task_id, const std::vector<std::vector<schedule_block>> &machine_schedules);
void saveScheduleAsCSV(schedule &schedule, std::string filename, JSSP *jssp);



#endif //PROJECT_3_SCHEDULE_BUILDER_H
