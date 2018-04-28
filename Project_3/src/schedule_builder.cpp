#include "schedule_builder.h"
using namespace std;

void buildSchedule(schedule &schedule, const vector<pair<task*, task*>> &path, JSSP *jssp) {

    // Add tasks to machines in correct order
    schedule.machine_schedules.resize((unsigned)jssp->getNumMachines());
    for(auto &machine_schedule: schedule.machine_schedules){
        machine_schedule.clear();
    }

    for (auto &edge: path){
        task* current_task = edge.second;
        schedule_block block;
        block.task = current_task;

        double job_end_time = 0;
        double machine_end_time = 0;

        bool earlier_slot_available = false;
        double available_start_time = 0;
        int available_start_pos = 0;
        double prev_task_end_time = 0;

        // If there are tasks that come before the current task in this job (current task is not at column 0 in job_tasks)
        if ((block.task->task_id - 1) % jssp->getNumMachines() > 0){
            // Save when the preceding task ended
            pair<int, int> job_position = findJobInMachineSchedules(block.task->task_id - 1, schedule.machine_schedules);
            job_end_time = schedule.machine_schedules[job_position.first][job_position.second].start_time +
                           schedule.machine_schedules[job_position.first][job_position.second].task->process_time;
        }
        // If there are tasks that come before the current task in this machine (current task is not at column 0 in machine_schedules)
        if(!schedule.machine_schedules[current_task->machine_no].empty()){
            // Save when the preceding task ended
            machine_end_time = schedule.machine_schedules[current_task->machine_no].back().start_time +
                               schedule.machine_schedules[current_task->machine_no].back().task->process_time;

            earlier_slot_available = false;
            for (int i = 0; i < schedule.machine_schedules[current_task->machine_no].size(); i++){
                if (i == 0){
                    prev_task_end_time = 0;
                }
                else {
                    prev_task_end_time = schedule.machine_schedules[current_task->machine_no][i - 1].start_time
                                         + schedule.machine_schedules[current_task->machine_no][i - 1].task->process_time;
                }
                if (job_end_time < schedule.machine_schedules[current_task->machine_no][i].start_time){
                    double window_size = min(schedule.machine_schedules[current_task->machine_no][i].start_time - prev_task_end_time,
                                             schedule.machine_schedules[current_task->machine_no][i].start_time - job_end_time);
                    if (window_size >= block.task->process_time){
                        available_start_time = max(prev_task_end_time, job_end_time);
                        available_start_pos = i;
                        earlier_slot_available = true;
                        break;
                    }
                }
            }
        }

        if (earlier_slot_available){
            block.start_time = available_start_time;
            auto it = schedule.machine_schedules[current_task->machine_no].begin();
            schedule.machine_schedules[current_task->machine_no].insert(it + available_start_pos, block);
        }
        else{
            block.start_time = max(job_end_time, machine_end_time);
            schedule.machine_schedules[current_task->machine_no].push_back(block);
        }
    }

//    for (auto &edge: path){
//        task* task = edge.second;
//        schedule_block block;
//        block.task = task;
//        if (schedule.machine_schedules[task->machine_no].empty()){
//            block.start_time = 0;
//        }
//        else{
//            block.start_time = schedule.machine_schedules[task->machine_no].back().start_time + schedule.machine_schedules[task->machine_no].back().task->process_time;
//        }
//        schedule.machine_schedules[task->machine_no].push_back(block);
//    }
//
//    // Avoid collisions between tasks in the same job
//    for (auto &edge: path){
//        task* current_task = edge.second;
//
//        double job_end_time = 0;
//        double machine_end_time = 0;
//        for (int block_no = 0; block_no < schedule.machine_schedules[current_task->machine_no].size(); block_no++){
//            // Find the current task in machine_schedules
//            if (schedule.machine_schedules[current_task->machine_no][block_no].task->task_id == current_task->task_id){
//
//                // If there are tasks that come before the current task in this job (current task is not at column 0 in job_tasks)
//                if ((schedule.machine_schedules[current_task->machine_no][block_no].task->task_id - 1) % jssp->getNumMachines() > 0){
//                    // Save when the preceding task ended
//                    pair<int, int> job_position = findJobInMachineSchedules(schedule.machine_schedules[current_task->machine_no][block_no].task->task_id - 1, schedule.machine_schedules);
//                    job_end_time = schedule.machine_schedules[job_position.first][job_position.second].start_time +
//                                   schedule.machine_schedules[job_position.first][job_position.second].task->process_time;
//                }
//                // If there are tasks that come before the current task in this machine (current task is not at column 0 in machine_schedules)
//                if (block_no > 0){
//                    // Save when the preceding task ended
//                    machine_end_time = schedule.machine_schedules[current_task->machine_no][block_no - 1].start_time +
//                                       schedule.machine_schedules[current_task->machine_no][block_no - 1].task->process_time;
//                }
//                schedule.machine_schedules[current_task->machine_no][block_no].start_time = max(job_end_time, machine_end_time);
//            }
//        }
//    }

    // Find makespan of schedule
    double latest_finish = 0;
    for (auto &machine: schedule.machine_schedules){
        schedule_block last_task = machine.back();
        if (last_task.start_time + last_task.task->process_time > latest_finish){
            latest_finish = last_task.start_time + last_task.task->process_time;
        }
    }
    schedule.makespan = latest_finish;
}

pair<int, int> findJobInMachineSchedules(int task_id, const vector<vector<schedule_block>> &machine_schedules) {
    // Returns the position of task_id in machine_schedules
    for (int machine_no = 0; machine_no < machine_schedules.size(); machine_no++){
        for (int block_no = 0; block_no < machine_schedules[machine_no].size(); block_no++){
            if (machine_schedules[machine_no][block_no].task->task_id == task_id){
                pair<int,int> position = make_pair(machine_no, block_no);
                return position;
            }
        }
    }
}

void saveScheduleAsCSV(schedule &schedule, string filename, JSSP *jssp) {
    ofstream file;
    file.open ("../solutions/" + filename + ".csv");
    file << "Task,Start,Finish,Job,Description\n";

    for (auto &machine_schedule: schedule.machine_schedules){
        for (auto &schedule_block: machine_schedule){
            file << "Machine " << schedule_block.task->machine_no << ",";
            file << schedule_block.start_time << "," << schedule_block.start_time + schedule_block.task->process_time << ",";
            file << "Job " << schedule_block.task->job_id << "," << schedule_block.task->job_id << "/" << (schedule_block.task->task_id - 1) % jssp->getNumMachines() << "\n";
        }
    }
    file.close();
}
