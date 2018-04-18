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