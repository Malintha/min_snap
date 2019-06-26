#include "solver.h"
#include <iostream>
#include "yaml_test.h"
#include "Trajectory.h"
#include <vector>
#include <algorithm>

using namespace std;

int main() {

    string yaml_fpath = "/home/malintha/Desktop/qptest/goals.yaml";
    char cstr[yaml_fpath.size()+1];
    copy(yaml_fpath.begin(), yaml_fpath.end(), cstr);
    cstr[yaml_fpath.size()] = '\0';
    vector<Trajectory> h1 = processYamlFile(cstr, 0);
    cout<<"h1: "<<h1[0].pos.size()<<endl;

    Solver s1(1, 4, 5, 0, 10);
    vector<Trajectory> res1 = s1.solve(h1, true);
    cout<<endl<<"####solved####"<<endl;

    vector<Trajectory> h2 = processYamlFile(cstr, 0);
    cout<<"h2: "<<h2[0].pos.size()<<endl;

    vector<Trajectory> res2 = s1.solve(h2, true);
    cout<<endl<<"####solved####"<<endl;
}