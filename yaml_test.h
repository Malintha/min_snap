#ifndef yamltest_h
#define yamltest_h

#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <yaml.h>
#include "Trajectory.h"

using namespace std;
using namespace Eigen;

vector<Trajectory> processYamlFile(char* fPath, int horizon_id) {
    FILE *fh = fopen(fPath, "r");
    yaml_parser_t parser;
    yaml_token_t token;
    yaml_event_t event;

    if (!yaml_parser_initialize(&parser)) {
        std::cout << "Failed to initialize parser!\n" << std::endl;
    }
    if (fh == NULL) {
        std::cout << "Failed to open file!\n" << std::endl;
    }

    yaml_parser_set_input_file(&parser, fh);

    int subgoal_id = 0;
    int drone_id = 0;
    bool drone_node = false;
    bool subgoalset_node = false;
    bool subgoal_node = false;
    bool subgoals_node = false;
    int subgoals = 0;
    bool times_node = false;
    bool timeset_node = false;
    int timeset_id = 0;

    vector<double> pos_sequence;
    vector<Trajectory> tr_list;
    vector<double> t_list;
  
    do {
        if (!yaml_parser_parse(&parser, &event)) {
        printf("Parser error %d\n", parser.error);
        exit(EXIT_FAILURE);
        }

        switch (event.type) {
            case YAML_SEQUENCE_START_EVENT:
                subgoal_node = true;

            break;
            
            case YAML_SEQUENCE_END_EVENT:
                subgoal_node = false;
                if(timeset_node)
                    timeset_node = false;
            break;

            case YAML_MAPPING_START_EVENT:
            break;
            
            case YAML_MAPPING_END_EVENT:
                if(subgoalset_node && drone_node) {
                    subgoalset_node = false;
                }
                else if(drone_node) {
                    drone_node = false;
                }
            break;

            case YAML_SCALAR_EVENT:
                char *c = (char *)event.data.scalar.value;
                string s = std::string(c);
                const std::regex drone_regex("drone(\\d)");
                const std::regex subgoalset_regex("subgoalset(\\d)");
                const std::regex subgoal_regex("subgoal");
                const std::regex subgoals_regex("subgoals");
                const std::regex times_regex("times");
                const std::regex timeset_regex("timeset(\\d)");
                std::smatch pieces_match;

                if(std::regex_match(s, pieces_match,drone_regex)) {
                    drone_node = true;
                    ssub_match sub_match = pieces_match[1];
                    std::string piece = sub_match.str();
                    drone_id = std::stoi(piece);
                    continue;
                }
                else if(std::regex_match(s, pieces_match,subgoalset_regex)) {
                    subgoalset_node = true;
                    ssub_match sub_match = pieces_match[1];
                    std::string piece = sub_match.str();
                    subgoal_id = std::stoi(piece);
                    continue;
                }
                else if(std::regex_match(s, subgoals_regex)) {
                    subgoals_node = true;
                    continue;
                }
                else if(std::regex_match(s, times_regex)) {
                    times_node = true;
                    continue;
                }
                else if(std::regex_match(s, pieces_match, timeset_regex)) {
                    timeset_node = true;
                    ssub_match sub_match = pieces_match[1];
                    std::string piece = sub_match.str();
                    timeset_id = std::stoi(piece);
                    continue;
                }

                if(subgoals_node) {
                    subgoals = std::stoi(s);
                    subgoals_node = false;
                }

                if(timeset_node && times_node) {
                    if(timeset_id == horizon_id) {
                        t_list.push_back(std::stod(s));
                    }
                }

                if(subgoalset_node && subgoal_node) {
                    if(subgoal_id == horizon_id) {

                        if(pos_sequence.size() < 3*subgoals) {
                            pos_sequence.push_back(std::stod(s));
                        }
                        if(pos_sequence.size() == 3*subgoals) {
                            Eigen::Vector3d p;
                            Trajectory tr;
                            for(int i=0;i<pos_sequence.size();i++) {
                                int d = i%3;
                                if(d <= 2) {
                                    p[d] = pos_sequence[i];
                                }
                                if(d==2) {
                                    tr.pos.push_back(p);
                                }
                            }
                            tr.tList = t_list;
                            tr_list.push_back(tr);
                            pos_sequence.clear();
                        }
                    }
                }
        break;
        }
        if (event.type != YAML_STREAM_END_EVENT) {
            yaml_event_delete(&event);
        }
    } 
    while (event.type != YAML_STREAM_END_EVENT);
    yaml_event_delete(&event);
    return tr_list;
    }

#endif