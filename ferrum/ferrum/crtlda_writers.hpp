#ifndef FERRUM_LIBNAR_CRTLDA_WRITERS_H_
#define FERRUM_LIBNAR_CRTLDA_WRITERS_H_

#include "ferrum/concrete.hpp"
#include "ferrum/crtlda_defs.hpp"
#include "ferrum/dmc.hpp"
#include "ferrum/mathops.hpp"
#include "ferrum/svi_util.hpp"
#include "ferrum/util.hpp"

#include <iostream>
#include <map>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ferrum {
  class DKVWriters {
private:
  ferrum::SmartWriter* sw_template_usage_ = NULL;
  ferrum::SmartWriter* sw_slot_usage_ = NULL;
  ferrum::SmartWriter* sw_frame_usage_ = NULL;
  ferrum::SmartWriter* sw_role_usage_ = NULL;
  ferrum::SmartWriter* sw_verb_usage_ = NULL;
  ferrum::SmartWriter* sw_arc_usage_ = NULL;
  //ferrum::SmartWriter* sw_labeler = NULL;
  ferrum::SmartWriter* sw_model_tsv_ = NULL;
  std::vector< ferrum::SmartWriter* > ptrs_;
  void reg_sw(ferrum::SmartWriter*& p, const std::string& name);
public:
  DKVWriters();
  ~DKVWriters();
  void create_template_usage(const std::string& name);
  void create_slot_usage(const std::string& name);
  void create_frame_usage(const std::string& name);
  void create_role_usage(const std::string& name);
  void create_verb_usage(const std::string& name); 
  void create_arc_usage(const std::string& name);
  void create_model_tsv(const std::string& name);
  bool to_file_template_usage();
  bool to_file_slot_usage();
  bool to_file_frame_usage();
  bool to_file_role_usage();
  bool to_file_verb_usage(); 
  bool to_file_arc_usage();
  bool to_model_tsv();
  ferrum::SmartWriter* sw_template_usage();
  ferrum::SmartWriter* sw_slot_usage();
  ferrum::SmartWriter* sw_frame_usage();
  ferrum::SmartWriter* sw_role_usage();
  ferrum::SmartWriter* sw_verb_usage(); 
  ferrum::SmartWriter* sw_arc_usage();
  ferrum::SmartWriter* sw_model_tsv();
};
}

#endif
