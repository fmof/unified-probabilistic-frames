#include "ferrum/crtlda_writers.hpp"

namespace ferrum {
  DKVWriters::DKVWriters() {
  }
  DKVWriters::~DKVWriters() {
    const size_t size = ptrs_.size();
    for(size_t i = 0; i < size; ++i) {
      delete (ptrs_[i]);
    }
  }
  void DKVWriters::reg_sw(ferrum::SmartWriter*& p, const std::string& name) {
    p = new ferrum::SmartWriter(name);
    ptrs_.push_back(p);
  }
  void DKVWriters::create_template_usage(const std::string& name) {
    return reg_sw(sw_template_usage_, name);
  }
  void DKVWriters::create_slot_usage(const std::string& name) {
    return reg_sw(sw_slot_usage_, name);
  }
  void DKVWriters::create_frame_usage(const std::string& name) {
    return reg_sw(sw_frame_usage_, name);
  }
  void DKVWriters::create_role_usage(const std::string& name) {
    return reg_sw(sw_role_usage_, name);
  }
  void DKVWriters::create_verb_usage(const std::string& name) {
    return reg_sw(sw_verb_usage_, name);
  }
  void DKVWriters::create_arc_usage(const std::string& name) {
    return reg_sw(sw_arc_usage_, name);
  }
  void DKVWriters::create_model_tsv(const std::string& name) {
    return reg_sw(sw_model_tsv_, name);
  }
  bool DKVWriters::to_file_template_usage() {
    return sw_template_usage_ != NULL && sw_template_usage_->to_file();
  }
  bool DKVWriters::to_file_frame_usage() {
    return sw_frame_usage_ != NULL && sw_frame_usage_->to_file();
  }
  bool DKVWriters::to_file_slot_usage() {
    return sw_slot_usage_ != NULL && sw_slot_usage_->to_file();
  }
  bool DKVWriters::to_file_role_usage() {
    return sw_role_usage_ != NULL && sw_role_usage_->to_file();
  }
  bool DKVWriters::to_file_verb_usage() {
    return sw_verb_usage_ != NULL && sw_verb_usage_->to_file();
  }
  bool DKVWriters::to_file_arc_usage() {
    return sw_arc_usage_ != NULL && sw_arc_usage_->to_file();
  }
  bool DKVWriters::to_model_tsv() {
    return sw_model_tsv_ != NULL && sw_model_tsv_->to_file();
  }
  ferrum::SmartWriter* DKVWriters::sw_template_usage() {
    return sw_template_usage_;
  }
  ferrum::SmartWriter* DKVWriters::sw_slot_usage() {
    return sw_slot_usage_;
  }
  ferrum::SmartWriter* DKVWriters::sw_frame_usage() {
    return sw_frame_usage_;
  }
  ferrum::SmartWriter* DKVWriters::sw_role_usage() {
    return sw_role_usage_;
  }
  ferrum::SmartWriter* DKVWriters::sw_verb_usage() {
    return sw_verb_usage_;
  }
  ferrum::SmartWriter* DKVWriters::sw_arc_usage() {
    return sw_arc_usage_;
  }
  ferrum::SmartWriter* DKVWriters::sw_model_tsv() {
    return sw_model_tsv_;
  }
}
