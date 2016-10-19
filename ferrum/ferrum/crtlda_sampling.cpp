#include "ferrum/crtlda_sampling.hpp"
#include "ferrum/version.hpp"

#include <Eigen/Dense>

namespace ferrum {
  void SamplingStrategy::heldout(bool b) {
    heldout_ = b;
  }
  bool SamplingStrategy::heldout() {
    return heldout_;
  }
  bool SamplingStrategy::reestimate_hyperparameters(int iteration) {
    return !heldout_ && (1+iteration) > burn_in && (1+iteration) % 200 == 0;
  }
  int SamplingStrategy::reestimate_template_usage_every() {
    return reestimate_template_usage_every_;
  }
  int SamplingStrategy::reestimate_slot_usage_every() {
    return reestimate_slot_usage_every_;
  }
  int SamplingStrategy::reestimate_gov_every() {
    return reestimate_gov_every_;
  }
  int SamplingStrategy::reestimate_rel_every() {
    return reestimate_rel_every_;
  }
  SamplingStrategy::SamplingStrategy(int num_iter, int burnin) :
    num_iterations((const int)num_iter), 
    burn_in((const int)burnin), heldout_(false) {
  }
  SamplingStrategy::~SamplingStrategy() {
  }
  int WithKindsSamplingStrategy::reestimate_gov_kind_every() {
    return reestimate_gov_every_;
  }
  int WithKindsSamplingStrategy::reestimate_rel_kind_every() {
    return reestimate_rel_every_;
  }
  WithKindsSamplingStrategy::WithKindsSamplingStrategy(int num_iter,
						       int burnin) :
    SamplingStrategy(num_iter, burnin) {
    reestimate_gov_kind_every_ = 100;
    reestimate_rel_kind_every_ = 100;
  }
  WithKindsSamplingStrategy::~WithKindsSamplingStrategy() {
  }
  SampleEveryIter::SampleEveryIter(int num_iter, int burnin) :
    SamplingStrategy(num_iter, burnin) {
    reestimate_template_usage_every_ = 100;
    reestimate_slot_usage_every_ = 100;
    reestimate_gov_every_ = 100;
    reestimate_rel_every_ = 100;
  }
  bool SampleEveryIter::sample_template(int i, int d, int e) {
    return true;
  }
  bool SampleEveryIter::sample_slot(int i, int d, int e) {
    return true;
  }
  bool SampleEveryIter::reestimate_template_usage(int iter_index) {
    return iter_index >= burn_in && (1+iter_index) % this->reestimate_template_usage_every() == 0;
  }
  bool SampleEveryIter::reestimate_slot_usage(int iter_index) {
    return !heldout_ && (1+iter_index) >= burn_in && (1+iter_index) % this->reestimate_slot_usage_every() == 0;
  }
  bool SampleEveryIter::reestimate_gov(int iter_index) {
    return !heldout_ && (1+iter_index) >= burn_in && (1+iter_index) % this->reestimate_gov_every() == 0;
  }
  bool SampleEveryIter::reestimate_rel(int iter_index) {
    return !heldout_ && (1+iter_index) >= burn_in && (1+iter_index) % this->reestimate_rel_every() == 0;
  }
  SampleWithKindsEveryIter::SampleWithKindsEveryIter(int num_iter, int burnin) : SamplingStrategy(num_iter, burnin), SampleEveryIter(num_iter, burnin), WithKindsSamplingStrategy(num_iter, burnin) {
  }
  bool SampleWithKindsEveryIter::sample_gov_kind(int i, int d, int e, int m, bool latent_kinds) {
    return latent_kinds && true;
  }
  bool SampleWithKindsEveryIter::sample_rel_kind(int i, int d, int e, int m, bool latent_kinds) {
    return latent_kinds && true;
  }
  bool SampleWithKindsEveryIter::sample_slot(int i, int d, int e) {
    return true;
  }
  bool SampleWithKindsEveryIter::reestimate_gov_kind(int iteration) {
    return !heldout_ && (1+iteration) >= burn_in && (1+iteration) % this->reestimate_gov_kind_every() == 0;
  }
  bool SampleWithKindsEveryIter::reestimate_rel_kind(int iteration) {
    return !heldout_ && (1+iteration) >= burn_in && (1+iteration) % this->reestimate_rel_kind_every() == 0;
  }
}
