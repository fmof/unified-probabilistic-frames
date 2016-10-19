#ifndef FERRUM_LIBNAR_CRTLDA_VARIATIONAL_H_
#define FERRUM_LIBNAR_CRTLDA_VARIATIONAL_H_

#include "ferrum/concrete.hpp"
#include "ferrum/crtlda_defs.hpp"
#include "ferrum/dmc.hpp"
#include "ferrum/mathops.hpp"
#include "ferrum/svi_util.hpp"
#include "ferrum/util.hpp"
#include "ferrum/crtlda_writers.hpp"

//#include "ferrum/minsky.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <fstream>
#include <gsl/gsl_rng.h>
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
  struct UniformHyperSeedWeightedInitializer {
  private:
    std::vector<double> from_vec(const std::vector<double>& hyper, int mult) {
      const double d = (double)mult/(double)hyper.size();
      std::vector<double> vec(hyper.size(), d);
      ferrum::sum_in_first(&vec, hyper);
      mathops::add_uniform_noise(&vec, -d, d);
      return vec;
    }
  public:
    UniformHyperSeedWeightedInitializer() {
    }
    std::vector<double> assignment(int support_size) {
      double inv_nt_ = 1.0/(double)support_size;
      std::vector<double> vec(support_size, inv_nt_);
      mathops::add_uniform_noise(&vec, -inv_nt_, inv_nt_);
      double norm = ferrum::sum(vec);
      ferrum::scalar_product(1.0/norm, &vec);
      return vec;
    }
    std::vector<double> slots(const std::vector<double>& hyper) {
      return from_vec(hyper, 1);
    }
    std::vector<double> frames(const std::vector<double>& hyper) {
      return from_vec(hyper, 1);
    }
    std::vector<double> verbs(const std::vector<double>& hyper) {
      return from_vec(hyper, 1);
    }
    std::vector<double> roles(const std::vector<double>& hyper) {
      return from_vec(hyper, 1);
    }
    std::vector<double> arcs(const std::vector<double>& hyper) {
      return from_vec(hyper, 1);
    }
    std::vector<double> words(const std::vector<double>& hyper) {
      return from_vec(hyper, 1);
    }
    std::vector<double> usage_template(const std::vector<double>& hyper, int num_words) {
      return from_vec(hyper, num_words);
    }
  };
    
  struct VStrategy {
    VStrategy() {
    }
    double em_frobenius_threshold = 1E-6;
    double eta_density_threshold = 1E-4;
    StepSizeUpdater ssu;
    int batch_size = -1;
    int num_learn_iters = 100;
    int num_e_iters = 25;
    int num_m_iters = 1;
    int num_e_threads = 1;
    int num_m_threads = 1;
    int hyper_update_min = 20;
    int hyper_update_iter = 5;
    int update_model_every = 5;
    int partial_restarts = 0;
    int num_learn_restart_iters = 25;
    int num_e_restart_iters = 25;
    int print_topics_every = 5;
    int print_topics_k = 10;
    int print_usage_every = 5;
    int em_verbosity = 1;
    int label_every = 5;
    bool heldout = false;
    bool never_update_model = false;
    bool compute_elbo = true;
    bool elbo_as_break = true;
  };
}

#endif
