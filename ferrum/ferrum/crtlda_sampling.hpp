#ifndef FERRUM_LIBNAR_CRTLDA_SAMPLING_H_
#define FERRUM_LIBNAR_CRTLDA_SAMPLING_H_

#include "ferrum/crtlda_defs.hpp"
#include "ferrum/minsky.hpp"
#include "ferrum/dmc.hpp"
#include "ferrum/util.hpp"
#include "ferrum/concrete.hpp"
#include "ferrum/crtlda_writers.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>

#include <fstream>
#include <iostream>
#include "stdlib.h"
#include <time.h>

// for pair
#include <map>
#include <utility>
#include <unordered_set>
#include <string>
#include <vector>

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

namespace ferrum {
  class SamplingStrategy {
  public:
    int num_iterations;
    int burn_in;
    int reestimate_template_usage_every_;
    int reestimate_slot_usage_every_;
    int reestimate_gov_every_;
    int reestimate_rel_every_;
    virtual void heldout(bool b);
    virtual bool heldout();
    virtual bool sample_template(int, int, int) = 0;
    virtual bool sample_slot(int, int, int) = 0;
    virtual bool reestimate_hyperparameters(int iteration);
    virtual int reestimate_template_usage_every();
    virtual int reestimate_slot_usage_every();
    virtual int reestimate_gov_every();
    virtual int reestimate_rel_every();
    virtual bool reestimate_template_usage(int iteration) = 0;
    virtual bool reestimate_gov(int iteration) = 0;
    virtual bool reestimate_rel(int iteration) = 0;
    virtual bool reestimate_slot_usage(int iteration) = 0;
    SamplingStrategy(int num_iter, int burnin);
    virtual ~SamplingStrategy();
  protected:
    bool heldout_;
  };
  class WithKindsSamplingStrategy : virtual public SamplingStrategy {
  public:
    virtual int reestimate_gov_kind_every();
    virtual int reestimate_rel_kind_every();
    virtual bool sample_gov_kind(int iteration, int d_index, int e_index, int m_index, bool latent_kinds) = 0;
    virtual bool sample_rel_kind(int iteration, int d_index, int e_index, int m_index, bool latent_kinds) = 0;
    virtual bool reestimate_gov_kind(int iteration) = 0;
    virtual bool reestimate_rel_kind(int iteration) = 0;
    WithKindsSamplingStrategy(int num_iter, int burnin);
    ~WithKindsSamplingStrategy();
    int reestimate_gov_kind_every_;
    int reestimate_rel_kind_every_;
  };
  class SampleEveryIter : virtual public SamplingStrategy {
  public:
    SampleEveryIter(int num_iter, int burnin);
    virtual bool sample_template(int i, int d, int e);
    virtual bool sample_slot(int i, int d, int e);
    virtual bool reestimate_template_usage(int iter_index);
    virtual bool reestimate_slot_usage(int iter_index);
    virtual bool reestimate_gov(int iter_index);
    virtual bool reestimate_rel(int iter_index);
  };
  class SampleWithKindsEveryIter : public SampleEveryIter, public WithKindsSamplingStrategy {
  public:
    SampleWithKindsEveryIter(int num_iter, int burnin);
    virtual bool sample_gov_kind(int i, int d, int e, int m, bool latent_kinds);
    virtual bool sample_rel_kind(int i, int d, int e, int m, bool latent_kinds);
    virtual bool sample_slot(int i, int d, int e);
    virtual bool reestimate_gov_kind(int iteration);
    virtual bool reestimate_rel_kind(int iteration);
  };
}

#endif // FERRUM_LIBNAR_CRTLDA_SAMPLING_H_
