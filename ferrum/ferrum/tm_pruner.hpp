#ifndef FERRUM_LIBNAR_TM_PRUNER_HPP_
#define FERRUM_LIBNAR_TM_PRUNER_HPP_

#include "ferrum/concrete.hpp"

namespace ferrum {
  template <typename M, typename C>
  class BOWPruner {
  public:
    //protected:
    const concrete::Communication& communication;
  public:
    BOWPruner<M, C>(const concrete::Communication& comm) : communication(comm) {
    }
    virtual M clause_create(const C&) const = 0;
    typedef C concrete_type;
  };
}

#endif
