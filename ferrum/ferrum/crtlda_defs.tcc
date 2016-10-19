#ifndef FERRUM_LIBNAR_CRTLDA_DEFS_TCC_
#define FERRUM_LIBNAR_CRTLDA_DEFS_TCC_

#include "ferrum/crtlda_defs.hpp"
#include "ferrum/minsky.hpp"

namespace ferrum {
  template <typename DocType>
  inline int get_num_mentions(const DocType& comm) {
    return 0;
  }

  template <typename P, template <typename> class TSW>
  template <typename... Args>
  TSW<P>* ConcreteSituationLabeler<P, TSW>::make_with_args(Args... args) {
    if(sw == NULL) {
      sw = new TSW<P>(args...);
    } else {
      ERROR << "Cannot recreate the CSW (yet)";
      sw = NULL;
    }
    return dynamic_cast<TSW<P>*>(sw);
  }
}

#endif
