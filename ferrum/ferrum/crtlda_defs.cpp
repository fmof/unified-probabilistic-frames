#include "ferrum/crtlda.hpp"
#include "ferrum/data_util.hpp" // for RedisThriftSmartWriter
#include "ferrum/concrete.hpp" // for ConcreteSmartWriter
#include "ferrum/thrift_protocol_defs.hpp"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

//#include "ferrum/dmc.hpp"
#include "ferrum/util.hpp"

namespace ferrum {
  /////////////////////
  // DISCRETE MODEL  //
  /////////////////////
  const int* const DiscreteModel::make_slot_arr(const int n_slots){
    const int* const arr = (int*)ferrum::MALLOC(sizeof(int)*n_slots);
    return arr;
  }

  void DiscreteModel::fill_hyper(double* &array, double value, const int size) {
    ferrum::allocate_1d(array, size);
    ferrum::fill_array1d<double>(array, value, size);
  }

  SymmetricHyperparams::SymmetricHyperparams(double base) : 
    h_theta(base), h_slot(base), h_gov(base), h_rel(base), h_gov_kind(base), h_rel_kind(base) {
  }
  SymmetricHyperparams::SymmetricHyperparams() : 
    SymmetricHyperparams(0.1) {
  }
  //////////////////////////////////////////////////////////////////////////////
  
  BaseSituationLabeler::BaseSituationLabeler(bool label) :
    BaseSituationLabeler(label, 1) {
  }
  BaseSituationLabeler::BaseSituationLabeler(bool label, int label_every) :
    sw(NULL),
    every(label_every),
    do_labeling(label) {
  }
  BaseSituationLabeler::~BaseSituationLabeler() {
    if(sw != NULL) {
      delete sw;
    }
  }
  ferrum::SmartWriter* BaseSituationLabeler::operator()() {
    return sw;
  }
  bool BaseSituationLabeler::perform_labeling() {
    return do_labeling;
  }
  bool BaseSituationLabeler::perform_labeling(int iter) {
    return perform_labeling() && (iter % every == 0);
  }
  template <typename P, template <typename> class TSW >
  ConcreteSituationLabeler<P, TSW>::ConcreteSituationLabeler(bool label, int levery) :
    BaseSituationLabeler(label, levery) {
  }
  template <typename P, template <typename> class TSW >
  ConcreteSituationLabeler<P, TSW>::ConcreteSituationLabeler(bool label) :
    ConcreteSituationLabeler(label, 1) {
  }
  template <typename P, template <typename> class TSW >
  ConcreteSituationLabeler<P, TSW>::~ConcreteSituationLabeler() {
  }
  template <typename P, template <typename> class TSW >
  void ConcreteSituationLabeler<P, TSW>::make(const std::string& base) {
    if(sw == NULL) {
      //sw = new concrete::util::ConcreteSmartWriter<P>(base);
      sw = new TSW<P>(base);
    } else {
      ERROR << "Cannot recreate the CSW (yet)";
    }
  }
  template <typename P, template <typename> class TSW >
  TSW<P>* ConcreteSituationLabeler<P, TSW>::operator()() {
    if(sw == NULL) {
      return NULL;
    }
    return dynamic_cast< TSW<P>* >(sw);
  }
  template <typename P, template <typename> class TSW >
  bool ConcreteSituationLabeler<P, TSW>::perform_labeling() {
    return BaseSituationLabeler::perform_labeling();
  }
  template <typename P, template <typename> class TSW >
  bool ConcreteSituationLabeler<P, TSW>::perform_labeling(int iter) {
    return BaseSituationLabeler::perform_labeling(iter);
  }
}

template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TBinaryProtocol, ::concrete::util::ConcreteSmartWriter >;
template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TCompactProtocol, ::concrete::util::ConcreteSmartWriter >;
template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TJSONProtocol, ::concrete::util::ConcreteSmartWriter >;

template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TBinaryProtocol, ::ferrum::db::RedisThriftSmartWriter>;
template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TCompactProtocol, ::ferrum::db::RedisThriftSmartWriter>;
template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TJSONProtocol, ::ferrum::db::RedisThriftSmartWriter>;

template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TBinaryProtocol, ::ferrum::TarThriftSmartWriter>;
template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TCompactProtocol, ::ferrum::TarThriftSmartWriter>;
template class ferrum::ConcreteSituationLabeler<ferrum::thrift::TJSONProtocol, ::ferrum::TarThriftSmartWriter>;
