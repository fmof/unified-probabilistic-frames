#ifndef FERRUM_OPTIMIZE_TCC_
#define FERRUM_OPTIMIZE_TCC_

#include "ferrum/optimize.hpp"
#include <memory>

namespace optimize {

  // template <typename... Args>
  // std::shared_ptr<Function>
  // Minimizer::get_function
  // (
  //  OptimizationMethod opt_how,
  //  Args... fn_args
  //  ) {
  //   std::shared_ptr<Function> func(NULL);
  //   if(opt_how == optimize::OptimizationMethod::LBFGS) {
  //     MinimizationFunctionCreator<optimize::OptimizationMethod::LBFGS> mfc;
  //     func = mfc(fn_args...);
  //     // std::shared_ptr<optimize::LibLBFGSFunction>
  //     // 	(new optimize::LibLBFGSFunction(fn_args...));
  //   } else if(opt_how == optimize::OptimizationMethod::ADAGRAD) {
  //     MinimizationFunctionCreator<optimize::OptimizationMethod::ADAGRAD> mfc;
  //     func = mfc(fn_args...);
  // 	// std::shared_ptr<optimize::AdaGradFunction>
  // 	// (new optimize::AdaGradFunction(fn_args...));
  //   } else {
  //     ERROR << "Cannot create function for unknown optimization method " << opt_how;
  //     throw 10;
  //   }
  //   return func;
  // }

  // template <typename... Args>
  // std::shared_ptr<Function>
  // Minimizer::get_function
  // (
  //  OptimizationMethod opt_how,
  //  Args... fn_args
  //  ) {
  //   std::shared_ptr<Function> func(NULL);
  //   if(opt_how == optimize::OptimizationMethod::LBFGS) {
  //     MinimizationFunctionCreator<optimize::OptimizationMethod::LBFGS> mfc;
  //     func = mfc(fn_args...);
  //     // std::shared_ptr<optimize::LibLBFGSFunction>
  //     // 	(new optimize::LibLBFGSFunction(fn_args...));
  //   } else if(opt_how == optimize::OptimizationMethod::ADAGRAD) {
  //     MinimizationFunctionCreator<optimize::OptimizationMethod::ADAGRAD> mfc;
  //     func = mfc(fn_args...);
  // 	// std::shared_ptr<optimize::AdaGradFunction>
  // 	// (new optimize::AdaGradFunction(fn_args...));
  //   } else {
  //     ERROR << "Cannot create function for unknown optimization method " << opt_how;
  //     throw 10;
  //   }
  //   return func;
  // }

#ifndef PUBLIC_UPF
  template <>
  template <typename... Args>
  std::shared_ptr<Function>
  MinimizationFunctionCreator<optimize::OptimizationMethod::LBFGS>::operator()
  (   
   Args... fn_args
   ) {
    std::shared_ptr<Function> func(new optimize::LibLBFGSFunction(fn_args...));
    return func;
  }
#endif
  template <>
  template <typename... Args>
  std::shared_ptr<Function>
  MinimizationFunctionCreator<optimize::OptimizationMethod::ADAGRAD>::operator()
  (   
   Args... fn_args
   ) {
    std::shared_ptr<Function> func(new optimize::AdaGradFunction(fn_args...));
    return func;
  }
}

#endif
