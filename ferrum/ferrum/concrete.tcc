#ifndef FERRUM_LIBNAR_CONCRETE_TCC_
#define FERRUM_LIBNAR_CONCRETE_TCC_

#include <array>
#include <tuple>

namespace concrete {
  namespace util {
    template <typename C, typename... ITs>
    StructIterator<C,ITs...>::ConcreteStructIteratorReturn::ConcreteStructIteratorReturn(StructIterator<C,ITs...>* p) :
      obj_ptr_(NULL),
      par_(p) {
      need_to_reset_.fill(true);
      init_iters_();
    }

    template <typename C, typename... ITs>
    template <typename It>
    StructIterator<C,ITs...>::IterStruct<It>::IterStruct() :
      c{}, e{} {
    }      
      
    template <typename C, typename... ITs>
    typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn::pointer
    StructIterator<C, ITs...>::ConcreteStructIteratorReturn::operator->() {
      return obj_ptr_;
    }

    template <typename C, typename... ITs>
    typename StructIterator<C, ITs...>::ConcreteStructIteratorReturn::reference
    StructIterator<C, ITs...>::ConcreteStructIteratorReturn::operator*() {
      return *obj_ptr_;
    }

    // template <typename C, typename... ITs>
    // typename StructIterator<C, ITs...>::ConcreteStructIteratorReturn&
    // StructIterator<C, ITs...>::ConcreteStructIteratorReturn::operator++() {
    //   return *obj_ptr_;
    // }
 
    // template <typename C, typename... ITs>   
    // bool operator==(const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i1,
    // 		    const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i2) {
    //   return i1.obj_ptr_ == i2.obj_ptr_;
    // }
    // template <typename C, typename... ITs>
    // bool operator!=(const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i1,
    // 		    const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i2) {
    //   return i1.obj_ptr_ != i2.obj_ptr_;
    // }
    template <typename C, typename... ITs>   
    bool StructIterator<C,ITs...>::ConcreteStructIteratorReturn::operator==(const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i1) const {
      return obj_ptr_ == i1.obj_ptr_;
    }
    template <typename C, typename... ITs>
    bool StructIterator<C,ITs...>::ConcreteStructIteratorReturn::operator!=(const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i1) const {
      return i1.obj_ptr_ != obj_ptr_;
    }

    template <typename C, typename... ITs>
    typename StructIterator<C,ITs...>::iterator
    StructIterator<C, ITs...>::begin() {
      iterator iter(this);
      return iter;
    }
    template <typename C, typename... ITs>
    typename StructIterator<C, ITs...>::iterator
    StructIterator<C, ITs...>::end() {
      return iterator::end();
    }

    template <typename C, typename... ITs>
    typename StructIterator<C,ITs...>::const_iterator
    StructIterator<C, ITs...>::begin() const {
      iterator iter(this);
      return iter;
    }
    template <typename C, typename... ITs>
    typename StructIterator<C, ITs...>::const_iterator
    StructIterator<C, ITs...>::end() const {
      return iterator::end();
    }

    template <typename C, typename... ITs>
    StructIterator<C, ITs...>::ConcreteStructIteratorReturn::ConcreteStructIteratorReturn() :
      obj_ptr_(NULL),
      par_(NULL) {
    }

    template <typename C, typename... ITs>
    typename StructIterator<C, ITs...>::ConcreteStructIteratorReturn
    StructIterator<C, ITs...>::ConcreteStructIteratorReturn::end() {
      return ConcreteStructIteratorReturn();
    }

    template <typename C, typename... ITs>
    StructIterator<C, ITs...>::StructIterator(const concrete::Communication& comm) : c_(&comm) {
    }

    template <typename C, typename... IterTypes>
    template <int which>
    int StructIterator<C, IterTypes...>::ConcreteStructIteratorReturn::check(bool check_current) {
      auto& curr_iter = std::get<which>(iters_);
      int rec_ret = which;
      if(check_current) {
	bool curr_at_end = (curr_iter.c == curr_iter.e);
	if(curr_at_end) {
	  goto check_reset_part;
	} else {
	  return rec_ret;
	}
      } else {
	bool next_at_end = (++(curr_iter.c) == curr_iter.e);
	if(next_at_end) {
	  goto check_reset_part;
	} else {
	  return rec_ret;
	}
      }
    check_reset_part:
      obj_ptr_ = NULL;
      need_to_reset_[which] = true;
      if(which == 0) {
	return which-1;
      }
      // add recursive base case
      /*int o*/rec_ret = check< (which > 0 ? which - 1 : which) >(check_current);
      //if(orec_ret != (which - 1)) rec_ret = orec_ret;
      return rec_ret;
    }

    template <typename C, typename... IterTypes>
    typename StructIterator<C, IterTypes...>::ConcreteStructIteratorReturn&
    StructIterator<C, IterTypes...>::ConcreteStructIteratorReturn::advance_(bool init, bool check_current) {
      constexpr int num_axes = (int)(sizeof...(IterTypes) - 1);
      if( !init && std::get<num_axes>(iters_).c == std::get<num_axes>(iters_).e ) {
      	obj_ptr_ = NULL;
      	goto ret;
      }
      {
	int end_idx = check< num_axes >(check_current);
	if(end_idx < 0) {
	  obj_ptr_ = NULL;
	  goto ret;
	}
	if(end_idx < num_axes) {
	  reset_inner_(index< num_axes >());
	} else {
	  auto& curr = std::get<num_axes>(iters_);
	  if(curr.c != curr.e) {
	    set_obj_ptr_();
	  } else {
	    ERROR << "Did not expect this";
	  }
	}
      }
    ret:
      return *this;
    }

    /**
     * Advance the iterator forward; guaranteed to point to
     * a (de)referencable address, or `::end()`.
     */
    template <typename C, typename... IterTypes>
    typename StructIterator<C, IterTypes...>::ConcreteStructIteratorReturn&
    StructIterator<C, IterTypes...>::ConcreteStructIteratorReturn::operator++() {
      if(par_ == NULL) {
	ERROR << "par_ is null!!!";
	throw 10;
      }
      return advance_(false, false);
    //   constexpr int num_axes = (int)(sizeof...(IterTypes) - 1);
    //   if( std::get<num_axes>(iters_).c == std::get<num_axes>(iters_).e ) {
    // 	obj_ptr_ = NULL;
    // 	goto ret;
    //   }
    //   {
    // 	int end_idx = check< num_axes >();
    // 	if(end_idx < 0) {
    // 	  obj_ptr_ = NULL;
    // 	  goto ret;
    // 	}
    // 	if(end_idx < num_axes) {
    // 	  reset_inner_(index< num_axes >());
    // 	}
    //   }
    // ret:
    //   return *this;
    }
  }
}
#endif
