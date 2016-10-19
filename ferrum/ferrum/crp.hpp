// Header for CRP-derived distributions

#ifndef FERRUM_LIBNAR_CRP_H_
#define FERRUM_LIBNAR_CRP_H_

#include "ferrum/dmc.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/mathops.hpp"

#include <gsl/gsl_sf_log.h>
#include <stdexcept>
#include <unordered_map>


namespace crp {
  
  /**
   * Class for maintaining occupancy histograms across tables.
   * In addition to maintaining the histogram, this maintains
   * the number of tables and the total number of customers
   * associated with these tables.
   */
  class OccupancyHistogram {
  public:
    typedef std::unordered_map<int, int> histogram_type;
  private:
    histogram_type histogram_;
    int tables_;
    int customers_;
  public:
    OccupancyHistogram() : tables_(0), customers_(0) {
    }
    /**
     * Return the number of customers. This is equivalent (but much faster)
     * to calling:
     *
     *   int total = 0;
     *   for(auto val_pair : this->histogram()) {
     *     total += val_pair.first * val_pair.second;
     *   }
     */
    inline const int customers() const {
      return customers_;
    }
    inline void add_customer() {
      ++customers_;
    }
    inline void remove_customer() {
      --customers_;
    }
    /**
     * Return the total number of tables. This is equivalent (but much faster)
     * to summing all values in this histogram:
     *
     *   int total = 0;
     *   for(auto val_pair : this->histogram()) {
     *     total += val_pair.second;
     *   }
     */
    inline const int tables() const {
      return tables_;
    }
    inline void create_table() {
      ++tables_;
      if(histogram_.find(1) == histogram_.end()) {
	histogram_[1] = 1;
      } else {
	++histogram_[1];
      }
    }

    /**
     * Add a customer to an existing table.
     */
    inline void seat_customer_at_existing(histogram_type::iterator itp) {
      ++histogram_[itp->first + 1];
      if((itp->second -= 1) == 0){
	histogram_.erase(itp->first);
      }
    }
    /**
     * Return the number of deleted tables: should be either 0 or 1
     */
    inline int unseat_customer(histogram_type::iterator itp) {
      int num_del = 0;
      // always decrease the current count
      --histogram_[ itp->first ];
      // now, if it's the case that, at any other table, there's another customer
      // transfer the count down one
      if(itp->first > 1) {
	++histogram_[ itp->first - 1 ];
      } else {
	// otherwise, we have to delete the table
	--tables_;
	++num_del;
      }
      // now, we need to adjust the frequencies
      if((itp->second -= 1) == 0) {
	histogram_.erase(itp->first);
      }
      return num_del;
    }
    inline void erase(histogram_type::iterator it) {
      histogram_.erase(it);
    }
    inline int& operator[](int table_occupancy) {
      return histogram_[table_occupancy];
    }
    inline const histogram_type& histogram() {
      return histogram_;
    }
    inline histogram_type::iterator begin() { 
      return histogram_.begin(); 
    }
    inline histogram_type::iterator end() { 
      return histogram_.end();
    }
  };

  /**
   * This class defines a Chinese Restaurant Process over
   * observable items of type D (default type: int). 
   */
  template <typename D = int, class Hash = std::hash<D> >
  class CRP {
  protected:
    typedef std::unordered_map<D, OccupancyHistogram, Hash> histo_type;
    // For every dish, keep a (sparse) histogram of table occupancies.
    histo_type table_hist_by_dish_;
    // For every dish, record the number of customers eating that dish
    std::unordered_map<D, int> customers_by_dish_;
    // the total number of customers being served
    int total_customers_;
    // the total number of tables
    int num_tables_;

    double concentration_;
    double log_concentration_;
    double discount_;
    double log_discount_;
    
  public:
    typedef typename histo_type::const_iterator const_iterator;

    CRP<D, Hash>(double concentration, double discount) : concentration_(concentration), discount_(discount) {
      log_concentration_ = gsl_sf_log(concentration_);
      if(discount == 0.0) {
	log_discount_ = mathops::NEGATIVE_INFINITY;
      } else if(discount >= 1.0) {
	ERROR << "Cannot create a CRP with discount >= 1.0";
	throw std::domain_error("Cannot create a CRP with discount >= 1.0");
      } else {
	log_discount_ = gsl_sf_log(discount_);
      }
    }

    inline const double concentration() {
      return concentration_;
    }
    inline const double discount() {
      return discount_;
    }
    inline const_iterator begin() { 
      return table_hist_by_dish_.begin(); 
    }
    inline const_iterator end() {
      return table_hist_by_dish_.end();
    }
    /**
     * Return the total number of occupied tables.
     */
    inline const int num_tables() {
      return num_tables_;
    }

    /**
     * Return the total number of tables serving a particular dish.
     */
    inline const int num_tables(const D& dish) {
      typename histo_type::const_iterator it = table_hist_by_dish_.find(dish);
      return it != table_hist_by_dish_.end() ? it->second.tables() : 0;
    }

    /**
     * Return the number of customers eating the "dish" type.
     */
    inline const int count(const D& dish) {
      typename histo_type::const_iterator it = table_hist_by_dish_.find(dish);
      return (it != table_hist_by_dish_.end()) ? (it->second.customers()) : 0;
    }

    /**
     * Return whether a dish label has been seen or not.
     */
    inline int novel_dish(const D& dish) {
      return table_hist_by_dish_.find(dish) == table_hist_by_dish_.end();
    }

    /**
     * Return the probability of a dish-type label, with prior 
     * probability p(label) = p0. The probability of label is:
     *
     *   if label is new:
     *        p0 * [ concentration + discount * num_tables ]
     *        ----------------------------------------------
     *               num_customers + concentration
     *   otherwise, if label == l:
     *        count(label) - discount * num_tables_serving(label) 
     *        --------------------------------------------------- +
     *                    num_customers + concentration
     *
     *           p0 * [ concentration + discount * num_tables ]
     *           ----------------------------------------------
     *                    num_customers + concentration
     */
    inline const double prob(const D& label, double p0) {
      double base_num = concentration_ + discount_ * num_tables_;
      double num = base_num * p0;
      if(!novel_dish(label)) {
	// add to the num
	num += count(label) - discount_ * num_tables(label);
      }
      double denom = total_customers_ + concentration_;
      return num / denom;
    }

    /**
     * Given a prior probabilty p(label) = p0, this returns the
     * probability of sharing an existing table:
     *
     */
    inline const double prob_share_table(const D& label, int table_index, double p0) {
      typename histo_type::const_iterator it = table_hist_by_dish_.find(label);
      if(it == table_hist_by_dish_.end()) {
	return 0.0;
      }
      double num = it->second.tables() - discount_;
      double denom = it->second.customers() + concentration_;
      return num / denom;
    }
    
    /**
     * Given a prior probabilty p(label) = p0, this returns the
     * probability of sharing an existing table. It is equivalent, but faster, to
     *
     *    sum_{t : Table} prob_share_table(label, t, p0)
     *
     */
    inline const double prob_share_any_table(const D& label, const double p0) {
      typename histo_type::const_iterator it = table_hist_by_dish_.find(label);
      if(it == table_hist_by_dish_.end()) {
	return 0.0;
      }
      const double num = it->second.customers() - discount_*num_tables_;
      const double denom = it->second.customers() + concentration_ * p0;
      return num / denom;
    }

    inline const double prob_new_table(const D& dish, const double p0) {
      //return 1.0 - prob_share_any_table(dish);
      typename histo_type::const_iterator it = table_hist_by_dish_.find(dish);
      if(it == table_hist_by_dish_.end()) {
	return 1.0;
      }
      const double num = (concentration_ + discount_*num_tables_) * p0;
      const double denom = it->second.customers() + concentration_ * p0;
      return num / denom;
    }

  protected:
    void create_table(OccupancyHistogram* ohp) {
      ++num_tables_;
      ohp->create_table();
    }

    // void remove_table(OccupancyHistogram* ohp) {
    // 	???;
    // 	--num_tables_;
    // }

    inline void add_customer(OccupancyHistogram* ohp) {
      ohp->add_customer();
      ++total_customers_;
    }

    inline void remove_customer(OccupancyHistogram* ohp) {
      ohp->remove_customer();
      --total_customers_;
    }

  public:
    void seat(const D& label, const double p0) {
      INFO << "Seating " << label;
      const double p_share = prob_share_any_table(label, p0);
      const double p_new = prob_new_table(label, p0);
      const int sampled_table = dmc::cat::u_sample(p_new, p_share);
      OccupancyHistogram& oh = table_hist_by_dish_[label];
      if(sampled_table == 0) {
	// create a new table, and assign!
	create_table(&oh);
      } else {
	// sample from the histogram
	double r = mathops::sample_uniform(count(label));
	for(std::unordered_map<int,int>::iterator hist_iter = oh.begin();
	    hist_iter != oh.end(); ++hist_iter) {
	  r -= hist_iter->first * hist_iter->second;
	  if(r <= 0) {
	    oh.seat_customer_at_existing(hist_iter);
	    // ++oh[hist_iter->first + 1];
	    // if((hist_iter->second -= 1) == 0){
	    // 	// delete the table
	    // 	oh.erase(hist_iter);
	    // }	      
	    break;
	  }
	}
      }
      add_customer(&oh);
    }
    void unseat(D label) {
      typename histo_type::iterator it = table_hist_by_dish_.find(label);
      if(it == table_hist_by_dish_.end()){
	ERROR << "Trying to unseat a customer of type " << label << ", but our table occupany histogram is empty!";
	return;
      }
      OccupancyHistogram& oh = it->second;
      double r = mathops::sample_uniform(oh.customers());
      //for(auto oh_record : oh) {
      for(std::unordered_map<int,int>::iterator oh_record = oh.begin();
	  oh_record != oh.end(); ++oh_record) {
	r -= oh_record->first * oh_record->second;
	if(r <= 0.0) {
	  num_tables_ -= oh.unseat_customer(oh_record);
	}
      }
      remove_customer(&oh);
    }
      
    // template <typename O> inline D sample(const O& observation, double p0) {
      
    // }

    // def _sample_table(self, k):
    //   if k not in self.tables: return -1
    //   p_new = (self.theta + self.d * self.ntables) * self.base.prob(k)
    //     norm = p_new + self.ncustomers[k] - self.d * len(self.tables[k])
    //     x = random.random() * norm
    //     for i, c in enumerate(self.tables[k]):
    // 	  if x < c - self.d: return i
    //         x -= c - self.d
    //     return -1

    /////////////////////////////////////////////////////////////////

    // inline const double log_prob(D dish, double log_p0) {
    //   return log_prob(dish, 0.0, 0.0, log_p0);
    // }

    // //TODO: fix
    // inline const double log_prob(D dish, double del_numer, double del_denom, double log_p0) {
    //   // FIXME: cache log(alpha)
    //   double num =  mathops::log_add(log(count(dish) + del_numer), log(concentration_) + log_p0);
    //   double denom = log(_total_customers + del_denom + _alpha);
    //   return num - denom;
    // }
  };

  /**
   * This class defines a Pitman-Yor Process over
   * observable items of type D (default type: int). 
   */
  template <typename D = int, class Hash = std::hash<D> >
  class PYP : public CRP<D, Hash> {
  public:
  };

  // class IntegerCRP : public CRP<int> {
  // };

}

#endif
