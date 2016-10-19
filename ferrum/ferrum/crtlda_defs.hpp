#ifndef FERRUM_LIBNAR_CRTLDA_DEFS_H_
#define FERRUM_LIBNAR_CRTLDA_DEFS_H_

#include "ferrum/dmc.hpp"
#include "ferrum/entities.hpp"
#include "ferrum/minsky.hpp"
#include "ferrum/lock.hpp"
#include "ferrum/util.hpp"
#include "ferrum/concrete.hpp"
#include "ferrum/words.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>

#include <fstream>
#include <iostream>
#include "stdlib.h"
#include <time.h>

#include <map>
#include <memory>
#include <mutex> // for std::lock_guard
#include <utility>
#include <unordered_set>
#include <string>
#include <thread>
#include <vector>

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

namespace ferrum {
  template <typename DocType>
  inline
  int get_num_mentions(const DocType& comm);

  template <typename G, typename R>
  class Mention {
  private:
    AnnotatedToken< G > gov_;
    AnnotatedToken< G > head_;
    R rel_;
    std::string rel_str_;
    std::vector< std::string > span_;
    std::string id_;
    bool is_latent_;
    friend std::ostream & operator<<(std::ostream &os, const Mention<G,R>& mention) {
      return os << "Mention[gov:" << mention.rel_ << "<" << mention.rel_str_ << ">(gov = " << mention.gov_ << ", head = " << mention.head_ <<  ")]";
    }
    
  public:
    Mention< G, R >() {      
    }
    void latent(bool b) {
      is_latent_ = b;
    }
    bool latent() const {
      return is_latent_;
    }
    void id(const std::string& id) {
      id_ = id;
    }
    const std::string& id() {
      return id_;
    }
    void gov(const AnnotatedToken<G>& gov) {
      gov_ = gov;
    }
    void head(const AnnotatedToken<G>& head) {
      head_ = head;
    }
    void rel(const R&& rel) {
      rel_ = rel;
    }
    void rel(const R& rel) {
      rel_ = rel;
    }
    void rel_str(const std::string&& rel_str) {
      rel_str_ = rel_str;
    }
    void rel_str(const std::string& rel_str) {
      rel_str_ = rel_str;
    }
    void span(const std::vector<std::string>& span) {
      span_ = span;
    }
    const AnnotatedToken<G>& gov() const {
      return gov_;
    }
    const AnnotatedToken<G>& head() const {
      return head_;
    }
    const R& rel() const {
      return rel_;
    }
    const std::string& rel_str() const {
      return rel_str_;
    }
    const std::vector< std::string >& span() const {
      return span_;
    }
  };

  template <typename G, typename R, typename C>
  class Entity {
  private:
    C canonical_name_;
    std::vector< Mention< G, R > > mentions_;
    std::string id_;
  public:
    typedef G predicate_t;
    typedef R relation_t;
    typedef C canon_t;
    const std::vector< Mention< G, R > >& mentions() {
      return mentions_;
    }
    const std::vector< Mention< G, R > >& mentions() const {
      return mentions_;
    }
    friend const std::vector<Mention<G,R> >& get_mentions(const Entity<G,R,C>& entity) {
      return entity.mentions_;
    }
    friend int num_mentions(const Entity<G,R,C>& entity) {
      return entity.mentions_.size();
    }
    // inline ferrum::db::RedisQuery form_redis_query() const {
      
    // }
    inline const std::string& id() const {
      return id_;
    }
    inline void id(const std::string& _id_) {
      id_ = _id_;
    }
    inline const std::string& redis_id() {
      return  "entity::" + id_;
    }

    inline C canonical_name() {
      return canonical_name_;
    }
    inline void canonical_name(C name) {
      canonical_name_ = name;
    }

    inline const int num_mentions() const {
      return mentions_.size();
    }
    inline const int num_mentions() {
      return const_cast<const Entity<G,R,C>* >(this)->num_mentions();
    }
    inline const Mention<G, R>& mention(int mention_index) const {
      return mentions_[mention_index];
    }
    inline const Mention<G, R>& operator[](std::size_t idx) const {
      return mentions_[idx];
    }
    
    inline void add_mentions(const std::vector< Mention< G, R> >& mentions) {
      mentions_.insert(mentions_.end(), mentions.begin(), mentions.end());
    }
    inline void add_mention(const Mention<G,R>& mention) {
      mentions_.push_back(mention);
    }

    inline ferrum::pair_icount<int> gov_view_lemma_histogram(const Vocabulary<G>& vocab,
							      const Vocabulary<std::string>& lemma_vocab) const {
      ferrum::pair_icount<int> hist;
      for(const Mention<G, R>& ment : mentions_) {
	if(ment.latent()) continue;
	int v = vocab.index( ment.gov().view() );
	int l = lemma_vocab.index(ment.gov().lemma());
	hist[ std::pair<int,int>(v, l) ] += 1;
      }
      return hist;
    }
    inline ferrum::pair_icount<int> rel_view_str_histogram(const Vocabulary<R>& vocab,
							    const Vocabulary<std::string>& lemma_vocab) const {
      ferrum::pair_icount<int> hist;
      for(const Mention<G, R>& ment : mentions_) {
	if(ment.latent()) continue;
	int v = vocab.index( ment.rel() );
	int l = lemma_vocab.index(ment.rel_str());
	hist[ std::pair<int,int>(v, l) ] += 1;
      }
      return hist;
    }
    inline std::map<int, int> gov_lemma_histogram(const Vocabulary<std::string>& vocab) const {
      std::map<int,int> hist;
      for(const auto& ment : mentions_) {
	hist[vocab.index( ment.gov().lemma() )] += 1;
      }
      return hist;
    }
    inline std::map<int, int> rel_str_histogram(const Vocabulary<std::string>& vocab) const {
      std::map<int,int> hist;
      for(const auto& ment : mentions_) {
	int i = vocab.index(ment.rel_str());
	hist[i] += 1;
      }
      return hist;
    }
    inline std::map<int, int> gov_histogram(const Vocabulary<G>& vocab) const {
      std::map<int,int> hist;
      for(const Mention<G, R>& ment : mentions_) {
	if(ment.latent()) continue;
	hist[vocab.index( ment.gov().view() )] += 1;
      }
      return hist;
    }
    inline std::map<int, int> rel_histogram(const Vocabulary<R>& vocab) const {
      std::map<int,int> hist;
      for(const Mention<G, R>& ment : mentions_) {
	if(ment.latent()) continue;
	int i = vocab.index(ment.rel());
	hist[i] += 1;
      }
      return hist;
    }
    inline std::map<int, int> gov_histogram(const Vocabulary<G>& vocab) {
      return const_cast<const Entity<G, R, C>* >(this)->gov_histogram(vocab);
    }
    inline std::map<int, int> rel_histogram(const Vocabulary<R>& vocab) {
      return const_cast<const Entity<G, R, C>* >(this)->rel_histogram(vocab);
    }

    inline std::map<G, int> gov_histogram() const {
      std::map<G,int> hist;
      for(const Mention<G, R>& ment : mentions_) {
	if(ment.latent()) continue;
	hist[ment.gov().view()] += 1;
      }
      return hist;
    }
    inline std::map<R, int> rel_histogram() const {
      std::map<R,int> hist;
      for(const Mention<G, R>& ment : mentions_) {
	if(ment.latent()) continue;
	hist[ment.rel()] += 1;
      }
      return hist;
    }
    std::map< std::pair<G, std::string> , int> compute_gov_view_lemma_counts() const {
      std::map< std::pair<G, std::string>, int> counts;
      for(const auto& mention : this->mentions()) {
	++counts[std::pair<G, std::string>(mention.gov().view(), 
					   mention.gov().lemma())];
      }
      return counts;
    }
    std::map< std::pair<R, std::string> , int> compute_rel_arc_counts() const {
      std::map< std::pair<R, std::string>, int> counts;
      for(const auto& mention : this->mentions()) {
	++counts[std::pair<G, std::string>(mention.rel(), 
					   mention.rel_str())];
      }
      return counts;
    }
  };

  // an entity counter for ferrum::Entity<G,R,C>
  template <typename G, typename R, typename C>
  class GRCEntityCounter : public VirtualEntityCounterCRTP< Entity<G, R, C> > {
  protected:
    const void* inner_pred_vocab(const VirtualVocabulary* vv) const {
      return reinterpret_cast<void* >(vv->downcast<G>() );
    }
    const void* inner_rel_vocab(const VirtualVocabulary* vv) const {
      return reinterpret_cast<void* >(vv->downcast<R>() );
    }
  public:
    typedef Entity<G,R,C> E;
    bool use_views;
    GRCEntityCounter<G,R,C>(bool uv) : use_views(uv) {
    }
    virtual ~GRCEntityCounter() {
    }
    void count_predicates(const VirtualVocabulary* vv, const E& entity, std::vector<double>& counts) const {
      const Vocabulary< G >* vocab = vv->downcast<G>();
      if(use_views) {
	for(const auto& ment : entity.mentions()) {
	  if(ment.latent()) continue;
	  counts[vocab->index( ment.gov().view() )] += 1;
	}
      } else {
	for(const auto& ment : entity.mentions()) {
	  if(ment.latent()) continue;
	  counts[vocab->index( ment.gov().lemma() )] += 1;
	}
      }
    }
    void count_relations(const VirtualVocabulary* vv, const E& entity, std::vector<double>& counts) const {
      const Vocabulary< R >* vocab = vv->downcast<R>();
      if(use_views) {
	for(const auto& ment : entity.mentions()) {
	  if(ment.latent()) continue;
	  counts[vocab->index( ment.rel() )] += 1;
	}
      } else {
	for(const auto& ment : entity.mentions()) {
	  if(ment.latent()) continue;
	  counts[vocab->index( ment.rel_str() )] += 1;
	}
      }
    }
  };

  template <typename G, typename R, typename C>
  class DocumentGRC {
  private:
    friend std::ostream & operator<<(std::ostream &os, const DocumentGRC<G,R,C>& doc) {
      return os << "Document[id:" << doc.id <<  "]";
    }
    std::vector< Entity<G, R, C> > entities_;
    typedef Mention<G,R> MentionGR;
  public:
    std::string id;
    DocumentGRC<G,R,C>() : id() {
    }
    DocumentGRC<G,R,C>(const std::string id_) : id(id_) {
    }
    inline const std::string& redis_id() {
      return  "document::" + id;
    }
    // void to_redis(ferrum::db::Redis& db) {
    //   std::stringstream stream;
    //   stream << std::to_string(entities_.size());
    //   for(const auto& ent : entities_) {
    // 	stream << ent.redis_id() << " ";
    //   }
    //   ferrum::db::RedisQuery rq
    // 	(
    // 	 redis_id(),
    // 	 "entities",
    // 	 doc_stream.str()
    // 	 )
    // 	.hset();
    //   for(const auto& ent : entities_) {
    // 	rq += ent.form_redis_query().hset();
    //   }
    //   rq.query_all(db);
    // }

    template < typename Pruner > // e.g., of type EntityMentionPruner< MentionGR >
    DocumentGRC<G,R,C>(const concrete::Communication& communication,
  		       const Pruner& mention_pruner) : id(communication.id) {
      const concrete::EntitySet* const esp = concrete_util::first_entity_set(communication,
  									     "Stanford");
      if(! esp) return;
      int num_mentions = 0;
      std::vector< concrete::Entity > c_entity_set = esp->entityList;
      concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention =
  	concrete_util::mention_id_to_mention(communication, "Stanford");
      for(std::vector< concrete::Entity >::const_iterator entity_iterator = c_entity_set.begin();
  	  entity_iterator != c_entity_set.end(); ++entity_iterator) {
  	concrete::Entity conc_entity = *entity_iterator;
  	Entity<G, R, C> entity;
  	entity.canonical_name(conc_entity.canonicalName);
  	entity.id(conc_entity.uuid.uuidString);
  	for(concrete::UUID mention_uuid : conc_entity.mentionIdList) {
  	  const concrete::EntityMention conc_entity_mention = mention_id_to_mention[mention_uuid];
  	  // create a mention based on the mention pruner
  	  std::vector<Mention<G,R> > created_mentions = mention_pruner.prune(conc_entity_mention);
  	  entity.add_mentions(created_mentions);
  	}
  	if(entity.num_mentions() > 0) {
  	  num_mentions += entity.num_mentions();
  	  this->add_entity(entity);
  	}
      }
      INFO << "Document " << this->id << " has " << num_entities() << " entities and " << num_mentions << " mentions";
    }

    const std::vector< Entity<G, R, C> >& entities() const {
      return entities_;
    }    
    const std::vector< Entity<G, R, C> >& entities() {
      return entities_;
    }
    friend const std::vector<Entity<G,R,C> >& get_entities(const DocumentGRC<G,R,C >& doc) {
      return doc.entities_;
    }
    friend int num_entities(const DocumentGRC<G,R,C >& doc) {
      return doc.num_entities();
    }
    friend int get_num_mentions(const DocumentGRC<G,R,C >& doc) {
      int nm = 0;
      for(const auto& ent : doc.entities()) {
	nm += ent.mentions().size();
      }
      return nm;
    }

    inline const int num_entities() const {
      return entities_.size();
    }
    inline const int num_entities() {
      return entities_.size();
    }

    inline void add_entity(const Entity<G, R, C>& entity) {
      entities_.push_back(entity);
    }

    inline const Entity<G, R, C>& operator[](const size_t& idx) const {
      return entities_[idx];
    }

    inline const G& governor(int ei, int mi) const {
      return entities_[ei].mention(mi).gov().view();
    }
    inline const R& relation(int ei, int mi) const {
      return entities_[ei].mention(mi).relation();
    }
    inline const C& canonical_name(int ei) const {
      return entities_[ei].canonical_name();
    }

    std::map<std::string, int> compute_gov_lemma_counts() const {
      std::map<std::string, int> counts;
      const int num_ent = this->num_entities();
      for(int ei = 0; ei < num_ent; ++ei) {
  	const auto& entity = entities_[ei];
  	const int num_ments = entity.num_mentions();
	for(int mi = 0; mi < num_ments; ++mi) {
	  ++counts[entity.mention(mi).gov().lemma()];
	}
      }
      return counts;
    }

    std::map<std::string, int> compute_gov_counts() const {
      std::map<std::string, int> counts;
      const int num_ent = this->num_entities();
      for(int ei = 0; ei < num_ent; ++ei) {
	const auto& entity = entities_[ei];
	const int num_ments = entity.num_mentions();
	for(int mi = 0; mi < num_ments; ++mi) {
	  ++counts[entity.mention(mi).gov().lemma()];
	}
      }
      return counts;
    }

    std::map<G, int> compute_gov_view_counts() const {
      std::map<G, int> counts;
      const int num_ent = this->num_entities();
      for(int ei = 0; ei < num_ent; ++ei) {
	const auto& entity = entities_[ei];
	const int num_ments = entity.num_mentions();
	for(int mi = 0; mi < num_ments; ++mi) {
	  ++counts[entity.mention(mi).gov().view()];
	}
      }
      return counts;
    }

    std::map< std::pair<G, std::string> , int> compute_gov_view_lemma_counts() const {
      std::map< std::pair<G, std::string>, int> counts;
      for(const auto& entity : entities_) {
	for(const auto& mention : entity.mentions()) {
	  ++counts[std::pair<G, std::string>(mention.gov().view(), 
					     mention.gov().lemma())];
	}
      }
      return counts;
    }

    std::map<R, int> compute_rel_counts() const {
      std::map<R, int> counts;
      const int num_ent = this->num_entities();
      for(int ei = 0; ei < num_ent; ++ei) {
	const auto& entity = entities_[ei];
	const int num_ments = entity.num_mentions();
	for(int mi = 0; mi < num_ments; ++mi) {
	  const R rel = entity.mention(mi).rel();	  
	  ++counts[rel];
	}
      }
      return counts;
    }

    std::map< std::pair<R, std::string> , int> compute_rel_arc_counts() const {
      std::map< std::pair<R, std::string>, int> counts;
      for(const auto& entity : entities_) {
	for(const auto& mention : entity.mentions()) {
	  ++counts[std::pair<G, std::string>(mention.gov().rel(), 
					     mention.gov().rel_str())];
	}
      }
      return counts;
    }


  };

  class BaseCorpusIterator {
  public:
    virtual ~BaseCorpusIterator() {
    }
  };

  template <typename VCI>
  class VirtualCorpusIterator : public BaseCorpusIterator {
  };

  // forward declare
  template <typename DocType>
  class InMemoryCorpusIterator;
  template <typename D /*, template <typename> class CIt = CorpusIterator*/ >
  class InMemoryCorpus;

  template <typename DocType>
  class InMemoryCorpusIterator : public VirtualCorpusIterator< InMemoryCorpusIterator<DocType> > {
  private:
    const std::vector<DocType>* documents_; // this is shared
    typename std::vector<DocType>::const_iterator it_;
    InMemoryCorpusIterator<DocType>() : document(), iteration_idx(0) {
    }
  public:
    InMemoryCorpusIterator<DocType>(const std::vector<DocType>* docs) :
      documents_(docs),
      it_(docs->begin()),
      iteration_idx(0) {
      document = const_cast<DocType*>(&*it_);
    }
    InMemoryCorpusIterator<DocType>(const InMemoryCorpus<DocType>* corp) : InMemoryCorpusIterator<DocType>(&(corp->documents)) {
    }
    static InMemoryCorpusIterator<DocType> end(const std::vector<DocType>* docs) {
      InMemoryCorpusIterator<DocType> ci;
      ci.documents_ = docs;
      ci.it_ = ci.documents_->end();
      ci.iteration_idx = ci.documents_->size();
      return ci;
    }
    static InMemoryCorpusIterator<DocType> end(const InMemoryCorpus<DocType>* corpus) {
      return InMemoryCorpusIterator<DocType>::end(&(corpus->documents));
    }
    typedef const InMemoryCorpusIterator<DocType>& reference;
    typedef const InMemoryCorpusIterator<DocType>* pointer;

    DocType* document;
    unsigned int iteration_idx;
    
    pointer operator->() const {
      return this;
    }
    reference operator*() const {
      return *this;
    }
    InMemoryCorpusIterator<DocType>& operator++() { //prefix
      ++it_;
      document = const_cast<DocType*>(&(*it_));
      ++iteration_idx;
      return *this;
    }
    friend bool operator==(const InMemoryCorpusIterator<DocType>& i1, const InMemoryCorpusIterator<DocType>& i2) {
      return i1.it_ == i2.it_;
    }
    friend bool operator!=(const InMemoryCorpusIterator<DocType>& i1, const InMemoryCorpusIterator<DocType>& i2) {
      return i1.it_ != i2.it_;
    }
  };

  struct BackgroundComputeParams {
  public:
    double min_freq = 0.0001;
    double min_log_prob = -200.0;
  };

  class AbstractCorpus {
  public:
    virtual ~AbstractCorpus() {
    }
    virtual const std::string& get_name() = 0;
    virtual const int num_docs() = 0;
    virtual const int num_docs() const = 0;

    //virtual const AbstractCorpus* downcast() const = 0;
    //virtual int num_mentions() const = 0;
  };

  template <typename DocType, template <typename> class CType, template <typename> class IteratorT >
  class CRTPCorpus : public AbstractCorpus {
  protected:
    std::string name_;
    ferrum::Mutex mut_;
    bool mt_;
  public:
    typedef IteratorT<DocType> const_iterator;
    CRTPCorpus<DocType, CType, IteratorT>(const std::string& name) : name_(name), mut_(), mt_() {
    }
    CRTPCorpus<DocType, CType, IteratorT>() : CRTPCorpus<DocType, CType, IteratorT>("") {
    }
    virtual ~CRTPCorpus() {
    }
    void multithreaded(bool b) {
      mt_ = b;
    }
    bool multithreaded() {
      return mt_;
    }
    const std::string& get_name() const {
      return name_;
    }
    const std::string& get_name() {
      //return const_cast< CRTPCorpus<DocType, CType, IteratorT>* >(this)->get_name();
      return name_;
    }
    void name(const std::string& n_name) {
      name_ = n_name;
    }

    // virtual CType<DocType>* base_clone() const {
    //   std::lock_guard<ferrum::Mutex> lock(mut_);
    //   CType<DocType>* clone = new CType<DocType>;
      
    //   return clone;
    // }

    virtual void copy_crucial_attributes_into(CType<DocType>* receiver) const {
      receiver->mt_ = mt_;
    }
    
    // synchronized
    CType<DocType>* create(const std::string& n_name) {
      std::lock_guard<ferrum::Mutex> lock(mut_);
      CType<DocType>* corpus = new CType<DocType>;
      this->copy_crucial_attributes_into(corpus);
      corpus->reset();
      corpus->name(n_name);
      return corpus;
    }

    virtual const DocType& operator[](const size_t idx) = 0;
    virtual void fill(const size_t idx, DocType* doc) = 0;

    // void __atomic_populate(CType<DocType>* ncorp, const std::vector<size_t>& chosen) {
    //   for(size_t d_idx : chosen) {
    // 	DocType doc;
    // 	// the following is synchronized
    // 	fill(d_idx, &doc);
    // 	ncorp->add_document( doc );
    //   }
    // }

    void populate_from(CType<DocType>* ncorp, const std::vector<size_t>& chosen) {
      for(size_t d_idx : chosen) {
	DocType doc;
	// the following is synchronized
	fill(d_idx, &doc);
	ncorp->add_document( doc );
      }
      // if(mt_) {
      // 	//INFO << "Attempting to acquire lock...";
      // 	//std::lock_guard<ferrum::Mutex> lock(mut_);
      // 	//INFO << "lock acquired";
      // 	this->__atomic_populate(ncorp, chosen);
      // 	//INFO << "Releasing lock";
      // } else {
      // 	this->__atomic_populate(ncorp, chosen);
      // }
    }

    virtual void reset() {
    }

    virtual CType<DocType>* random_subset(size_t num_to_pick, unsigned int iseed = 0) {
      std::vector<size_t> chosen(num_to_pick);
      size_t nd = this->num_docs();
      // the following do/while is for scoping memory management
      typedef unsigned long int ul;
      ul rand_num = 0;
      std::hash<std::thread::id> hasher;
      ul seed = hasher(std::this_thread::get_id());
      seed = (seed & 0xFFFF) | ((ul)iseed << 16);
      INFO << "Using seed " << seed;
      do {
	std::vector<size_t> all(nd);
	for(size_t i = 0; i < nd; ++i) {
	  all[i] = i;
	}
	gsl_rng *rg = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rg, seed);
	rand_num = gsl_rng_get(rg);
	mathops::random_choose(chosen.data(),
			       num_to_pick,
			       all.data(), nd,
			       sizeof(size_t), rg);
	gsl_rng_free(rg);
      } while(0);
      
      std::string r_name(name_ + ":random" + std::to_string(num_to_pick) + "::"+ std::to_string(rand_num) + "::" + std::to_string(seed));
      //CType<DocType>* ret_c = NULL;
      // if(mt_) {
      // 	INFO << "Attempting to acquire lock...";
      // 	std::lock_guard<ferrum::Mutex> lock(mut_);
      // 	INFO << "lock acquired";

      // 	ret_c = create( r_name );
      // 	populate_from(ret_c, chosen);
      // 	INFO << "Releasing lock";
      // } else {
      CType<DocType>* ret_c = create( r_name );
      populate_from(ret_c, chosen);
	//}
      return ret_c;
    }

    const CType<DocType>* downcast() const {
      return static_cast<const CType<DocType>*>(this);
    }

    virtual const_iterator begin() const {
      const_iterator ci( downcast());
      return ci;
    }
    virtual const_iterator end() const {
      return const_iterator::end( downcast() );
    }
    virtual void add_document(const DocType& document) = 0;
    virtual bool has_document(const std::string&) {
      return false;
    }
  };

  template <typename DocType, template <typename> class CType, template <typename> class IteratorT >
  inline int num_mentions(CRTPCorpus<DocType, CType, IteratorT>* corpus) {
    int nm = 0;
    for(IteratorT<DocType> it = corpus->begin(); it != corpus->end(); ++it) {
      const DocType& doc = *(it->document);
      nm += get_num_mentions(doc);
    }
    return nm;
  }

  template <typename DocType, template <typename> class CType, template <typename> class IteratorT >
  inline int num_words(CRTPCorpus<DocType, CType, IteratorT>* corpus) {
    int nm = 0;
    for(IteratorT<DocType> it = corpus->begin(); it != corpus->end(); ++it) {
      const DocType& doc = *(it->document);
      nm += num_words(doc);
    }
    return nm;
  }

  template <typename CCorpus, typename Vocab>
  class BackgroundComputer {
  public:
    virtual ~BackgroundComputer() {
    }
    void complete_background_compute(std::vector<double>* background, BackgroundComputeParams bcp) {
      // make sure there are no zero counts
      ferrum::ensure_min(bcp.min_freq, background);
      // divide by the total number of words in the corpus
      double num_words = ferrum::sum(*background);
      ferrum::scalar_product(1.0/num_words, background);
      // now log it all
      ferrum::log(background);
      // truncate to min_log_prob
      ferrum::ensure_min(bcp.min_log_prob, background);
      // renormalize
      dmc::cat::log_renormalize(background);
    }
    virtual void before() {
    }
    virtual void after(const Vocab*, CCorpus*, const std::vector<double>*, const std::string&) {
    }
  };

  /**
   * An iterable, and batch-iterable, collection.
   */
  template <typename D>
  class InMemoryCorpus : public CRTPCorpus< D, InMemoryCorpus, InMemoryCorpusIterator > {
  protected:
    std::vector<D> documents;
  public:
    friend class InMemoryCorpusIterator<D>;
    typedef D DocType;
    typedef InMemoryCorpusIterator<DocType> const_iterator;
    template <typename DT> using CItPartial = InMemoryCorpusIterator<DT>;
    InMemoryCorpus<D>(const std::string& corpus_name) : CRTPCorpus< D, ::ferrum::InMemoryCorpus, ::ferrum::InMemoryCorpusIterator >(corpus_name) {
    }
    InMemoryCorpus<D>() : InMemoryCorpus<D>("") {
    }
    InMemoryCorpus<D>(const CRTPCorpus< D, ::ferrum::InMemoryCorpus, ::ferrum::InMemoryCorpusIterator >* corpus) : InMemoryCorpus<D>( *(corpus->downcast()) ) {
    }

    virtual ~InMemoryCorpus() {
    }

    virtual void reset() {
      documents.clear();
    }
    
    // virtual const InMemoryCorpus<D>* downcast() const {
    //   return static_cast<const InMemoryCorpus<D>*>(this);
    // }

    // virtual const_iterator begin() const {
    //   const_iterator ci(this);
    //   return ci;
    // }
    // virtual const_iterator end() const {
    //   return const_iterator::end(this);
    // }
    virtual void add_document(const D& document) {
      documents.push_back(document);
    }
    virtual int num_mentions() const {
      int nm = 0;
      for(const auto& it : *this) {
	const DocType& doc = *(it.document);
	nm += get_num_mentions(doc);
      }
      return nm;
    }

    virtual const std::vector<D>& get_corpus() const {
      return documents;
    }
    virtual const D& operator[](const size_t idx) {
      return documents[idx];
    }
    virtual void fill(const size_t idx, DocType* doc) {
      *doc = documents[idx];
    }
    // virtual const ??? operator[](const std::pair<size_t, size_t>& endpoints) const {
    // }
    virtual const int num_docs() {
      return documents.size();
    }
    virtual const int num_docs() const {
      return documents.size();
    }

    virtual InMemoryCorpus<D> subset(size_t start_i, size_t end_i) const {
      InMemoryCorpus<D> sub(this);
      sub.reset();
      sub.documents =
	std::vector<D>
	(
	 documents.begin() + start_i,
	 documents.begin() + ( end_i > documents.size() ? documents.size() : end_i )
	 );
      return sub;
    }

    std::vector<int> vec_entity_count() {
      const int nd = num_docs();
      std::vector<int> e(nd);
      for(int i = 0; i<nd; i++) {
	e[i] = documents[i].num_entities();
      }
      return e;
    }

    /**
     * Return a list of stop words for mention govs.
     * This is computed by IDF; the top **percentage** are returned.
     */
    template <typename W>
    inline std::vector<std::pair<W, double> > gov_stopwords_by_idf(double percentage) {
      std::map<W, int> doc_occurs;
      for(const D& doc : documents) {
	std::map<W, int> doc_verb_idf = doc.compute_gov_counts();
	for(const auto& entry : doc_verb_idf) {
	  doc_occurs[entry.first] += 1;
	}
      }
      int v_size = doc_occurs.size();
      std::vector< std::pair<W, double> > idf(v_size);
      const double N = (double)(this->num_docs());
      int i = 0;
      for(const auto& entry : doc_occurs) {
	idf[i++] = std::pair<W, double>(entry.first, gsl_sf_log(N / (entry.second + 1.0)));
      }
      // sort it in ascending mode
      std::sort(idf.begin(), idf.end(),
		[](const std::pair<W,double> &left, const std::pair<W,double> &right) {
		  return left.second < right.second;
		});
      // then figure out how many to prune off
      int num_to_return = v_size * percentage;
      typename std::vector< std::pair<W, double> >::const_iterator first = idf.begin();
      typename std::vector< std::pair<W, double> >::const_iterator last = idf.begin() + num_to_return;
      std::vector< std::pair<W, double> > newVec(first, last);
      return newVec;
    }

    template <typename W>
    inline std::map< std::pair<int, int>, int> gov_view_doc_cooccur(const Vocabulary<W>& vocab) const {
      std::map< std::pair<int, int>, int> doc_occurs;
      for(const D& doc : documents) {
	std::map<W, int> doc_verb_idf = doc.compute_gov_view_counts();
	std::vector<int> seen;
	for(const auto& entry : doc_verb_idf) {
	  seen.push_back(vocab.index(entry.first));
	}
	const int num_types_seen = seen.size();
	for(int i = 0; i < num_types_seen; ++i) {
	  for(int j = 0; j < i; ++j) {
	    const int w_i = seen[i];
	    const int w_j = seen[j];
	    doc_occurs[ std::pair<int, int>(w_i, w_j) ] += 1;
	    //doc_occurs[ std::pair<int, int>(w_j, w_i) ] += 1;
	  }
	}
      }
      return doc_occurs;
    }
    template <typename W>
    inline std::map<int, int> gov_view_doc_occur(const Vocabulary<W>& vocab) const {
      std::map<int, int> doc_occurs;
      for(const D& doc : documents) {
	std::map<W, int> doc_verb_idf = doc.compute_gov_view_counts();
	for(const auto& entry : doc_verb_idf) {
	  doc_occurs[ vocab.index(entry.first) ] += 1;
	}
      }
      return doc_occurs;
    }

    inline std::map< std::pair<int, int>, int> gov_lemma_doc_cooccur(const Vocabulary<std::string>& vocab) const {
      std::map< std::pair<int, int>, int> doc_occurs;
      for(const D& doc : documents) {
	std::map<std::string, int> doc_verb_idf = doc.compute_gov_lemma_counts();
	std::vector<int> seen;
	for(const auto& entry : doc_verb_idf) {
	  seen.push_back(vocab.index(entry.first));
	}
	const int num_types_seen = seen.size();
	for(int i = 0; i < num_types_seen; ++i) {
	  const int w_i = seen[i];
	  for(int j = 0; j < i; ++j) {
	    const int w_j = seen[j];
	    std::pair<int,int> p(w_i, w_j);
	    int prev = doc_occurs[p];
	    doc_occurs[ p ] = prev + 1;
	    //doc_occurs[ std::pair<int, int>(w_j, w_i) ] += 1;
	  }
	}
      }
      return doc_occurs;
    }
    inline std::map<int, int> gov_lemma_doc_occur(const Vocabulary<std::string>& vocab) const {
      std::map<int, int> doc_occurs;
      for(const D& doc : documents) {
	std::map<std::string, int> doc_verb_idf = doc.compute_gov_lemma_counts();
	for(const auto& entry : doc_verb_idf) {
	  doc_occurs[ vocab.index(entry.first) ] += 1;
	}
      }
      return doc_occurs;
    }

}; // ends InMemoryCorpus

  // TODO: implement this as callers to doc.${appropriate histogram function}
  // template <typename G, typename R, typename C>
  // std::map<int, int> doc_multinomial(const ferrum::DocumentGRC<G,R,C>& doc,
  // 				     const ferrum::AnnotationLevel& level,
  // 				     const ferrum::StructureType& structure) {
    
  // };

  class SymmetricHyperparams {
  public:
    double h_theta;
    double h_slot;
    double h_gov;
    double h_rel;
    double h_gov_kind;
    double h_rel_kind;
    SymmetricHyperparams();
    SymmetricHyperparams(double symmetric_base);
  };

  /**
   * A Discrete Chain-Restricted Template admixture model.
   * This is derived from LDA, where observations are
   * generated by discrete distributions with Dirichlet priors.
   * The basic assumptions are that a document is a
   * collection of entities (corefence chains), and each
   * entity is a collection of individual mentions.
   * Both templates and slots are assigned at the
   * entity level, and some combination of these assignments
   * is responsible for every mention (or entity) observable.
   * A mention, or even an entity, may have any number
   * of discrete observations associated with it: it is
   * up to the inference algorithm to properly handle
   * those observations.
   */
  class DiscreteModel {
  private:
    // setting model parameters
    int num_templates_;
    // this array is of size num_templates:
    // how many slots for each template
    std::vector<int> num_slots_;

    //// hyperparameters
    // how to use each template
    std::vector<double> hyper_theta_;
    std::vector<double> hyper_slot_;
    // how often governors show
    std::vector<double> hyper_gov_;
    bool h_gov_set_;
    // how often relations show
    std::vector<double> hyper_rel_;
    bool h_rel_set_;

    //// priors
    // Each document has a prior probability of using
    // any particular template
    std::vector< std::vector<double> > prior_template_;
    // Each template has a prior probability of using
    // any particular slot
    std::vector< std::vector<double> > prior_slot_;
    // Each template has a prior probability of generating 
    // any particular governor
    std::vector< std::vector<double> > prior_governor_;
    // Each template-specific slot has a prior probability 
    // of generating any particular relation
    std::vector< std::vector< std::vector<double> > > prior_relation_;

    const int* const make_slot_arr(const int n_slots);
    void fill_hyper(double* &array, double value, const int size);

  public:
    DiscreteModel(int n_templates, int n_slots,
		  SymmetricHyperparams* shp) : num_templates_(n_templates),
					       num_slots_(std::vector<int>(n_templates, n_slots)),
					       hyper_theta_(std::vector<double>(n_templates, shp->h_theta)), 
					       hyper_slot_(std::vector<double>(n_slots, shp->h_slot)) {
    }
    DiscreteModel(int n_templates, const std::vector<int>& n_slots) : 
      num_templates_(n_templates), num_slots_(n_slots),
      hyper_theta_(std::vector<double>(n_templates)), 
      hyper_slot_(std::vector<double>(n_slots.size())) {
    }
    ~DiscreteModel() {
    }

    int num_templates() {
      return num_templates_;
    }
    std::vector<int> num_slots() {
      return num_slots_;
    }

    void hyper_theta(const std::vector<double>& h) {
      hyper_theta_ = h;
    }
    void hyper_slot(const std::vector<double>& h) {
      hyper_slot_ = h;
    }
    void hyper_rel(const std::vector<double>& h) {
      hyper_rel_ = h;
      h_rel_set_ = true;
    }
    void hyper_gov(const std::vector<double>& h) {
      hyper_gov_ = h;
      h_gov_set_ = true;
    }

    std::vector<double> hyper_theta() {
      return hyper_theta_;
    }
    std::vector<double> hyper_slot() {
      return hyper_slot_;
    }
    std::vector<double> hyper_gov() {
      if(! h_gov_set_) {
	ERROR << "Hyperparameters for governors are not set";
      }
      return hyper_gov_;
    }
    std::vector<double> hyper_rel() {
      if(! h_rel_set_) {
	ERROR << "Hyperparameters for relations are not set";
      }
      return hyper_rel_;
    }
    DiscreteModel& hyper_gov(double num, double val) {
      hyper_gov_ = std::vector<double>(num, val);
      h_gov_set_ = true;
      return *this;
    }
    DiscreteModel& hyper_rel(double num, double val) {
      hyper_rel_ = std::vector<double>(num, val);
      h_rel_set_ = true;
      return *this;
    }

    //transfer functions 
    inline void prior_template_usage(const std::vector< std::vector< double > >& p) {
      prior_template_ = p;
    }
    inline void prior_gov(const std::vector< std::vector< double > >& p) {
      prior_governor_ = p;
    }
    inline void prior_slot_usage(const std::vector< std::vector< double > >& p) {
      prior_slot_ = p;
    }
    inline void prior_rel(const std::vector< std::vector< std::vector< double > > >& p) {
      prior_relation_ = p;
    }

    inline const std::vector< std::vector< double> >& prior_template_usage() {
      return prior_template_;
    }
    inline const std::vector< std::vector< double> >& prior_governor() {
      return prior_governor_;
    }
    inline const std::vector< std::vector< double> >& prior_slot_usage() {
      return prior_slot_;
    }
    inline const std::vector<std::vector< std::vector< double> > >& prior_relation() {
      return prior_relation_;
    }

    template <typename W> void print_governors(const int num_per, const Vocabulary<W>& vocab) {
      int topic_idx = 0;
      for(const std::vector<double> topic : prior_governor_) {
	INFO << "Template " << topic_idx;
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	for(size_t item_idx = 0; item_idx < num_per; ++item_idx) {
	  size_t which = sorted_topic[item_idx];
	  INFO << "\t" << topic[which] << "\t" << vocab.word(which);
	}
	++topic_idx;
      }
    }
    
    template <typename G, typename R, typename C> 
    double latent_marginalized_loglikelihood(const DocumentGRC<G,R,C>&  doc,
					     int doc_id, int* num_obs_,
					     const Vocabulary<G>& gov_vocab,
					     const Vocabulary<R>& rel_vocab,
					     const std::vector<std::vector<double> >& doc_template_usage,
					     bool normalized = true) {
      double ll = 0.0;
      int num_obs = 0;
      const std::vector<double>& d_t_usage = doc_template_usage[doc_id];
      const int num_entities = doc.num_entities();
      for(int ei = 0; ei < num_entities; ++ei) {
	double ll_t = mathops::NEGATIVE_INFINITY;
	std::map<int, int> floating_gov_hist_ = doc[ei].gov_histogram(gov_vocab);
	std::map<int, int> floating_rel_hist_ = doc[ei].rel_histogram(rel_vocab);
	num_obs += (doc[ei].num_mentions()*2);
	for(int ti = 0; ti < num_templates_; ++ti) {
	  double x0 = gsl_sf_log(d_t_usage[ti]);
	  for(auto ghiter : floating_gov_hist_) {
	    const int gov_index = ghiter.first;
	    const int gov_count = ghiter.second;
	    x0 += gov_count*gsl_sf_log(prior_governor_[ti][gov_index]);
	  }
	  double ll_s = mathops::NEGATIVE_INFINITY;
	  for(int si = 0; si < num_slots_[ti]; ++si) {
	    double x = gsl_sf_log(prior_slot_[ti][si]);
	    for(auto rhiter : floating_rel_hist_) {
	      const int rel_index = rhiter.first;
	      const int rel_count = rhiter.second;
	      x += rel_count*gsl_sf_log(prior_relation_[ti][si][rel_index]);
	    }
	    ll_s = mathops::log_add(ll_s, x);
	  }
	  ll_t = mathops::log_add(ll_t, x0 + ll_s);
	}
	ll += ll_t;
      }
      if(num_obs_ != NULL) {
	*num_obs_ += num_obs;
      }
      // for every chain
      return normalized ? ll/(double)num_obs : ll;
    }


    template <typename G, typename R, typename C> 
    double latent_marginalized_loglikelihood(const InMemoryCorpus<DocumentGRC<G, R,C> >&  corpus,
					     const Vocabulary<G>& gov_vocab,
					     const Vocabulary<R>& rel_vocab,
					     const std::vector<std::vector<double> >& doc_template_usage,
					     bool normalized = true) {
      double ll = 0.0;
      int doc_id = 0;
      int num_obs = 0;
      for(auto doc : corpus.get_corpus()) {
	ll += latent_marginalized_loglikelihood<G, R,C>(doc, doc_id++, &num_obs,
							gov_vocab, rel_vocab,
							doc_template_usage,
							false);
      }
      return normalized ? ll/(double)num_obs : ll;
    }

    template <typename D, typename WW>
    std::vector<double> compute_coherences(const int M, const std::string& which, 
					   const InMemoryCorpus<D>& corpus,
					   const Vocabulary<WW>& vocab) {
      if(which == "gov") {
	std::map<int, int> gk_occur = 
	  corpus.template gov_view_doc_occur< WW >(vocab);
	std::map<std::pair<int, int>, int> gk_cooc = 
	  corpus.template gov_view_doc_cooccur< WW >(vocab);
	return ferrum::compute_coherences(M, prior_governor_, gk_occur, gk_cooc);
      }
      return std::vector<double>();
    }
  };

  ////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  template <typename GK, typename RK, typename G, typename R>
  struct DiscreteKindPrinter {
    GK* gk_vocab = NULL;
    RK* rk_vocab = NULL;
    G* g_vocab = NULL;
    R* r_vocab = NULL;

    int num_per_gov_kind = 10;
    int num_per_rel_kind = 10;
    int num_per_gov = 10;
    int num_per_rel = 10;

    struct print {
      int tg = 5;
      int usage_t = 5;
    } print;
  };

  class StringDiscreteKindPrinter :
    public DiscreteKindPrinter<
    Vocabulary<std::string>,
    Vocabulary<std::string>,
    Vocabulary<std::string>,
    Vocabulary<std::string> > {
  };

  /**
   * A ferrum::SmartWriter-based way to handle concrete::Situation and
   * concrete::SituationMention situation serializations.
   */
  class BaseSituationLabeler {
  protected:
    ferrum::SmartWriter* sw;
    void* unsafe_corpus;
    int every;
    bool do_labeling;
  public:
    BaseSituationLabeler(bool dl);
    BaseSituationLabeler(bool dl, int every);
    BaseSituationLabeler(const BaseSituationLabeler& bl) = delete;
    BaseSituationLabeler(BaseSituationLabeler&& bl) = delete;
    BaseSituationLabeler& operator=(const BaseSituationLabeler& bl) = delete;
    BaseSituationLabeler& operator=(BaseSituationLabeler&& bl) = delete;
    virtual ~BaseSituationLabeler();
    /**
     * Create the basic object to handle concrete::Situation{,Mention}
     * serialization. Access it through operator().
     */
    virtual void make(const std::string& base) = 0;
    template <typename CorpusT> CorpusT corpus() {
      return reinterpret_cast<CorpusT>(unsafe_corpus);
    }
    template <typename CorpusT> void corpus(const CorpusT* corp) {
      unsafe_corpus =
	reinterpret_cast<void*>
	(
	 const_cast< CorpusT* >(corp)
	 );
    }
    virtual ferrum::SmartWriter* operator()();
    virtual bool perform_labeling();
    virtual bool perform_labeling(int iter);
  };

  template <typename P, template <typename> class TSW >
  class ConcreteSituationLabeler : public BaseSituationLabeler {
  public:
    ConcreteSituationLabeler<P, TSW>(bool dl);
    ConcreteSituationLabeler<P, TSW>(bool dl, int every);
    ConcreteSituationLabeler<P, TSW>(const ConcreteSituationLabeler& sl) = delete;
    ConcreteSituationLabeler<P, TSW>(ConcreteSituationLabeler&& sl) = delete;
    ConcreteSituationLabeler<P, TSW>& operator=(const ConcreteSituationLabeler& sl) = delete;
    ConcreteSituationLabeler<P, TSW>& operator=(ConcreteSituationLabeler&& sl) = delete;
    ~ConcreteSituationLabeler();
    virtual void make(const std::string& base);
    template <typename... Args> TSW<P>* make_with_args(Args... args);
    virtual TSW<P>* operator()();
    virtual bool perform_labeling();
    virtual bool perform_labeling(int iter);
  };

  /**
   * A Discrete Chain-Restricted Template admixture model.
   * This is derived from LDA, where observations are
   * generated by discrete distributions with Dirichlet priors.
   * The basic assumptions are that a document is a
   * collection of entities (corefence chains), and each
   * entity is a collection of individual mentions.
   * Both templates and slots are assigned at the
   * entity level, and some combination of these assignments
   * is responsible for every mention (or entity) observable.
   * A mention, or even an entity, may have any number
   * of discrete observations associated with it: it is
   * up to the inference algorithm to properly handle
   * those observations.
   */
  class DiscreteModelWithKinds {
  private:
    // setting model parameters
    int num_templates_;
    // this array is of size num_templates:
    // how many slots for each template
    std::vector<int> num_slots_;

    //// hyperparameters
    // how to use each template
    std::vector<double> hyper_theta_;
    std::vector<double> hyper_slot_;
    // how often governors show
    std::vector<double> hyper_gov_;
    bool h_gov_set_;
    // how often relations show
    std::vector<double> hyper_rel_;
    bool h_rel_set_;
    // how often governors show
    std::vector<double> hyper_gov_kind_;
    bool h_gov_kind_set_;
    // how often relations show
    std::vector<double> hyper_rel_kind_;
    bool h_rel_kind_set_;

    //// priors
    // Each document has a prior probability of using
    // any particular template
    std::vector< std::vector<double> > prior_template_;
    // Each template has a prior probability of using
    // any particular slot
    std::vector< std::vector<double> > prior_slot_;
    // Each template has a prior probability of generating 
    // any particular governor *kind*
    std::vector< std::vector<double> > prior_governor_kind_;
    // Each template-specific slot has a prior probability 
    // of generating any particular relation *kind*
    std::vector< std::vector< std::vector<double> > > prior_relation_kind_;

    // Each governor kind has a prior probability of generating 
    // any particular lexical governor
    std::vector< std::vector<double> > prior_governor_;
    // Each relation kind has a prior probability of generating
    // any particular lexical relation
    std::vector< std::vector<double> > prior_relation_;

  public:
    DiscreteModelWithKinds(int n_templates, int n_slots) : num_templates_(n_templates),
							num_slots_(std::vector<int>(n_templates, n_slots)),
							hyper_theta_(std::vector<double>(n_templates)), 
							hyper_slot_(std::vector<double>(n_slots)) {
    }
    DiscreteModelWithKinds(int n_templates, const std::vector<int>& n_slots) : num_templates_(n_templates),
							num_slots_(n_slots),
							hyper_theta_(std::vector<double>(n_templates)), 
									       hyper_slot_(std::vector<double>(n_slots.size())) {
    }
    DiscreteModelWithKinds(int n_templates, int n_slots,
			   SymmetricHyperparams* shp) : num_templates_(n_templates),
							num_slots_(std::vector<int>(n_templates, n_slots)),
							hyper_theta_(std::vector<double>(n_templates, shp->h_theta)), 
							hyper_slot_(std::vector<double>(n_slots, shp->h_slot)) {
    }
    ~DiscreteModelWithKinds() {
    }

    int num_templates() {
      return num_templates_;
    }
    std::vector<int> num_slots() {
      return num_slots_;
    }
    int num_gov_kinds() {
      if(! h_gov_kind_set_) {
	ERROR << "Hyperparameters for governor kinds are not set";
	return -1;
      }
      return hyper_gov_kind_.size();
    }
    int num_rel_kinds() {
      if(! h_rel_kind_set_) {
	ERROR << "Hyperparameters for relation kinds are not set";
	return -1;
      }
      return hyper_rel_kind_.size();
    }

    std::vector<double> hyper_theta() {
      return hyper_theta_;
    }
    std::vector<double> hyper_slot() {
      return hyper_slot_;
    }
    std::vector<double> hyper_gov_kind() {
      if(! h_gov_kind_set_) {
	ERROR << "Hyperparameters for governor kinds are not set";
      }
      return hyper_gov_kind_;
    }
    std::vector<double> hyper_rel_kind() {
      if(! h_rel_kind_set_) {
	ERROR << "Hyperparameters for relation kinds are not set";
      }
      return hyper_rel_kind_;
    }
    DiscreteModelWithKinds& hyper_gov_kind(double num, double val) {
      hyper_gov_kind_ = std::vector<double>(num, val);
      h_gov_kind_set_ = true;
      return *this;
    }
    DiscreteModelWithKinds& hyper_rel_kind(double num, double val) {
      hyper_rel_kind_ = std::vector<double>(num, val);
      h_rel_kind_set_ = true;
      return *this;
    }
    std::vector<double> hyper_gov() {
      if(! h_gov_set_) {
	ERROR << "Hyperparameters for governors are not set";
      }
      return hyper_gov_;
    }
    std::vector<double> hyper_rel() {
      if(! h_rel_set_) {
	ERROR << "Hyperparameters for relations are not set";
      }
      return hyper_rel_;
    }
    void hyper_theta(const std::vector<double>& h) {
      hyper_theta_ = h;
    }
    void hyper_slot(const std::vector<double>& h) {
      hyper_slot_ = h;
    }
    void hyper_rel(const std::vector<double>& h) {
      hyper_rel_ = h;
      h_rel_set_ = true;
    }
    void hyper_rel_kind(const std::vector<double>& h) {
      hyper_rel_kind_ = h;
      h_rel_kind_set_ = true;
    }
    void hyper_gov(const std::vector<double>& h) {
      hyper_gov_ = h;
      h_gov_set_ = true;
    }
    void hyper_gov_kind(const std::vector<double>& h) {
      hyper_gov_kind_ = h;
      h_gov_kind_set_ = true;
    }
    DiscreteModelWithKinds& hyper_gov(double num, double val) {
      hyper_gov_ = std::vector<double>(num, val);
      h_gov_set_ = true;
      return *this;
    }
    DiscreteModelWithKinds& hyper_rel(double num, double val) {
      hyper_rel_ = std::vector<double>(num, val);
      h_rel_set_ = true;
      return *this;
    }

    //transfer functions 
    inline void prior_template_usage(const std::vector< std::vector< double > >& p) {
      prior_template_ = p;
    }
    inline void prior_gov_kind(const std::vector< std::vector< double > >& p) {
      prior_governor_kind_ = p;
    }
    inline void prior_gov(const std::vector< std::vector< double > >& p) {
      prior_governor_ = p;
    }
    inline void prior_slot_usage(const std::vector< std::vector< double > >& p) {
      prior_slot_ = p;
    }
    inline void prior_rel_kind(const std::vector< std::vector< std::vector< double > > >& p) {
      prior_relation_kind_ = p;
    }
    inline void prior_rel(const std::vector< std::vector< double > >& p) {
      prior_relation_ = p;
    }

    inline const std::vector< std::vector< double> >& prior_template_usage() {
      return prior_template_;
    }
    inline const std::vector< std::vector< double> >& prior_governor_kind() {
      return prior_governor_kind_;
    }
    inline const std::vector< std::vector< double> >& prior_governor() {
      return prior_governor_;
    }
    inline const std::vector< std::vector< double> >& prior_slot_usage() {
      return prior_slot_;
    }
    inline const std::vector<std::vector< std::vector< double> > >& prior_relation_kind() {
      return prior_relation_kind_;
    }
    inline const std::vector< std::vector< double> >& prior_relation() {
      return prior_relation_;
    }

    void print_2d_distribution(const std::vector<std::vector<double> >& vec_dist,
			       std::ostream & outter = std::cout) {
      for(const auto& dist : vec_dist) {
	std::stringstream stream;
	for(const auto& p : dist) {
	  stream << p << " ";
	}
	outter << stream.str();
	outter << std::endl;
      }
    }
    

    void print_template_usage(std::ostream & outter = std::cout) {
      print_2d_distribution(prior_template_, outter);
    }
    void print_slot_usage(std::ostream & outter = std::cout) {
      print_2d_distribution(prior_slot_, outter);
    }

    template <typename W> void print_governors(const int num_per, const Vocabulary<W>& vocab) {
      int topic_idx = 0;
      for(const std::vector<double> topic : prior_governor_kind_) {
	INFO << "Template " << topic_idx;
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	for(size_t item_idx = 0; item_idx < num_per; ++item_idx) {
	  size_t which = sorted_topic[item_idx];
	  INFO << "\t" << topic[which] << "\t" << vocab.word(which);
	}
	++topic_idx;
      }
    }

    template <typename PrintingStruct>
    void print_templates(const std::string& f_name,
			 const PrintingStruct* const ps) {
      std::ofstream myfile;
      myfile.open(f_name);
      myfile << "template\tslot\tframe\trole\tgov\twhich_type\twhich_val\tprob\n";
      for(int idx_template = 0; idx_template < num_templates_; ++idx_template) {
	// print slot
	std::vector<size_t> sorted = ferrum::sort_indices(prior_slot_[idx_template], false);
	for(const auto& item_idx : sorted) {
	  myfile << idx_template << "\t-\t-\t-\t-\tslot\t" << item_idx << "\t" << prior_slot_[idx_template][item_idx] << "\n";
	}
	// print frames (gov kinds)
	{
	  const std::vector<double>& topic = prior_governor_kind_[idx_template];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), ps->num_per_gov_kind);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << idx_template << "\t-\t-\t-\t-\tframe\t" << ps->gk_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
	// print roles (rel kinds)
	{
	  const int ns = num_slots_[idx_template];
	  for(int slot = 0; slot < ns; ++slot) {
	    const std::vector<double>& topic = prior_relation_kind_[idx_template][slot];
	    std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	    const size_t max_iter = ferrum::min(sorted.size(), ps->num_per_rel_kind);
	    for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	      size_t which = sorted[item_idx];
	      myfile << idx_template << "\t" << slot << "\t-\t-\t-\trole\t" << ps->rk_vocab->word(which) << "\t" << topic[which] << "\n";
	    }
	  }
	}
      }
      // print govs
      {
	const int num_frames = prior_governor_.size();
	for(int f_idx = 0; f_idx < num_frames; ++f_idx) {
	  const std::vector<double>& topic = prior_governor_[f_idx];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), ps->num_per_gov);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << "-\t-\t" << f_idx << "\t-\t-\tgov\t" << ps->g_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
      }
      // print rels
      {
	const int num_roles = prior_relation_.size();
	for(int r_idx = 0; r_idx < num_roles; ++r_idx) {
	  const std::vector<double>& topic = prior_relation_[r_idx];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), ps->num_per_rel);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << "-\t-\t-\t" << r_idx << "\t-\trel\t" << ps->r_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
      }      
      myfile.close();
      INFO << "Wrote templates to " << f_name;
    }
    
    template <typename G, typename R, typename C> 
    double latent_marginalized_loglikelihood(const DocumentGRC<G,R,C>&  doc,
					     int doc_id,
					     const Vocabulary<G>& gov_vocab,
					     const Vocabulary<R>& rel_vocab,
					     bool normalized = true) {
      double ll = 0.0;
      int num_obs = 0;
      WARN << "THIS FUNCTION IS NOT YET DONE!!!";
      const std::vector<double>& d_t_usage = prior_template_[doc_id];
      const int num_entities = doc.num_entities();
      for(int ei = 0; ei < num_entities; ++ei) {
	double ll_t = mathops::NEGATIVE_INFINITY;
	std::map<int, int> floating_gov_hist_ = doc[ei].gov_histogram(gov_vocab);
	std::map<int, int> floating_rel_hist_ = doc[ei].rel_histogram(rel_vocab);
	for(int ti = 0; ti < num_templates_; ++ti) {
	  double x0 = gsl_sf_log(d_t_usage[ti]);
	  for(auto ghiter : floating_gov_hist_) {
	    const int gov_index = ghiter.first;
	    const int gov_count = ghiter.second;
	    x0 += gov_count*gsl_sf_log(prior_governor_kind_[ti][gov_index]);
	    num_obs += gov_count;
	  }
	  double ll_s = mathops::NEGATIVE_INFINITY;
	  for(int si = 0; si < num_slots_[ti]; ++si) {
	    double x = gsl_sf_log(prior_slot_[ti][si]);
	    for(auto rhiter : floating_rel_hist_) {
	      const int rel_index = rhiter.first;
	      const int rel_count = rhiter.second;
	      x += rel_count*gsl_sf_log(prior_relation_kind_[ti][si][rel_index]);
	      num_obs += rel_count;
	    }
	    ll_s = mathops::log_add(ll_s, x);
	  }
	  ll_t = mathops::log_add(ll_t, x0 + ll_s);
	}
	ll += ll_t;
      }
      // for every chain
      return normalized ? ll/(double)num_obs : ll;
    }


    template <typename G, typename R, typename C> 
    double latent_marginalized_loglikelihood(const InMemoryCorpus<DocumentGRC<G, R,C> >&  corpus,
					     const Vocabulary<G>& gov_vocab,
					     const Vocabulary<R>& rel_vocab,
					     bool normalized = true) {
      double ll = 0.0;
      int doc_id = 0;
      for(auto doc : corpus.get_corpus()) {
	ll += latent_marginalized_loglikelihood<G, R,C>(doc, doc_id++,
							gov_vocab, rel_vocab,
							normalized);
      }
      return ll;
    }

    double latent_and_kind_marginalized_loglikelihood(const int num_kinds, const int obs_index,
						      const std::vector<std::vector<double> >& obs_prior,
						      const std::vector<double>& kind_prior) {
      //marginalize over the latent gov_kind
      double k_ll = mathops::NEGATIVE_INFINITY;
      for(int k = 0; k < num_kinds; ++k) {
	double lp = gsl_sf_log(obs_prior[k][obs_index]) + gsl_sf_log(kind_prior[k]);
	k_ll = mathops::log_add(k_ll, lp);
      }
      return k_ll;
    }

    template <typename E> 
    double latent_and_kind_marginalized_loglikelihood(const E&  entity,
						      const Vocabulary<std::string>& gov_vocab,
						      const Vocabulary<std::string>& rel_vocab,
						      const std::vector<double>& d_t_usage) {
      double ll_t = mathops::NEGATIVE_INFINITY;
      for(int ti = 0; ti < num_templates_; ++ti) {
	double x0 = gsl_sf_log(d_t_usage[ti]);
	double ll_m_for_gov = 0.0;
	for(int mi = 0; mi < entity.num_mentions(); ++mi) {
	  const int num_gov_kinds = hyper_gov_kind_.size();
	  const int gov_index = gov_vocab.index(entity.mention(mi).gov().lemma());
	  ll_m_for_gov += latent_and_kind_marginalized_loglikelihood(num_gov_kinds, gov_index,
								     prior_governor_, prior_governor_kind_[ti]);
	}
	x0 += ll_m_for_gov;
	double ll_s = mathops::NEGATIVE_INFINITY;
	for(int si = 0; si < num_slots_[ti]; ++si) {
	  double x = gsl_sf_log(prior_slot_[ti][si]);
	    double ll_m_for_rel = 0.0;
	    for(int mi = 0; mi < entity.num_mentions(); ++mi) {
	      const int num_rel_kinds = hyper_rel_kind_.size();
	      const int rel_index = rel_vocab.index(entity.mention(mi).rel_str());
	      ll_m_for_rel += latent_and_kind_marginalized_loglikelihood(num_rel_kinds, rel_index,
									 prior_relation_, prior_relation_kind_[ti][si]);
	    }
	    x += ll_m_for_rel;
	  ll_s = mathops::log_add(ll_s, x);
	}
	ll_t = mathops::log_add(ll_t, x0 + ll_s);
      }
      return ll_t;
    }


    template <typename G, typename R, typename C> 
    double latent_and_kind_marginalized_loglikelihood(const DocumentGRC<G,R,C>&  doc,
						      int doc_id, int* num_obs,
						      const Vocabulary<std::string>& gov_vocab,
						      const Vocabulary<std::string>& rel_vocab,
						      std::vector< std::vector<double> > corp_specific_prior_template,
						      bool normalized = true) {
      double ll = 0.0;
      int num_obs_ = 0;
      const std::vector<double>& d_t_usage = corp_specific_prior_template[doc_id];
      const int num_entities = doc.num_entities();
      for(int ei = 0; ei < num_entities; ++ei) {
	const Entity<G,R,C>& entity = doc[ei];
	double ll_t = latent_and_kind_marginalized_loglikelihood< Entity<G,R,C> >(entity, gov_vocab, rel_vocab, d_t_usage);
	num_obs_ += (entity.num_mentions()*2);
	ll += ll_t;
      }
      *num_obs += num_obs_;
      // for every chain
      return normalized ? ll/(double)num_obs_ : ll;
    }


    template <typename G, typename R, typename C> 
    double latent_and_kind_marginalized_loglikelihood(const InMemoryCorpus<DocumentGRC<G, R,C> >&  corpus,
						      const Vocabulary<std::string>& gov_vocab,
						      const Vocabulary<std::string>& rel_vocab,
						      std::vector< std::vector<double> > prior_template,
						      bool normalized = true) {
      double ll = 0.0;
      int doc_id = 0;
      int num_obs = 0;
      for(auto doc : corpus.get_corpus()) {
	ll += latent_and_kind_marginalized_loglikelihood<G, R,C>(doc, doc_id++, &num_obs,
								 gov_vocab, rel_vocab,
								 prior_template,
								 false);
      }
      return normalized ? ll/(double)num_obs : ll;
    }

    template <typename D, typename WW>
    std::vector<double> compute_coherences(const int M, const std::string& which, 
					   const InMemoryCorpus<D>& corpus,
					   const Vocabulary<WW>& vocab) {
      if(which == "gov_kind") {
	std::map<int, int> gk_occur = 
	  corpus.template gov_view_doc_occur< WW >(vocab);
	std::map<std::pair<int, int>, int> gk_cooc = 
	  corpus.template gov_view_doc_cooccur< WW >(vocab);
	return ferrum::compute_coherences(M, prior_governor_kind_, gk_occur, gk_cooc);
      }
      if(which == "gov_obs_kind_marginalized") {
	std::map<int, int> gk_occur = 
	  corpus.template gov_lemma_doc_occur(vocab);
	std::map<std::pair<int, int>, int> gk_cooc = 
	  corpus.template gov_lemma_doc_cooccur(vocab);
	//now comes the tricky part: we need to marginalize
	std::vector<std::vector< double> > marginalized;
	int n_t = prior_governor_kind_.size();
	for(int i = 0; i < n_t; ++i) {
	  std::vector<double> X;
	  const std::vector<double>& gk_i_dist = prior_governor_kind_[i];
	  int dim = gk_i_dist.size();
	  for(int obs_type = 0; obs_type < hyper_gov_.size(); ++obs_type) {
	    double p = 0.0;
	    for(int j = 0; j < dim; ++j) {
	      p += gk_i_dist[j] * prior_governor_[j][obs_type];
	    }
	    X.push_back(p);
	  }
	  marginalized.push_back(X);
	}
	return ferrum::compute_coherences(M, marginalized, gk_occur, gk_cooc);
      }
      if(which == "gov_obs") {
	std::map<int, int> gk_occur = 
	  corpus.template gov_lemma_doc_occur(vocab);
	std::map<std::pair<int, int>, int> gk_cooc = 
	  corpus.template gov_lemma_doc_cooccur(vocab);
	return ferrum::compute_coherences(M, prior_governor_, gk_occur, gk_cooc);
      }
      return std::vector<double>();
    }
  };

  /**
   * A Discrete Chain-Restricted Template admixture model.
   * This is derived from LDA, where observations are
   * generated by discrete distributions with Dirichlet priors.
   * The basic assumptions are that a document is a
   * collection of entities (corefence chains), and each
   * entity is a collection of individual mentions.
   * Both templates and slots are assigned at the
   * entity level, and some combination of these assignments
   * is responsible for every mention (or entity) observable.
   * A mention, or even an entity, may have any number
   * of discrete observations associated with it: it is
   * up to the inference algorithm to properly handle
   * those observations.
   */
  class DiscreteModelWithKindsGlobalSlots {
  private:
    // setting model parameters
    int num_templates_;
    // how many slots (in total)
    int num_slots_;

    //// hyperparameters
    // how to use each template
    std::vector<double> hyper_theta_;
    std::vector<double> hyper_slot_;
    // how often governors show
    std::vector<double> hyper_gov_;
    bool h_gov_set_;
    // how often relations show
    std::vector<double> hyper_rel_;
    bool h_rel_set_;
    // how often governors show
    std::vector<double> hyper_gov_kind_;
    bool h_gov_kind_set_;
    // how often relations show
    std::vector<double> hyper_rel_kind_;
    bool h_rel_kind_set_;

    //// priors
    // Each document has a prior probability of using
    // any particular template
    std::vector< std::vector<double> > prior_template_;
    // Each template has a prior probability of using
    // any particular slot
    std::vector< std::vector<double> > prior_slot_;
    // Each template has a prior probability of generating 
    // any particular governor *kind*
    std::vector< std::vector<double> > prior_governor_kind_;
    // Each template-specific slot has a prior probability 
    // of generating any particular relation *kind*
    std::vector< std::vector<double> > prior_relation_kind_;

    // Each governor kind has a prior probability of generating 
    // any particular lexical governor
    std::vector< std::vector<double> > prior_governor_;
    // Each relation kind has a prior probability of generating
    // any particular lexical relation
    std::vector< std::vector<double> > prior_relation_;

  public:
    DiscreteModelWithKindsGlobalSlots(int n_templates, int n_slots) : num_templates_(n_templates),
								      num_slots_(n_slots),
								      hyper_theta_(std::vector<double>(n_templates, 0.1)), 
								      hyper_slot_(std::vector<double>(n_slots, 0.1)) {
    }
    DiscreteModelWithKindsGlobalSlots
    (int n_templates, int n_slots,
     SymmetricHyperparams* shp) : num_templates_(n_templates),
				  num_slots_(n_slots),
				  hyper_theta_(std::vector<double>(n_templates, shp->h_theta)), 
				  hyper_slot_(std::vector<double>(n_slots, shp->h_slot)) {
    }
    ~DiscreteModelWithKindsGlobalSlots() {
    }

    int num_templates() {
      return num_templates_;
    }
    int num_slots() {
      return num_slots_;
    }
    int num_gov_kinds() {
      if(! h_gov_kind_set_) {
	ERROR << "Hyperparameters for governor kinds are not set";
	return -1;
      }
      return hyper_gov_kind_.size();
    }
    int num_rel_kinds() {
      if(! h_rel_kind_set_) {
	ERROR << "Hyperparameters for relation kinds are not set";
	return -1;
      }
      return hyper_rel_kind_.size();
    }

    std::vector<double> hyper_theta() {
      return hyper_theta_;
    }
    std::vector<double> hyper_usage() {
      return hyper_theta();
    }
    std::vector<double> hyper_slot() {
      return hyper_slot_;
    }
    std::vector<double> hyper_gov_kind() {
      if(! h_gov_kind_set_) {
	ERROR << "Hyperparameters for governor kinds are not set";
      }
      return hyper_gov_kind_;
    }
    std::vector<double> hyper_frame() {
      return hyper_gov_kind();
    }
    std::vector<double> hyper_rel_kind() {
      if(! h_rel_kind_set_) {
	ERROR << "Hyperparameters for relation kinds are not set";
      }
      return hyper_rel_kind_;
    }
    std::vector<double> hyper_role() {
      return hyper_rel_kind();
    }
    DiscreteModelWithKindsGlobalSlots& hyper_gov_kind(double num, double val) {
      hyper_gov_kind_ = std::vector<double>(num, val);
      h_gov_kind_set_ = true;
      return *this;
    }
    DiscreteModelWithKindsGlobalSlots& hyper_rel_kind(double num, double val) {
      hyper_rel_kind_ = std::vector<double>(num, val);
      h_rel_kind_set_ = true;
      return *this;
    }
    std::vector<double> hyper_gov() {
      if(! h_gov_set_) {
	ERROR << "Hyperparameters for governors are not set";
      }
      return hyper_gov_;
    }
    std::vector<double> hyper_verb() {
      return hyper_gov();
    }
    std::vector<double> hyper_rel() {
      if(! h_rel_set_) {
	ERROR << "Hyperparameters for relations are not set";
      }
      return hyper_rel_;
    }
    std::vector<double> hyper_arc() {
      return hyper_rel();
    }
    void hyper_theta(const std::vector<double>& h) {
      hyper_theta_ = h;
    }
    void hyper_slot(const std::vector<double>& h) {
      hyper_slot_ = h;
    }
    void hyper_rel(const std::vector<double>& h) {
      hyper_rel_ = h;
      h_rel_set_ = true;
    }
    void hyper_rel_kind(const std::vector<double>& h) {
      hyper_rel_kind_ = h;
      h_rel_kind_set_ = true;
    }
    void hyper_gov(const std::vector<double>& h) {
      hyper_gov_ = h;
      h_gov_set_ = true;
    }
    void hyper_gov_kind(const std::vector<double>& h) {
      hyper_gov_kind_ = h;
      h_gov_kind_set_ = true;
    }
    DiscreteModelWithKindsGlobalSlots& hyper_gov(double num, double val) {
      hyper_gov_ = std::vector<double>(num, val);
      h_gov_set_ = true;
      return *this;
    }
    DiscreteModelWithKindsGlobalSlots& hyper_rel(double num, double val) {
      hyper_rel_ = std::vector<double>(num, val);
      h_rel_set_ = true;
      return *this;
    }

    //transfer functions 
    inline void prior_template_usage(const std::vector< std::vector< double > >& p) {
      prior_template_ = p;
    }
    inline void prior_gov_kind(const std::vector< std::vector< double > >& p) {
      prior_governor_kind_ = p;
    }
    inline void prior_gov(const std::vector< std::vector< double > >& p) {
      prior_governor_ = p;
    }
    inline void prior_slot_usage(const std::vector< std::vector< double > >& p) {
      prior_slot_ = p;
    }
    inline void prior_rel_kind(const std::vector< std::vector< double > >& p) {
      prior_relation_kind_ = p;
    }
    inline void prior_rel(const std::vector< std::vector< double > >& p) {
      prior_relation_ = p;
    }

    inline const std::vector< std::vector< double> >& prior_template_usage() {
      return prior_template_;
    }
    inline const std::vector< std::vector< double> >& prior_governor_kind() {
      return prior_governor_kind_;
    }
    inline const std::vector< std::vector< double> >& prior_governor() {
      return prior_governor_;
    }
    inline const std::vector< std::vector< double> >& prior_slot_usage() {
      return prior_slot_;
    }
    inline const std::vector< std::vector< double> >& prior_relation_kind() {
      return prior_relation_kind_;
    }
    inline const std::vector< std::vector< double> >& prior_relation() {
      return prior_relation_;
    }

    void print_2d_distribution(const std::vector<std::vector<double> >& vec_dist,
			       std::ostream & outter = std::cout) {
      for(const auto& dist : vec_dist) {
	std::stringstream stream;
	for(const auto& p : dist) {
	  stream << p << " ";
	}
	outter << stream.str();
	outter << std::endl;
      }
    }
    void print_template_usage(std::ostream & outter = std::cout) {
      print_2d_distribution(prior_template_, outter);
    }
    void print_slot_usage(std::ostream & outter = std::cout) {
      print_2d_distribution(prior_slot_, outter);
    }

    template <typename W>
    void print_vocab_dist(const size_t num_per, const Vocabulary<W>& vocab, 
			  const std::vector< std::vector<double> >& topics) {
      int topic_idx = 0;
      for(const auto& topic : topics) {
	std::stringstream stream;
	stream << "Topic " << topic_idx << " ::: ";
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	for(size_t item_idx = 0; item_idx < num_per && item_idx < vocab.num_words(); ++item_idx) {
	  size_t which = sorted_topic[item_idx];
	  stream << vocab.word(which) << " [" << topic[which] << "] ";
	}
	INFO << stream.str();
	++topic_idx;
      }
    }
    template <typename W> void print_frames(const int num_per, const Vocabulary<W>& vocab) {
      print_vocab_dist(num_per, vocab, prior_governor_kind_);
    }
    template <typename W> void print_roles(const int num_per, const Vocabulary<W>& vocab) {
      print_vocab_dist(num_per, vocab, prior_relation_kind_);
    }
    template <typename W> void print_verbs(const int num_per, const Vocabulary<W>& vocab) {
      print_vocab_dist(num_per, vocab, prior_governor_);
    }
    template <typename W> void print_arcs(const int num_per, const Vocabulary<W>& vocab) {
      print_vocab_dist(num_per, vocab, prior_relation_);
    }
    template <typename W> void print_frames(const int num_per, const Vocabulary<W>& vocab, const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, dists);
    }
    template <typename W> void print_roles(const int num_per, const Vocabulary<W>& vocab, const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, dists);
    }
    template <typename W> void print_verbs(const int num_per, const Vocabulary<W>& vocab, const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, dists);
    }
    template <typename W> void print_arcs(const int num_per, const Vocabulary<W>& vocab, const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, dists);
    }


    template <typename W>
    void print_vocab_dist(const size_t num_per, const Vocabulary<W>& vocab,
			  std::ostream& stream,
			  const std::vector< std::vector<double> >& topics) {
      for(const auto& topic : topics) {
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	for(size_t item_idx = 0; item_idx < num_per && item_idx < vocab.num_words(); ++item_idx) {
	  size_t which = sorted_topic[item_idx];
	  stream << topic[which] << " ";
	}
	stream << "\n";
      }
    }
    template <typename W> void print_frames(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream) {
      print_vocab_dist(num_per, vocab, stream, prior_governor_kind_);
    }
    template <typename W> void print_roles(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream) {
      print_vocab_dist(num_per, vocab, stream, prior_relation_kind_);
    }
    template <typename W> void print_verbs(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream) {
      print_vocab_dist(num_per, vocab, stream, prior_governor_);
    }
    template <typename W> void print_arcs(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream) {
      print_vocab_dist(num_per, vocab, stream, prior_relation_);
    }
    template <typename W> void print_frames(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream,
					    const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, stream, dists);
    }
    template <typename W> void print_roles(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream,
					   const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, stream, dists);
    }
    template <typename W> void print_verbs(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream,
					   const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, stream, dists);
    }
    template <typename W> void print_arcs(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream,
					  const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, stream, dists);
    }

    template <typename PrintingStruct>
    void print_templates(std::ostream& myfile,
			 const PrintingStruct* const ps) {
      myfile << "template\tslot\tframe\trole\tgov\twhich_type\twhich_val\tprob\n";
      for(int idx_template = 0; idx_template < num_templates_; ++idx_template) {
	// print slot
	std::vector<size_t> sorted = ferrum::sort_indices(prior_slot_[idx_template], false);
	for(const auto& item_idx : sorted) {
	  myfile << idx_template << "\t-\t-\t-\t-\tslot\t" << item_idx << "\t" << prior_slot_[idx_template][item_idx] << "\n";
	}
	// print frames (gov kinds)
	{
	  const std::vector<double>& topic = prior_governor_kind_[idx_template];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), (unsigned long)ps->num_per_gov_kind);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << idx_template << "\t-\t-\t-\t-\tframe\t" << ps->gk_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
      }
      // print roles (rel kinds)
      {
	for(int slot = 0; slot < num_slots_; ++slot) {
	  const std::vector<double>& topic = prior_relation_kind_[slot];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), (unsigned long)ps->num_per_rel_kind);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << "-\t" << slot << "\t-\t-\t-\trole\t" << ps->rk_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
      }
      // print govs
      {
	const int num_frames = prior_governor_.size();
	for(int f_idx = 0; f_idx < num_frames; ++f_idx) {
	  const std::vector<double>& topic = prior_governor_[f_idx];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), (unsigned long)ps->num_per_gov);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << "-\t-\t" << f_idx << "\t-\t-\tgov\t" << ps->g_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
      }
      // print rels
      {
	const int num_roles = prior_relation_.size();
	for(int r_idx = 0; r_idx < num_roles; ++r_idx) {
	  const std::vector<double>& topic = prior_relation_[r_idx];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), (unsigned long)ps->num_per_rel);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << "-\t-\t-\t" << r_idx << "\t-\trel\t" << ps->r_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
      }
    }
    template <typename PrintingStruct>
    void print_templates(const std::string& f_name,
			 const PrintingStruct* const ps) {
      std::ofstream myfile;
      myfile.open(f_name);
      this->print_templates(myfile, ps);
      myfile.close();
      INFO << "Wrote templates to " << f_name;
    }
    
    template <typename G, typename R, typename C> 
    double latent_marginalized_loglikelihood(const DocumentGRC<G,R,C>&  doc,
					     int doc_id,
					     const Vocabulary<G>& gov_vocab,
					     const Vocabulary<R>& rel_vocab,
					     bool normalized = true) {
      double ll = 0.0;
      int num_obs = 0;
      WARN << "THIS FUNCTION IS NOT YET DONE!!!";
      const std::vector<double>& d_t_usage = prior_template_[doc_id];
      const int num_entities = doc.num_entities();
      for(int ei = 0; ei < num_entities; ++ei) {
	double ll_t = mathops::NEGATIVE_INFINITY;
	std::map<int, int> floating_gov_hist_ = doc[ei].gov_histogram(gov_vocab);
	std::map<int, int> floating_rel_hist_ = doc[ei].rel_histogram(rel_vocab);
	for(int ti = 0; ti < num_templates_; ++ti) {
	  double x0 = gsl_sf_log(d_t_usage[ti]);
	  for(auto ghiter : floating_gov_hist_) {
	    const int gov_index = ghiter.first;
	    const int gov_count = ghiter.second;
	    x0 += gov_count*gsl_sf_log(prior_governor_kind_[ti][gov_index]);
	    num_obs += gov_count;
	  }
	  double ll_s = mathops::NEGATIVE_INFINITY;
	  for(int si = 0; si < num_slots_; ++si) {
	    double x = gsl_sf_log(prior_slot_[ti][si]);
	    for(auto rhiter : floating_rel_hist_) {
	      const int rel_index = rhiter.first;
	      const int rel_count = rhiter.second;
	      x += rel_count*gsl_sf_log(prior_relation_kind_[si][rel_index]);
	      num_obs += rel_count;
	    }
	    ll_s = mathops::log_add(ll_s, x);
	  }
	  ll_t = mathops::log_add(ll_t, x0 + ll_s);
	}
	ll += ll_t;
      }
      // for every chain
      return normalized ? ll/(double)num_obs : ll;
    }


    template <typename G, typename R, typename C> 
    double latent_marginalized_loglikelihood(const InMemoryCorpus<DocumentGRC<G, R,C> >&  corpus,
					     const Vocabulary<G>& gov_vocab,
					     const Vocabulary<R>& rel_vocab,
					     bool normalized = true) {
      double ll = 0.0;
      int doc_id = 0;
      for(auto doc : corpus.get_corpus()) {
	ll += latent_marginalized_loglikelihood<G, R,C>(doc, doc_id++,
							gov_vocab, rel_vocab,
							normalized);
      }
      return ll;
    }

    double latent_and_kind_marginalized_loglikelihood(const int num_kinds, const int obs_index,
						      const std::vector<std::vector<double> >& obs_prior,
						      const std::vector<double>& kind_prior) {
      //marginalize over the latent gov_kind
      double k_ll = mathops::NEGATIVE_INFINITY;
      for(int k = 0; k < num_kinds; ++k) {
	double lp = gsl_sf_log(obs_prior[k][obs_index]) + gsl_sf_log(kind_prior[k]);
	k_ll = mathops::log_add(k_ll, lp);
      }
      return k_ll;
    }

    template <typename E> 
    double latent_and_kind_marginalized_loglikelihood(const E&  entity,
						      const Vocabulary<std::string>& gov_vocab,
						      const Vocabulary<std::string>& rel_vocab,
						      const std::vector<double>& d_t_usage) {
      double ll_t = mathops::NEGATIVE_INFINITY;
      for(int ti = 0; ti < num_templates_; ++ti) {
	double x0 = gsl_sf_log(d_t_usage[ti]);
	double ll_m_for_gov = 0.0;
	for(int mi = 0; mi < entity.num_mentions(); ++mi) {
	  const int num_gov_kinds = hyper_gov_kind_.size();
	  const int gov_index = gov_vocab.index(entity.mention(mi).gov().lemma());
	  ll_m_for_gov += latent_and_kind_marginalized_loglikelihood(num_gov_kinds, gov_index,
								     prior_governor_, prior_governor_kind_[ti]);
	}
	x0 += ll_m_for_gov;
	double ll_s = mathops::NEGATIVE_INFINITY;
	for(int si = 0; si < num_slots_; ++si) {
	  double x = gsl_sf_log(prior_slot_[ti][si]);
	    double ll_m_for_rel = 0.0;
	    for(int mi = 0; mi < entity.num_mentions(); ++mi) {
	      const int num_rel_kinds = hyper_rel_kind_.size();
	      const int rel_index = rel_vocab.index(entity.mention(mi).rel_str());
	      ll_m_for_rel += latent_and_kind_marginalized_loglikelihood(num_rel_kinds, rel_index,
									 prior_relation_, prior_relation_kind_[si]);
	    }
	    x += ll_m_for_rel;
	  ll_s = mathops::log_add(ll_s, x);
	}
	ll_t = mathops::log_add(ll_t, x0 + ll_s);
      }
      return ll_t;
    }


    template <typename G, typename R, typename C> 
    double latent_and_kind_marginalized_loglikelihood(const DocumentGRC<G,R,C>&  doc,
						      int doc_id, int* num_obs,
						      const Vocabulary<std::string>& gov_vocab,
						      const Vocabulary<std::string>& rel_vocab,
						      std::vector< std::vector<double> > corp_specific_prior_template,
						      bool normalized = true) {
      double ll = 0.0;
      int num_obs_ = 0;
      const std::vector<double>& d_t_usage = corp_specific_prior_template[doc_id];
      const int num_entities = doc.num_entities();
      for(int ei = 0; ei < num_entities; ++ei) {
	const Entity<G,R,C>& entity = doc[ei];
	double ll_t = latent_and_kind_marginalized_loglikelihood< Entity<G,R,C> >(entity, gov_vocab, rel_vocab, d_t_usage);
	num_obs_ += (entity.num_mentions()*2);
	ll += ll_t;
      }
      *num_obs += num_obs_;
      // for every chain
      return normalized ? ll/(double)num_obs_ : ll;
    }


    template <typename G, typename R, typename C> 
    double latent_and_kind_marginalized_loglikelihood(const InMemoryCorpus<DocumentGRC<G, R,C> >&  corpus,
						      const Vocabulary<std::string>& gov_vocab,
						      const Vocabulary<std::string>& rel_vocab,
						      std::vector< std::vector<double> > prior_template,
						      bool normalized = true) {
      double ll = 0.0;
      int doc_id = 0;
      int num_obs = 0;
      for(auto doc : corpus.get_corpus()) {
	ll += latent_and_kind_marginalized_loglikelihood<G, R,C>(doc, doc_id++, &num_obs,
								 gov_vocab, rel_vocab,
								 prior_template,
								 false);
      }
      return normalized ? ll/(double)num_obs : ll;
    }

    template <typename D, typename WW>
    std::vector<double> compute_coherences(const int M, const std::string& which, 
					   const InMemoryCorpus<D>& corpus,
					   const Vocabulary<WW>& vocab) {
      if(which == "gov_kind") {
	std::map<int, int> gk_occur = 
	  corpus.template gov_view_doc_occur< WW >(vocab);
	std::map<std::pair<int, int>, int> gk_cooc = 
	  corpus.template gov_view_doc_cooccur< WW >(vocab);
	return ferrum::compute_coherences(M, prior_governor_kind_, gk_occur, gk_cooc);
      }
      if(which == "gov_obs_kind_marginalized") {
	std::map<int, int> gk_occur = 
	  corpus.template gov_lemma_doc_occur(vocab);
	std::map<std::pair<int, int>, int> gk_cooc = 
	  corpus.template gov_lemma_doc_cooccur(vocab);
	//now comes the tricky part: we need to marginalize
	std::vector<std::vector< double> > marginalized;
	int n_t = prior_governor_kind_.size();
	for(int i = 0; i < n_t; ++i) {
	  std::vector<double> X;
	  const std::vector<double>& gk_i_dist = prior_governor_kind_[i];
	  int dim = gk_i_dist.size();
	  for(int obs_type = 0; obs_type < hyper_gov_.size(); ++obs_type) {
	    double p = 0.0;
	    for(int j = 0; j < dim; ++j) {
	      p += gk_i_dist[j] * prior_governor_[j][obs_type];
	    }
	    X.push_back(p);
	  }
	  marginalized.push_back(X);
	}
	return ferrum::compute_coherences(M, marginalized, gk_occur, gk_cooc);
      }
      if(which == "gov_obs") {
	std::map<int, int> gk_occur = 
	  corpus.template gov_lemma_doc_occur(vocab);
	std::map<std::pair<int, int>, int> gk_cooc = 
	  corpus.template gov_lemma_doc_cooccur(vocab);
	return ferrum::compute_coherences(M, prior_governor_, gk_occur, gk_cooc);
      }
      return std::vector<double>();
    }
  };

  class DiscreteModelGlobalSlots {
  private:
    // setting model parameters
    int num_templates_;
    // how many slots (in total)
    int num_slots_;

    //// hyperparameters
    // how to use each template
    std::vector<double> hyper_theta_;
    std::vector<double> hyper_slot_;
    // how often governors show
    std::vector<double> hyper_gov_;
    bool h_gov_set_;
    // how often relations show
    std::vector<double> hyper_rel_;
    bool h_rel_set_;

    //// priors
    // Each document has a prior probability of using
    // any particular template
    std::vector< std::vector<double> > prior_template_;
    // Each template has a prior probability of using
    // any particular slot
    std::vector< std::vector<double> > prior_slot_;

    // Each governor kind has a prior probability of generating 
    // any particular lexical governor
    std::vector< std::vector<double> > prior_governor_;
    // Each relation kind has a prior probability of generating
    // any particular lexical relation
    std::vector< std::vector<double> > prior_relation_;

  public:
    DiscreteModelGlobalSlots(int n_templates, int n_slots) : num_templates_(n_templates),
							     num_slots_(n_slots),
							     hyper_theta_(std::vector<double>(n_templates, 0.1)), 
							     hyper_slot_(std::vector<double>(n_slots, 0.1)) {
    }
    DiscreteModelGlobalSlots
    (int n_templates, int n_slots,
     SymmetricHyperparams* shp) : num_templates_(n_templates),
				  num_slots_(n_slots),
				  hyper_theta_(std::vector<double>(n_templates, shp->h_theta)), 
				  hyper_slot_(std::vector<double>(n_slots, shp->h_slot)) {
    }
    ~DiscreteModelGlobalSlots() {
    }

    int num_templates() {
      return num_templates_;
    }
    int num_slots() {
      return num_slots_;
    }

    std::vector<double> hyper_theta() {
      return hyper_theta_;
    }
    std::vector<double> hyper_usage() {
      return hyper_theta();
    }
    std::vector<double> hyper_slot() {
      return hyper_slot_;
    }
    std::vector<double> hyper_gov() {
      if(! h_gov_set_) {
	ERROR << "Hyperparameters for governors are not set";
      }
      return hyper_gov_;
    }
    std::vector<double> hyper_verb() {
      return hyper_gov();
    }
    std::vector<double> hyper_rel() {
      if(! h_rel_set_) {
	ERROR << "Hyperparameters for relations are not set";
      }
      return hyper_rel_;
    }
    inline std::vector<double> hyper_arc() {
      return hyper_rel();
    }
    void hyper_theta(const std::vector<double>& h) {
      hyper_theta_ = h;
    }
    void hyper_slot(const std::vector<double>& h) {
      hyper_slot_ = h;
    }
    void hyper_rel(const std::vector<double>& h) {
      hyper_rel_ = h;
      h_rel_set_ = true;
    }
    void hyper_gov(const std::vector<double>& h) {
      hyper_gov_ = h;
      h_gov_set_ = true;
    }
    DiscreteModelGlobalSlots& hyper_gov(double num, double val) {
      hyper_gov_ = std::vector<double>(num, val);
      h_gov_set_ = true;
      return *this;
    }
    DiscreteModelGlobalSlots& hyper_rel(double num, double val) {
      hyper_rel_ = std::vector<double>(num, val);
      h_rel_set_ = true;
      return *this;
    }

    //transfer functions 
    inline void prior_template_usage(const std::vector< std::vector< double > >& p) {
      prior_template_ = p;
    }
    inline void prior_gov(const std::vector< std::vector< double > >& p) {
      prior_governor_ = p;
    }
    inline void prior_slot_usage(const std::vector< std::vector< double > >& p) {
      prior_slot_ = p;
    }
    inline void prior_rel(const std::vector< std::vector< double > >& p) {
      prior_relation_ = p;
    }

    inline const std::vector< std::vector< double> >& prior_template_usage() {
      return prior_template_;
    }
    inline const std::vector< std::vector< double> >& prior_governor() {
      return prior_governor_;
    }
    inline const std::vector< std::vector< double> >& prior_slot_usage() {
      return prior_slot_;
    }
    inline const std::vector< std::vector< double> >& prior_relation() {
      return prior_relation_;
    }

    void print_template_usage(std::ostream & outter = std::cout) {
      ferrum::print_2d_distribution(prior_template_, outter);
    }
    void print_slot_usage(std::ostream & outter = std::cout) {
      ferrum::print_2d_distribution(prior_slot_, outter);
    }

    template <typename W>
    void print_vocab_dist(const size_t num_per, const Vocabulary<W>& vocab, 
			  const std::vector< std::vector<double> >& topics) {
      int topic_idx = 0;
      for(const auto& topic : topics) {
	std::stringstream stream;
	stream << "Topic " << topic_idx << " ::: ";
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	for(size_t item_idx = 0; item_idx < num_per && item_idx < vocab.num_words(); ++item_idx) {
	  size_t which = sorted_topic[item_idx];
	  stream << vocab.word(which) << " [" << topic[which] << "] ";
	}
	INFO << stream.str();
	++topic_idx;
      }
    }
    template <typename W>
    void print_vocab_dist(const size_t num_per, const Vocabulary<W>& vocab,
			  std::ostream& stream,
			  const std::vector< std::vector<double> >& topics) {
      for(const auto& topic : topics) {
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	for(size_t item_idx = 0; item_idx < num_per && item_idx < vocab.num_words(); ++item_idx) {
	  size_t which = sorted_topic[item_idx];
	  stream << topic[which] << " ";
	}
	stream << "\n";
      }
    }
    template <typename W> void print_verbs(const int num_per, const Vocabulary<W>& vocab) {
      print_vocab_dist(num_per, vocab, prior_governor_);
    }
    template <typename W> void print_arcs(const int num_per, const Vocabulary<W>& vocab) {
      print_vocab_dist(num_per, vocab, prior_relation_);
    }
    template <typename W> void print_verbs(const int num_per, const Vocabulary<W>& vocab, const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, dists);
    }
    template <typename W> void print_arcs(const int num_per, const Vocabulary<W>& vocab, const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, dists);
    }
    template <typename W> void print_verbs(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream) {
      print_vocab_dist(num_per, vocab, stream, prior_governor_);
    }
    template <typename W> void print_arcs(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream) {
      print_vocab_dist(num_per, vocab, stream, prior_relation_);
    }
    template <typename W> void print_verbs(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream,
					   const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, stream, dists);
    }
    template <typename W> void print_arcs(const int num_per, const Vocabulary<W>& vocab, std::ostream& stream,
					  const std::vector<std::vector<double> >& dists) {
      print_vocab_dist(num_per, vocab, stream, dists);
    }

    template <typename PrintingStruct>
    void print_templates(std::ostream& myfile,
			 const PrintingStruct* const ps) {
      myfile << "template\tslot\twhich_type\twhich_val\tprob\n";
      for(int idx_template = 0; idx_template < num_templates_; ++idx_template) {
	// print slot
	std::vector<size_t> sorted = ferrum::sort_indices(prior_slot_[idx_template], false);
	for(const auto& item_idx : sorted) {
	  myfile << idx_template << "\t-\tslot\t" << item_idx << "\t" << prior_slot_[idx_template][item_idx] << "\n";
	}
	// print verbs (govs)
	{
	  const std::vector<double>& topic = prior_governor_[idx_template];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), (unsigned long)ps->num_per_gov);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << idx_template << "\t-\tverb\t" << ps->g_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
      }
      // print arcs (relations)
      {
	for(int slot = 0; slot < num_slots_; ++slot) {
	  const std::vector<double>& topic = prior_relation_[slot];
	  std::vector<size_t> sorted = ferrum::sort_indices(topic, false);
	  const size_t max_iter = ferrum::min(sorted.size(), (unsigned long)ps->num_per_rel);
	  for(size_t item_idx = 0; item_idx < max_iter; ++item_idx) {
	    size_t which = sorted[item_idx];
	    myfile << "-\t" << slot << "\trel\t" << ps->r_vocab->word(which) << "\t" << topic[which] << "\n";
	  }
	}
      }
    }
    template <typename PrintingStruct>
    void print_templates(const std::string& f_name,
			 const PrintingStruct* const ps) {
      std::ofstream myfile;
      myfile.open(f_name);
      this->print_templates(myfile, ps);
      myfile.close();
      INFO << "Wrote templates to " << f_name;
    }
  };
}

#include "ferrum/crtlda_defs.tcc"

#endif // FERRUM_LIBNAR_CRTLDA_DEFS_H_
