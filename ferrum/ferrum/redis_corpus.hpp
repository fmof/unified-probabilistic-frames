#ifndef FERRUM_REDIS_CORPUS_H_
#define FERRUM_REDIS_CORPUS_H_

#include "ferrum/crtlda_defs.hpp"
#include "ferrum/data_util.hpp"
#include "ferrum/lock.hpp"

#include <boost/algorithm/string.hpp>
#include <string>
#include <utility>
#include <vector>

namespace ferrum {
  extern const std::string REDIS_CORPUS_DOC_LIST; // document_list
  extern const std::string REDIS_CORPUS_DOC_HASH_PREFIX; // document:
  extern const std::string REDIS_CORPUS_DOC_FIELD; // docstr

  static inline std::string create_doc_hset_id(const std::string& id) {
    return REDIS_CORPUS_DOC_HASH_PREFIX + id;
  }

  // forward declare this
  template <typename Doc>
  class RedisCorpus;

  template <typename DocType>
  class RedisCorpusIterator {
  private:
    const std::vector<std::string>* doc_names_; //shared
    std::vector<std::string>::const_iterator it_;
    std::shared_ptr<ferrum::db::Redis> redis_;
    std::string doc_str_;

    void update_document(const std::string& name) {
      ferrum::db::RedisQuery rq
	(
	 create_doc_hset_id(name),
	 doc_str_
	 );
      rq.hget();
      try {
	(*redis_)(rq); // submit the query
      } catch(const ferrum::db::DBError& e) {
	ERROR << "Exception in trying to query redis server, with a name of \"" << name << "\", using the document list \"" << create_doc_hset_id(name) << "\", and with the field \"" << doc_str_ << "\"";
	throw e;
      }
      delete document;
      document = new DocType;
      ferrum::thrift::thrift_struct_from_string<ferrum::thrift::TCompactProtocol>
	(
	 rq.value(),
	 document
	 );
    }
    RedisCorpusIterator<DocType>(const std::string& ds = REDIS_CORPUS_DOC_FIELD);
  public:
    RedisCorpusIterator<DocType>(const std::vector<std::string>* doc_names, std::shared_ptr<ferrum::db::Redis> db, const std::string& ds = REDIS_CORPUS_DOC_FIELD);
    RedisCorpusIterator<DocType>(const RedisCorpus<DocType>* corpus);

    RedisCorpusIterator<DocType>(const CRTPCorpus< DocType, RedisCorpus, ::ferrum::RedisCorpusIterator >* corpus);

    RedisCorpusIterator<DocType>(const RedisCorpusIterator<DocType>& oit);
    RedisCorpusIterator<DocType>& operator=(const RedisCorpusIterator<DocType>& oit);
    RedisCorpusIterator<DocType>(RedisCorpusIterator<DocType>&& oit);
    RedisCorpusIterator<DocType>& operator=(RedisCorpusIterator<DocType>&& oit);
    ~RedisCorpusIterator();
    static RedisCorpusIterator<DocType> end(const std::vector<std::string>* doc_names, const std::string& ds) {
      RedisCorpusIterator<DocType> ci;
      ci.doc_names_ = doc_names;
      ci.it_ = ci.doc_names_->end();
      ci.iteration_idx = ci.doc_names_->size();
      ci.doc_str_ = ds;
      return ci;
    }
    static RedisCorpusIterator<DocType> end(const RedisCorpus<DocType>* corpus, const std::string& ds = REDIS_CORPUS_DOC_FIELD) {
      return RedisCorpusIterator<DocType>::end(&(corpus->doc_names_), ds);
    }
    static RedisCorpusIterator<DocType> end(const CRTPCorpus<DocType, RedisCorpus, ::ferrum::RedisCorpusIterator >* corpus, const std::string& ds = REDIS_CORPUS_DOC_FIELD) {
      return RedisCorpusIterator<DocType>::end(static_cast< const RedisCorpus<DocType>* >(corpus));
    }
    typedef const RedisCorpusIterator<DocType>& reference;
    typedef const RedisCorpusIterator<DocType>* pointer;

    DocType* document;
    unsigned int iteration_idx;
    
    pointer operator->() const {
      return this;
    }
    reference operator*() const {
      return *this;
    }
    RedisCorpusIterator<DocType>& operator++() { //prefix
      ++it_;
      if( it_ == doc_names_->end()) {
	return *this;
      }
      update_document(*it_);
      ++iteration_idx;
      return *this;
    }
    friend bool operator==(const RedisCorpusIterator<DocType>& i1, const RedisCorpusIterator<DocType>& i2) {
      return i1.it_ == i2.it_;
    }
    friend bool operator!=(const RedisCorpusIterator<DocType>& i1, const RedisCorpusIterator<DocType>& i2) {
      return i1.it_ != i2.it_;
    }
  };

  // template <typename Doc>
  // class RedisBackgroundComputer;

  template <typename Doc>
  class RedisCorpus : public CRTPCorpus< Doc, RedisCorpus, RedisCorpusIterator > {
  private:
    Doc document_;
    Mutex unify_mutex_;
    std::set<std::string> doc_set_;
  protected:
    std::shared_ptr<ferrum::db::Redis> redis_; // Crucial attribute
    std::vector<std::string> doc_names_;
    int load_batch_size_; // crucial attribute
    std::string doc_str_; // crucial attribute
    friend class RedisCorpusIterator<Doc>;
    class RedisCorpusStats {
    public:
      int num_docs_;
      int num_atoms_;
      bool atoms_computed;
      RedisCorpusStats() : num_docs_(0), num_atoms_(0), atoms_computed(false) {
      }
    };
    RedisCorpusStats rcs_;

    // synchronized
    void update_doc_stats(const std::string& d_name) {
      std::lock_guard<ferrum::Mutex> lock(this->mut_);
      doc_names_.push_back(d_name);
      ++(rcs_.num_docs_);
      doc_set_.emplace(d_name);
    }
  public:
    RedisCorpus(const std::string& name, std::shared_ptr<ferrum::db::Redis> redis_db);
    RedisCorpus(const std::string& name, const ferrum::db::Address& addr);
    RedisCorpus(const ferrum::db::Address& addr);
    RedisCorpus();
    RedisCorpus(const CRTPCorpus< Doc, ::ferrum::RedisCorpus, ::ferrum::RedisCorpusIterator >* corpus);
    /**
     * Shuffle the documents around. This invalidates any iterators.
     */
    void shuffle();
    RedisCorpus<Doc>& doc_str(const std::string& s);
    const std::string& doc_str() const;

    typedef typename CRTPCorpus< Doc, ::ferrum::RedisCorpus, ::ferrum::RedisCorpusIterator >::const_iterator const_iterator;
    typedef Doc DocType;

    //////////////////////////////////

    /**
     * Return the number of documents loaded
     */
    size_t load_from_key(const std::string& corpus_name, bool use_prefix = true, bool print_doc_ids = true, int num = -1);
    /**
     * This resets the vocabulary.
     */
    template <typename W>
    size_t _load_vocab(const std::string& base, ferrum::Vocabulary<W>& vocab) {
      if(! vocab.allow_new_words() ) {
	ERROR << "To load a vocabulary, it must allow new words.";
	return vocab.num_words();
      }
      vocab = ferrum::Vocabulary<W>();
      ferrum::db::RedisQuery nw(base, "num_words");
      nw.hget();
      (*redis_)(nw);
      size_t expected_num_words = (size_t)std::stoi( nw.value() );
      INFO << "Loading vocab " << base << ", expecting " << expected_num_words;
      vocab.__force_set_size(expected_num_words);

      ferrum::db::RedisQuery oov_index(base, "oov_index");
      oov_index.hget();
      (*redis_)(oov_index);
      int oov_index_i = std::stoi(oov_index.value());
      if(oov_index_i >= 0) {
	ferrum::db::RedisQuery oov_val(base, "oov_value");
	oov_val.hget();
	(*redis_)(oov_val);
	vocab.__force_set_oov( W(oov_val.value().c_str(), oov_val.value().size()),
			       oov_index_i);
	INFO << "Loading vocab " << base << ", using oov symbol " << oov_val.value() << " with index " << oov_index_i;
      } else {
	INFO << "Loading vocab " << base << ", using default oov symbol " << vocab.word(vocab.oov_index()) << " with index " << vocab.oov_index();
      }

      ferrum::db::RedisQuery word_save(base, "word_list");
      word_save.hget();
      (*redis_)(word_save);
      
      std::vector<std::string> words;
      boost::split(words, word_save.value(), boost::is_any_of(" "));
      if(words.size() != expected_num_words) {
	ERROR << "We expect " << expected_num_words << " but got back " << words.size() ;
	throw 5;
      }
      for(size_t i = 0; i < expected_num_words; ++i) {
	if((int)i == oov_index_i) continue;
	vocab.__force_set_word(words[i], i);
      }

      return (size_t)vocab.num_words();
    }
    template <typename W>
    size_t load_vocab(const std::string& corpus_name, const std::string& vocab_name, ferrum::Vocabulary<W>& vocab) {
      std::string base(	 vocab_name + ":" + corpus_name + ":" + vocab_name );
      return _load_vocab(base, vocab);
    }
    template <typename W>
    void load_vocab(const std::string& vocab_name, ferrum::Vocabulary<W>& vocab) {
      load_vocab(this->name_, vocab_name, vocab);
    }
    template <typename W>
    void _save_vocab(const std::string& base, const ferrum::Vocabulary<W>& vocab) {
      size_t i = 0;
      ferrum::db::RedisQuery query
	(
	 base,
	 "num_words",
	 std::to_string(vocab.all_words().size())
	 );
      query.hset();

      ferrum::db::RedisQuery oov_index
	(
	 base,
	 "oov_index",
	 std::to_string(vocab.oov_index())
	 );
      oov_index.hset();
      query += oov_index;

      std::string oov_val_str = (vocab.oov_index() >= 0) ? std::string( vocab.word(vocab.oov_index()) ) : "";
      ferrum::db::RedisQuery oov_val
	(
	 base,
	 "oov_value",
	 oov_val_str.c_str(),
	 oov_val_str.size()
	 );
      oov_val.hset();
      (*redis_)(oov_val);

      std::stringstream ss;
      for(const W& word : vocab.all_words() ) {
	ss << word;
	++i;
	if(i < vocab.all_words().size()) {
	  ss << " ";
	}
      }
      ferrum::db::RedisQuery word_save
	(
	 base,
	 "word_list",
	 ss.str().c_str(),
	 ss.str().size()
	 );
      word_save.hset();
      query += word_save;
      (*redis_)(query);
    }
    template <typename W>
    void save_vocab(const std::string& corpus_name, const std::string& vocab_name, const ferrum::Vocabulary<W>& vocab) {
      std::string base(	 vocab_name + ":" + corpus_name + ":" + vocab_name );
      _save_vocab(base, vocab);
    }
    template <typename W>
    void save_vocab(const std::string& vocab_name, const ferrum::Vocabulary<W>& vocab) {
      save_vocab(this->name_, vocab_name, vocab);
    }

    /**
     * Unify the provided vocabularies against the stored vocabulary,
     * and change the documents according a given DocMapper instancer.
     * DocMapper must have operator()(Doc&, const std::vector< std::vector<int> >&) defined
     */
    template <typename W, typename DocMapper>
    void unify_vocabs(const std::string& unification_name, const std::vector<std::string>& vnames, std::vector<ferrum::Vocabulary<W>* >& vocabs, DocMapper dm, const std::string& final_docstr = REDIS_CORPUS_DOC_FIELD);

    // load the next batch from the DB
    virtual RedisCorpus<Doc>& operator++() {
      return *this;
    }

    virtual bool has_document(const std::string& id) {
      return doc_set_.count(id) > 0;
    }

    virtual void copy_crucial_attributes_into(RedisCorpus<Doc>* receiver) const {
      CRTPCorpus<Doc, ::ferrum::RedisCorpus, ::ferrum::RedisCorpusIterator>::copy_crucial_attributes_into(receiver);
      receiver->redis_ = redis_;
      receiver->doc_str_ = doc_str_;
      receiver->load_batch_size_ = load_batch_size_;
    }

    virtual void add_document(const Doc& document) {
      ferrum::db::RedisQuery query
	(
	 REDIS_CORPUS_DOC_LIST + ":" + this->name_,
	 document.id
	 );
      query.rpush();
      std::string safe_thrift_string( ferrum::thrift::thrift_struct_to_string<ferrum::thrift::TCompactProtocol>(document) );
      ferrum::db::RedisQuery doc_str
	(
	 create_doc_hset_id(document.id),
	 doc_str_,
	 safe_thrift_string.c_str(),
	 safe_thrift_string.size()
	 );
      query += (doc_str.hset());
      (*redis_)(query);
      update_doc_stats(document.id);
    }

    // TODO: this *should* be synchronized, but what about if
    // it's called from fill/2?
    virtual const Doc& operator[](const size_t idx) {
      const std::string& name = doc_names_.at(idx);
      ferrum::db::RedisQuery doc_str
	(
	 create_doc_hset_id(name),
	 doc_str_
	 );
      doc_str.hget();
      (*redis_)(doc_str);
      ferrum::thrift::thrift_struct_from_string<ferrum::thrift::TCompactProtocol>
	(
	 doc_str.value(),
	 &document_
	 );
      return document_;
    }

    // synchronize
    virtual void fill(const size_t idx, Doc* doc) {
      if(this->mt_) {
	std::lock_guard<ferrum::Mutex> lock(this->mut_);
	*doc = const_cast< RedisCorpus<Doc>* >(this)->operator[](idx);
      } else {
	*doc = const_cast< RedisCorpus<Doc>* >(this)->operator[](idx);
      }
    }
    virtual const int num_docs() {
      return rcs_.num_docs_;
    }
    virtual const int num_docs() const {
      return const_cast< RedisCorpus<Doc>* >(this)->num_docs();
    }
    const std::vector<std::string>& doc_names() {
      return doc_names_;
    }

    virtual const_iterator begin() const {
      const_iterator ci(&doc_names_, redis_, doc_str_);
      return ci;
    }
    virtual const_iterator end() const {
      return const_iterator::end(&doc_names_, doc_str_);
    }

    virtual RedisCorpus<Doc> subset(size_t start_i, size_t end_i) const {
      RedisCorpus<Doc> sub(this);
      sub.reset();
      for(size_t i = start_i; i < end_i && i < doc_names_.size(); ++i) {
	const std::string& dname = this->doc_names_.at(i);
	sub.update_doc_stats( dname );
      }
      return sub;
    }
    virtual RedisCorpus<Doc> subset(size_t start_i, size_t end_i) {
      return const_cast<const RedisCorpus<Doc>*>(this)->subset(start_i, end_i);
    }

    // synchronize
    virtual void reset() {
      std::lock_guard<ferrum::Mutex> lock(this->mut_);
      doc_names_.clear();
      doc_set_.clear();
      rcs_ = RedisCorpusStats();
      document_ = Doc();
    }

    std::shared_ptr<ferrum::db::Redis> get_db() const {
      return redis_;
    }

    //friend class RedisBackgroundComputer<Doc>;
  };

  template <typename DocType>
  inline int get_num_mentions(RedisCorpus<DocType>* corpus) {
    int nm = 0;
    for(auto it = corpus->begin(); it != corpus->end(); ++it) {
      const DocType& document = *(it->document);
      nm += get_num_mentions(document);
    }
    return nm;
  }

  template <typename Doc, typename Vocab>
  class RedisBackgroundComputer : public BackgroundComputer< RedisCorpus<Doc>, Vocab > {
  public:
    RedisBackgroundComputer();
    void serialize_background(RedisCorpus<Doc>* rcdb, const std::stringstream& ss, unsigned int expected_size, const std::vector<double>& counts, const std::string& bck_name);

    virtual void after(const Vocab*, RedisCorpus<Doc>*, const std::vector<double>*, const std::string&);
  };

}

#endif

#include "ferrum/redis_corpus.tcc"
