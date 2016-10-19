#ifndef FERRUM_REDIS_CORPUS_TCC_
#define FERRUM_REDIS_CORPUS_TCC_

namespace ferrum {
  template <typename Doc, typename V>
  RedisBackgroundComputer<Doc, V>::RedisBackgroundComputer() {
  }

  template <typename Doc, typename V>
  void RedisBackgroundComputer<Doc, V>::serialize_background
  (
   RedisCorpus<Doc>* rcdb,
   const std::stringstream& ss,
   unsigned int expected_size,
   const std::vector<double>& counts,
   const std::string& bck_name
   ) {
    ferrum::db::RedisQuery bckgrnd
      (
       rcdb->get_name(),
       bck_name,
       ss.str().c_str(),
       (size_t)expected_size
       );
    bckgrnd.hmset();
    rcdb->get_db()->operator()(bckgrnd);
  }

  template <typename Doc, typename V>
  void RedisBackgroundComputer<Doc, V>::after
  (
   const V* vocab,
   RedisCorpus<Doc>* rcdb,
   const std::vector<double>* counts,
   const std::string& counts_name
) {
    std::stringstream ss;      
    unsigned int expected_size = vocab->serialize_string(*counts, ss);
    serialize_background(rcdb, ss, expected_size, *counts, counts_name);
  }

  /////////////////////////////////////////////////////////////////////////

  template <typename DocType>
  RedisCorpusIterator<DocType>::RedisCorpusIterator(const std::string& ds) :
    doc_str_(ds), document(NULL) {
  }

  template <typename DocType>
  RedisCorpusIterator<DocType>::RedisCorpusIterator(const RedisCorpus<DocType>* corpus) :
    RedisCorpusIterator<DocType>(&(corpus->doc_names_), corpus->redis_, corpus->doc_str()) {
  }

  template <typename DocType>
  RedisCorpusIterator<DocType>::RedisCorpusIterator(const std::vector<std::string>* doc_names, std::shared_ptr<ferrum::db::Redis> db, const std::string& ds) :
    doc_names_(doc_names),
    it_(doc_names->begin()),
    redis_(db),
    doc_str_(ds),
    document(new DocType),
    iteration_idx(0) {
    if(it_ != doc_names->end()) {
      update_document(*it_);
    } else {
      *this = RedisCorpusIterator<DocType>::end(doc_names, ds);
    }
  }

  template <typename DocType>
  RedisCorpusIterator<DocType>::RedisCorpusIterator(const CRTPCorpus< DocType, RedisCorpus, ::ferrum::RedisCorpusIterator >* corpus) :
    RedisCorpusIterator<DocType>(static_cast< const RedisCorpus<DocType>* >(corpus) ) {
  }

  template <typename DocType>
  RedisCorpusIterator<DocType>::RedisCorpusIterator(const RedisCorpusIterator<DocType>& oit) :
    doc_names_( oit.doc_names_),
    it_(oit.it_),
    redis_(oit.redis_),
    doc_str_(oit.doc_str_),
    document(NULL),
    iteration_idx(oit.iteration_idx) {
    if(oit.document != NULL) {
      document = new DocType(*(oit.document));
    }
  }

  template <typename DocType>
  RedisCorpusIterator<DocType>& RedisCorpusIterator<DocType>::operator=(const RedisCorpusIterator<DocType>& oit) {
    if(this == &oit) {
      return *this;
    }
    doc_names_ = oit.doc_names_;
    it_ = oit.it_; 
    redis_ = oit.redis_;
    doc_str_ = oit.doc_str_;
    iteration_idx = oit.iteration_idx;
    if(oit.document != NULL) {
      document = new DocType(*(oit.document));
    } else {
      document = NULL;
    }
    return *this;
  }

  template <typename DocType>
  RedisCorpusIterator<DocType>::RedisCorpusIterator(RedisCorpusIterator<DocType>&& oit) :
    doc_names_( oit.doc_names_),
    it_( std::move(oit.it_) ),
    redis_( std::move(oit.redis_) ),
    doc_str_( std::move(oit.doc_str_) ),
    document(NULL),
    iteration_idx(oit.iteration_idx) {
    if(oit.document != NULL) {
      document = new DocType(std::move( *(oit.document) ) );
      oit.document = NULL;
    }
  }

  template <typename DocType>
  RedisCorpusIterator<DocType>& RedisCorpusIterator<DocType>::operator=(RedisCorpusIterator<DocType>&& oit) {
    if(this != &oit) {
      doc_names_ = oit.doc_names_;
      it_ = std::move(oit.it_); 
      redis_ = std::move(oit.redis_);
      doc_str_ = std::move( oit.doc_str_ ),
	iteration_idx = oit.iteration_idx;
      if(oit.document != NULL) {
	document = new DocType( std::move( *(oit.document) ) );
	delete oit.document;
      } else {
	document = NULL;
      }
    }
    return *this;
  }

  template <typename DocType>
  RedisCorpusIterator<DocType>::~RedisCorpusIterator() {
    delete document;
  }

  /////////////////////////////////////////////////////////////////////////////

  template <typename Doc>
  RedisCorpus<Doc>::RedisCorpus(const std::string& name, std::shared_ptr<ferrum::db::Redis> redis_db) :
    CRTPCorpus< Doc, ::ferrum::RedisCorpus, ::ferrum::RedisCorpusIterator >(name),
    redis_(redis_db), doc_str_(REDIS_CORPUS_DOC_FIELD) {
  }
  template <typename Doc>
  RedisCorpus<Doc>::RedisCorpus(const std::string& name, const ferrum::db::Address& addr) :
    CRTPCorpus< Doc, ::ferrum::RedisCorpus, ::ferrum::RedisCorpusIterator >(name),
    redis_(new ferrum::db::Redis(addr)),
    doc_str_(REDIS_CORPUS_DOC_FIELD)  {
  }

  template <typename Doc>
  RedisCorpus<Doc>::RedisCorpus(const ferrum::db::Address& addr) :
    RedisCorpus<Doc>("", addr) {
  }

  template <typename Doc>
  RedisCorpus<Doc>::RedisCorpus() :
    CRTPCorpus< Doc, ::ferrum::RedisCorpus, ::ferrum::RedisCorpusIterator >(""),
    redis_(), doc_str_(REDIS_CORPUS_DOC_FIELD) {
  }

  template <typename Doc>
  RedisCorpus<Doc>::RedisCorpus(const CRTPCorpus< Doc, ::ferrum::RedisCorpus, ::ferrum::RedisCorpusIterator >* corpus) :
    RedisCorpus<Doc>( *(corpus->downcast()) ) {
  }

  template <typename Doc>
  void RedisCorpus<Doc>::shuffle() {
    std::lock_guard<ferrum::Mutex> lock(this->mut_);
    std::random_shuffle(doc_names_.begin(), doc_names_.end());
  }

  template <typename Doc>
  RedisCorpus<Doc>& RedisCorpus<Doc>::doc_str(const std::string& s) {
    doc_str_ = s;
    return *this;
  }

  template <typename Doc>
  const std::string& RedisCorpus<Doc>::doc_str() const {
    return doc_str_;
  }

  /**
   * Return the number of documents loaded
   */
  template <typename Doc>
  size_t RedisCorpus<Doc>::load_from_key(const std::string& corpus_name, bool use_prefix, bool print_doc_ids, int num) {
    ferrum::db::RedisQuery query
      (
       use_prefix ? (REDIS_CORPUS_DOC_LIST + ":" + corpus_name) : corpus_name,
       "0 -1"
       );
    query.lrange();
    (*redis_)(query);
    int i = 0;
    for(const std::string& d_name : query.vec_value() ) {
      if(i == num) break;
      update_doc_stats(d_name);
      if(print_doc_ids) {
	INFO << d_name;
      }
      ++i;
    }
    return rcs_.num_docs_;
  }

  /**
   * Unify the provided vocabularies against the stored vocabulary,
   * and change the documents according a given DocMapper instancer.
   * DocMapper must have operator()(Doc&, const std::vector< std::vector<int> >&) defined
   */
  template <typename Doc>
  template <typename W, typename DocMapper>
  void RedisCorpus<Doc>::unify_vocabs
  (
   const std::string& unification_name,
   const std::vector<std::string>& vnames,
   std::vector<ferrum::Vocabulary<W>* >& vocabs,
   DocMapper dm,
   const std::string& final_docstr
   ) {
    if(vnames.size() != vocabs.size()) {
      ERROR << "Cannot unify " << vocabs.size() << " vocabularies with " << vnames.size() << " names.";
      return;
    }
    const size_t size = vnames.size();
    std::vector<ferrum::Vocabulary<W>*> nvocabs;
    std::vector<std::vector<int> > word_maps;
    enum Branches {
      UNSET = -1,
      NON_EXIST,
      ALREADY_EXIST
    };
    Branches branch = Branches::UNSET;
    bool error = false;

    ferrum::db::MultiClientRedisLock mcrl(this->name_, redis_);
    for(size_t i = 0; i < size; ++i) {
      const std::string& name = vnames[i];
      ferrum::Vocabulary<W>* vocab = vocabs[i];
      std::string base(name + ":" + unification_name);
      if(! redis_->has_key(base) ) {
	// doesn't exist yet, so save the passed-in-vocab
	_save_vocab(base, *(const_cast<ferrum::Vocabulary<W>*>(vocab)));
	nvocabs.push_back(vocab);
	if(branch == Branches::UNSET) {
	  branch = Branches::NON_EXIST;
	} else if(branch != Branches::NON_EXIST) {
	  ERROR << "Cannot have different branches";
	  error = true;
	  break;
	}
      } else {
	if(branch == Branches::UNSET) {
	  branch = Branches::ALREADY_EXIST;
	} else if(branch != Branches::ALREADY_EXIST) {
	  ERROR << "Cannot have different branches";
	  error = true;
	  break;
	}
	// make a new vocab
	ferrum::Vocabulary<W>* nvocab = new ferrum::Vocabulary<W>();
	size_t num_loaded = _load_vocab(base, *nvocab);
	INFO << "During vocab unification, loaded already set vocab with " << num_loaded << " words";
	std::vector<int> wmap = nvocab->update_with(*vocab);
	word_maps.push_back(wmap);
	nvocabs.push_back(nvocab);
	_save_vocab(base, *nvocab);
      }
    }
    if(error) {
      ERROR << "Encountered an error; returning now";
      return;
    }
    // only update if we had to unify
    if(branch == Branches::ALREADY_EXIST) {
      for(const std::string& dname : doc_names_) {
	INFO << "Grabbing " << dname << " under " << create_doc_hset_id(dname) << " and field " << doc_str_;
	ferrum::db::RedisQuery doc_str
	  (
	   create_doc_hset_id(dname),
	   doc_str_
	   );
	doc_str.hget();
	(*redis_)(doc_str);
	Doc doc_to_change;
	ferrum::thrift::thrift_struct_from_string<ferrum::thrift::TCompactProtocol>
	  (
	   doc_str.value(),
	   &doc_to_change
	   );
	dm(doc_to_change, word_maps);

	std::string safe_thrift_string( ferrum::thrift::thrift_struct_to_string<ferrum::thrift::TCompactProtocol>(doc_to_change) );
	ferrum::db::RedisQuery changed_doc_str
	  (
	   create_doc_hset_id(doc_to_change.id),
	   final_docstr,
	   safe_thrift_string.c_str(),
	   safe_thrift_string.size()
	   );
	(*redis_)(changed_doc_str.hset());
	  

	if(final_docstr != doc_str_) {
	  ferrum::db::RedisQuery hdel_q
	    (
	     create_doc_hset_id(dname),
	     doc_str_
	     );
	  hdel_q.hdel();
	  (*redis_)(hdel_q);
	}
      }
      // swap vocabs, if applicable, and free
      for(size_t i = 0; i < size; ++i) {
	if(vocabs[i] == nvocabs[i]) {
	  continue;
	}
	*(vocabs[i]) = *(nvocabs[i]);
	delete nvocabs[i];
      }
    } else if(branch == Branches::NON_EXIST) {
      if(final_docstr != doc_str_) {
	for(const std::string& dname : doc_names_) {
	  INFO << "Remapping " << dname << " under " << create_doc_hset_id(dname) << " from field " << doc_str_ << " to " << final_docstr;
	  ferrum::db::RedisQuery doc_str
	    (
	     create_doc_hset_id(dname),
	     doc_str_
	     );
	  doc_str.hget();
	  (*redis_)(doc_str);
	  ferrum::db::RedisQuery changed_doc_str
	    (
	     create_doc_hset_id(dname),
	     final_docstr,
	     doc_str.value().c_str(),
	     doc_str.value().size()
	     );
	  (*redis_)(changed_doc_str.hset());

	  ferrum::db::RedisQuery hdel_q
	    (
	     create_doc_hset_id(dname),
	     doc_str_
	     );
	  hdel_q.hdel();
	  (*redis_)(hdel_q);
	}
      }
    }
  }
}

#endif
