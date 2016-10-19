// Header for word-based topic models

#ifndef FERRUM_LIBNAR_WTM_H_
#define FERRUM_LIBNAR_WTM_H_

#include "ferrum/concrete.hpp"
#include "ferrum/dmc.hpp"
#include "ferrum/mathops.hpp"
#include "ferrum/util.hpp"

#include <iostream>
#include "stdlib.h"

// for pair
#include <utility>
#include <string>
#include <vector>

namespace wtm {
  template <typename W>
  class Vocabulary {
  private:
    std::unordered_map<W, int> word_idx;
    std::vector<W> words;
  public:
    Vocabulary<W>(const W& oov) {
      words.push_back(oov);
      word_idx["__OOV__"] = 0;
    };
    void make_word(const W& word) {
      if(word_idx.find(word) == word_idx.end()) {
	const int prev_size = words.size();
	words.push_back(word);
	word_idx[word] = prev_size;
      }    
    }
    inline const W& word(const size_t& idx) const {
      return words.at(idx);
    }
    inline const W& word(const size_t& idx) {
      return words[idx];
    }
    inline const int index(const W& word) {
      return (word_idx.find(word) == word_idx.end()) ? 0 : word_idx[word];
    }
    inline const int index(const W& word) const {
      auto finder = word_idx.find(word);
      return (finder == word_idx.end()) ? 0 : finder->second;
    }

    inline const int num_words() {
      return words.size();
    }
  };

  template <typename T>
  class AnnotatedToken {
  private:
    std::string lemma_;
    std::string original_;
    std::string pos_;
    T view_;
  public:
    AnnotatedToken<T>(){
    }
    AnnotatedToken<T>(const std::string& lem, const std::string& orig, const std::string& pos) : lemma_(lem), original_(orig), pos_(pos) {
    }
    void transfer(AnnotatedToken<T>& other) {
      lemma_ = other.lemma_;
      original_ = other.original_;
      pos_ = other.pos_;
      view_ = other.view_;
    }
    void lemma(const std::string& l) {
      lemma_ = l;
    }
    void original(const std::string& orig) {
      original_ = orig;
    }
    void pos(const std::string& p) {
      pos_ = p;
    }
    void view(const T& view) {
      view_ = view;
    }
    const std::string& lemma() const {
      return lemma_;
    }
    const std::string& original() const {
      return original_;
    }
    const std::string& pos() const {
      return pos_;
    }
    const T& view() const {
      return view_;
    }
  };

  template <typename W>
  class WordPruner {
  protected:
    Vocabulary<W>* const vocab_ptr_;
    //protected:
    //   const concrete::Communication& communication;
  public:
    // WordPruner< W >(const concrete::Communication& comm) : communication(comm) {
    // }
    WordPruner< W >(const concrete::Communication& comm) {
    }
    WordPruner< W >(Vocabulary<W>* const vp) : vocab_ptr_(vp) {
    }
    virtual const W make_word_view(const std::string& word) const {
      return word;
    }
    virtual std::vector<W> prune(const concrete::Tokenization& tokenization) const {
      std::vector<W> words;
      for(concrete::Token token : tokenization.tokenList.tokenList) {
	W word = this->make_word_view(token.text);
	this->vocab_ptr_->make_word(word);
	words.push_back( word );
      }
      return words;
    }
  };

  template <typename W>
  class VerbPruner : public WordPruner<W> {
    //protected:
    //   const concrete::Communication& communication;
  public:
    // WordPruner< W >(const concrete::Communication& comm) : communication(comm) {
    // }
    VerbPruner< W >(const concrete::Communication& comm) : WordPruner<W>(comm) {
    }
    VerbPruner< W >(Vocabulary<W>* const vp) : WordPruner<W>(vp) {
    }
    virtual const W make_word_view(const std::string& word) const {
      return word;
    }
    virtual std::vector<W> prune(const concrete::Tokenization& tokenization) const {
      std::vector<W> words;
      const concrete::TokenTagging* pos = concrete_util::first_pos_tagging(tokenization, "Stanford");
      if(pos == NULL) {
	return words;
      }
      std::vector<concrete::TaggedToken> pos_tags = pos->taggedTokenList;
      for(concrete::Token token : tokenization.tokenList.tokenList) {
	if(pos_tags[token.tokenIndex].tag[0] != 'V') continue;
	W word = this->make_word_view(token.text);
	this->vocab_ptr_->make_word(word);
	words.push_back( word );
      }
      return words;
    }
  };

  template <typename W>
  class Document {
  private:
    std::vector< W > words_;
  public:
    const std::string id;
    Document< W >(const std::string id_) : id(id_) {
    }

    Document< W >(const concrete::Communication& communication,
		  const WordPruner< W >& word_pruner) : id(communication.id) {
      for(concrete::Section section : communication.sectionList) {
	if(! section.__isset.sentenceList) continue;
	for(concrete::Sentence sentence : section.sentenceList) {
	  if(! sentence.__isset.tokenization) continue;
	  concrete::Tokenization tokenization = sentence.tokenization;
	  std::vector< W > toks_to_add = word_pruner.prune(tokenization);
	  if(toks_to_add.size() > 0) {
	    words_.insert(words_.end(), toks_to_add.begin(), toks_to_add.end());
	  }
	}
      }
      INFO << "Document " << this->id << " has " << num_words() << " words";
    }

    inline const int num_words() const {
      return words_.size();
    }

    inline void add_word(const W& word) {
      words_.push_back(word);
    }

    inline const W& operator[](const size_t idx) const {
      return words_[idx];
    }
  };

  template <typename D> 
  class Corpus {
  private:
    std::vector<D> documents;
    std::string name;
    
  public:
    Corpus<D>(std::string corpus_name) : name(corpus_name) {
    }

    void add_document(D& document) {
      documents.push_back(document);
    }

    inline const std::vector<D>& get_corpus() const {
      return documents;
    }
    inline const D& operator[](const size_t idx) const {
      return documents[idx];
    }
    inline const int num_docs() {
      return documents.size();
    }
    /**
     * Return an array indicating the number of 
     * entities in each document.
     * This method allocates memory, but it is the 
     * callers responsibility to free it.
     */
    int* word_count() {
      const int nd = num_docs();
      int* e = (int*)ferrum::MALLOC(sizeof(int)*nd);
      for(int i = 0; i<nd; i++) {
	e[i] = documents[i]->num_words();
      }
      return e;
    }
  };

  class SymmetricHyperparams {
  public:
    double h_theta;
    double h_word;
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
  template <typename W>
  class DiscreteLDA {
  private:
    // setting model parameters
    const int num_topics_;

    //// hyperparameters
    // how to use each template
    double* hyper_theta_;
    // how often governors show
    double* hyper_word_;

    //// priors
    // Each document has a prior probability of using
    // any particular template
    std::vector< std::vector<double> > prior_topic_;
    // Each template has a prior probability of generating 
    // any particular governor
    std::vector< std::vector<double> > prior_word_;

    void fill_hyper(double* &array, double value, const int size) {
      ferrum::allocate_1d(array, size);
      ferrum::fill_array1d<double>(array, value, size);
    }

  public:
    DiscreteLDA<W>(int n_topic, SymmetricHyperparams* shp,
		   Vocabulary<W> *vocabs) : num_topics_(n_topic) {
      fill_hyper(hyper_theta_, shp->h_theta, num_topics_);
      fill_hyper(hyper_word_, shp->h_word, vocabs->num_words());
    }
    ~DiscreteLDA() {
      // free a bunch of stuff
    }
    inline const int num_topics() {
      return num_topics_;
    }
    inline const double* hyper_word() {
      return hyper_word_;
    }
    inline const double* hyper_theta() {
      return hyper_theta_;
    }
    
    //transfer functions
    inline void prior_topic(const std::vector< std::vector< double > >& posterior_topic) {
      prior_topic_ = posterior_topic;
    }
    inline void prior_word(const std::vector< std::vector< double > >& posterior_word) {
      prior_word_ = posterior_word;
    }

    inline const std::vector< std::vector< double> >& prior_topic() {
      return prior_topic_;
    }
    inline const std::vector< std::vector< double> >& prior_word() {
      return prior_word_;
    }

    void print_topics(const int num_per, const Vocabulary<W>& vocab) {
      int topic_idx = 0;
      for(const std::vector<double> topic : prior_word_) {
	INFO << "Topic " << topic_idx;
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	for(size_t item_idx = 0; item_idx < num_per; ++item_idx) {
	  size_t which = sorted_topic[item_idx];
	  INFO << "\t" << topic[which] << "\t" << vocab.word(which);
	}
	++topic_idx;
      }
    }
  };

  class SamplingStrategy {
  private:
    int reestimate_usage_every_;
    int reestimate_topics_every_;
  public:
    const int num_iterations;
    const int burn_in;
    virtual bool sample_topic(int, int, int) = 0;
    virtual bool reestimate_topics(int) = 0;
    virtual bool reestimate_usage(int) = 0;
    SamplingStrategy(int num_iter, int burnin) : num_iterations((const int)num_iter), 
						 burn_in((const int)burnin) {
    }
    inline const int reestimate_usage_every() {
      return reestimate_usage_every_;
    }
    inline const int reestimate_topics_every() {
      return reestimate_topics_every_;
    }
    inline void reestimate_topics_every(int i) {
      reestimate_topics_every_ = i;
    }
    inline void reestimate_usage_every(int i) {
      reestimate_usage_every_ = i;
    }
  };
  class SampleEveryIter : public SamplingStrategy {
  public:
    SampleEveryIter(int num_iter, int burnin) : SamplingStrategy(num_iter, burnin){
      this->reestimate_topics_every(100);
      this->reestimate_usage_every(100);
    }
    bool sample_topic(int i, int d, int e) {
      return true;
    }
    bool reestimate_usage(int iter_index) {
      return iter_index >= burn_in && iter_index % this->reestimate_usage_every() == 0;
    }
    bool reestimate_topics(int iter_index) {
      return iter_index >= burn_in && iter_index % this->reestimate_topics_every() == 0;
    }
  };

  /**
   * A collapsed Gibbs sampler for a discrete model.
   * This assumes that each observation is generated by a
   * categorical (discrete) distribution that is endowed
   * with a Dirichlet prior.
   *
   * The M type should be a DiscreteModel.
   */
  template <typename M, typename D, typename V>
  class CollapsedGibbsDMC {
  private:
    const int num_docs;
    SamplingStrategy* sample_strategy;

    std::vector< std::vector< int> >  assignments;
    // count tables
    std::vector< std::vector< int> >  c_doc_topic;
    // num_docs * num_topic;
    std::vector< std::vector< int> >  c_topic;
    // num_templates
    std::vector<int> c_topic_sums;

    std::vector<int>* c_doc_topic_ptr;

    // the model itself
    M* model;
    Corpus<D>* corpus_;

    int num_topics_;

    std::vector<int> num_words;

    dmc::gdmc topic_dmc;
    dmc::gdmc word_dmc;
    
    int current_word_index = -1;

  public:
    CollapsedGibbsDMC<M, D, V>(M* model, Corpus<D>* corpus, V* vocab) :
    num_docs(corpus->num_docs()), model(model), corpus_(corpus), num_topics_(model->num_topics()) {
      int num_docs = corpus->num_docs();
      for(int d = 0; d < num_docs; d++) {
	const int nw = ((*corpus)[d]).num_words();
	num_words.push_back(nw);
	c_doc_topic.push_back(std::vector<int>(num_topics_));
	assignments.push_back(std::vector<int>(nw));
      }
      const int vocab_size = vocab->num_words();
      for(int t = 0; t < model->num_topics(); t++) {
	c_topic_sums.push_back(0);
	c_topic.push_back(std::vector<int>(vocab_size));
      }
      word_dmc = dmc::gdmc(vocab_size, model->hyper_word(), model->num_topics());
      topic_dmc = dmc::gdmc(model->num_topics(), model->hyper_theta(), num_docs);
    }
    ~CollapsedGibbsDMC() {
    }
    int sample_topic() {
      std::vector<double> lps(num_topics_);
      for(int topic_idx = 0; topic_idx < num_topics_; topic_idx++) {
	double inner = word_dmc.log_u_conditional(current_word_index, c_topic[topic_idx], c_topic_sums[topic_idx], 1);
	inner += topic_dmc.log_u_conditional(topic_idx, *c_doc_topic_ptr, num_docs, 1);
	lps[topic_idx] = inner;
      }
      return dmc::cat::log_u_sample(lps);
    }

    void init(const Corpus<D>& corpus,
	      const V& vocab) {
      for(int di = 0; di < num_docs; di++){
	const int num_w = num_words[di];
	D doc = corpus[di];
	for(int wi = 0; wi < num_w; wi++){
	  // randomly assign to a topic
	  int topic = wi % model->num_topics();
	  assignments[di][wi] = topic;
	  ++c_doc_topic[di][topic];
	  ++c_topic[topic][ vocab.index(doc[wi]) ];
	  ++c_topic_sums[topic];
	}
      }
    }

    void learn(const V& vocab) {
      for(int iteration = 0; iteration < sample_strategy->num_iterations; iteration++) {
	for(int di = 0; di < num_docs; di++){
	  D doc = (*corpus_)[di];
	  c_doc_topic_ptr = &(c_doc_topic[di]);
	  const int num_w = num_words[di];
	  for(int wi = 0; wi < num_w; wi++){
	    if(sample_strategy->sample_topic(iteration, di, wi)) {
	      //get word
	      const int word = vocab.index(doc[wi]);
	      current_word_index = word;
	      //remove counts
	      const int prev = assignments[di][wi];
	      --c_doc_topic[di][prev];
	      --c_topic[prev][word];
	      --c_topic_sums[prev];
	      int sampled = sample_topic();
	      TRACE << "sampled topic:" << sampled;
	      //add counts back
	      ++c_topic[sampled][word];
	      ++c_topic_sums[sampled];
	      ++c_doc_topic[di][sampled];
	      if(assignments[di][wi] != sampled) {
		assignments[di][wi] = sampled;
	      }
	    }
	  }
	}
	if(sample_strategy->reestimate_topics(iteration)) {
	  word_dmc.reestimate_collapsed_parameters(c_topic);
	}
	if(sample_strategy->reestimate_usage(iteration)) {
	  topic_dmc.reestimate_collapsed_parameters(c_doc_topic);
	}
      }
    }

    // transfer the learned parameters back to the model
    void transfer_learned_parameters() {
      model->prior_topic(topic_dmc.collapsed_params());
      model->prior_word(word_dmc.collapsed_params());      
    }

    //setters
    void sampling_strategy(SamplingStrategy* ss) {
      sample_strategy = ss;
    }
  };

}

#endif // FERRUM_LIBNAR_CRTLDA_H_
