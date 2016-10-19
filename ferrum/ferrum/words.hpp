#ifndef FERRUM_WORDS_HPP_
#define FERRUM_WORDS_HPP_

#include "ferrum/minsky.hpp"
#include "ferrum/util.hpp"

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
  class StopWordList {
  private:
    typedef std::string string;
    std::unordered_set<string> words;
    void init(const string& fpath, const string& delim);
  public:
    StopWordList();
    StopWordList(const string& fpath);
    const bool contains(const string& word) const;
    const int num_words() const;
  };

  class WordMapper {
  private:
    typedef std::string string;
    std::unordered_map<string, string> words;
    string OOV;
    void init(const string& fpath, const string& delim,
	      int orig_word_col, int mapped_word_col);
  public:
    WordMapper(const string& oov);
    WordMapper(const string& oov, const string& fpath, const string& delim, int orig_col, int mapped_col);
    const string representation(const string& word);
    const string representation(const string& word) const;
    void add_word(const string& orig, const string& mapped);
    const int num_words() const;
  };

  // forward declare
  template <typename W> class Vocabulary;

  class VirtualVocabulary {
  public:
    virtual ~VirtualVocabulary() {
    }
    template <typename W> const Vocabulary<W>* downcast() const {
      return static_cast< const Vocabulary<W>* >(this);
    }
  };

  template <typename W>
  class Vocabulary : public VirtualVocabulary {
  private:
    int oov_index_;
    std::unordered_map<W, int> word_idx;
    std::vector<W> words;
    bool allow_new_words_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      ar & oov_index_;
      ar & word_idx;
      ar & words;
      ar & allow_new_words_;
    }
  public:
    typedef W Element;
    Vocabulary<W>() : oov_index_(-1), allow_new_words_(true) {
    }
    Vocabulary<W>(const W& oov) : oov_index_(0), allow_new_words_(true) {
      words.push_back(oov);
      word_idx[oov] = 0;
    };
    Vocabulary<W>(const W& oov, bool allow_new) : oov_index_(0), allow_new_words_(allow_new) {
      words.push_back(oov);
      word_idx[oov] = 0;
    };
    void clear();
    int set_oov(const W& oov) {
      if(allow_new_words_) {
	const int prev_size = words.size();
	words.push_back(oov);
	word_idx[oov] = prev_size;
	oov_index_ = prev_size;
	return 0;
      } else {
	ERROR << "Vocabulary does not allow new words -- NOT adding oov symbol " << oov;
	return 1;
      }
    }
    /**
     */
    void __force_set_oov(const W& oov, int oov_i) {
      if( words.size() <= (size_t)oov_i) {
	ERROR << "`words` is not large enough to accomodate oov index " << oov_i << "; resizing now";
	words.resize(oov_i);
      }
      words[oov_i] = oov;
      word_idx[oov] = oov_i;
      oov_index_ = oov_i;
    }
    void __force_set_word(const W& word, int word_i) {
      if( words.size() <= (size_t)word_i) {
	ERROR << "`words` is not large enough to accomodate word index " << word_i << "; resizing now";
	words.resize(word_i);
      }
      words[word_i] = word;
      word_idx[word] = word_i;
    }
    void __force_set_size(size_t expected_size) {
      words.resize(expected_size);
    }
    int oov_index() const {
      return oov_index_;
    }
    void allow_new_words(bool b) {
      allow_new_words_ = b;
    }
    bool allow_new_words() {
      return allow_new_words_;
    }
    void make_word(const W& word) {
      if(allow_new_words_ && word_idx.find(word) == word_idx.end()) {
	const int prev_size = words.size();
	words.push_back(word);
	word_idx[word] = prev_size;
      }    
    }
    int operator()(const W& word) {
      int ret = oov_index_;
      if(allow_new_words_) {
	typename std::unordered_map<W, int>::iterator it;
	if((it = word_idx.find(word)) == word_idx.end()) {
	  const int prev_size = words.size();
	  words.push_back(word);
	  word_idx[word] = prev_size;
	  ret = prev_size;
	} else {
	  ret = it->second;
	}
      }
      return ret;
    }
    /**
     * Return a mapping of old_index --> nex_index
     */
    std::vector<int> update_with(const Vocabulary<W>& other) {
      const size_t ow = (size_t)other.num_words();
      std::vector<int> mapper(ow);
      const size_t num_existing_words = num_words();
      // pivot the other vocab's indices with this vocab
      for(size_t i = 0; i < ow; ++i) {
	int nindex = -1;
	if((int)i == other.oov_index_ ) {
	  nindex = oov_index_;
	  INFO << "Vocabulary updating: " << i << " --> OOV --> " << nindex;
	} else {
	  nindex = this->operator()( other.words[i] );
	  INFO << "Vocabulary updating: " << i << " --> " << other.words[i] << " --> " << nindex;
	}
	if(nindex < 0 && oov_index_ != nindex) {
	  ERROR << "There's an indexing error going on!";
	  throw 25;
	}
	mapper[i] = nindex;
      }
      INFO << "Updated vocabulary has " << num_words() << " words (the union of vocabs with " << num_existing_words << " and " << ow << " words)";
      return mapper;
    }
    inline const W& word(const size_t& idx) const {
      return words.at(idx);
    }
    inline const W& word(const size_t& idx) {
      return words[idx];
    }
    inline const int index(const W& word) {
      return (word_idx.find(word) == word_idx.end()) ? oov_index_ : word_idx[word];
    }
    inline const int index(const W& word) const {
      auto finder = word_idx.find(word);
      return (finder == word_idx.end()) ? oov_index_ : finder->second;
    }

    inline const int num_words() const {
      return words.size();
    }
    const std::vector<W>& all_words() const {
      return words;
    }
    minsky::Vocab minskify() const;
    size_t reset_with(const minsky::Vocab&);
    template <typename V> friend bool operator==(const Vocabulary<V>& v1, const Vocabulary<V>& v2);

    unsigned int serialize_string(const std::vector<double>& counts, std::stringstream& ss) const {
      unsigned int expected_size = 0;
      for(size_t i = 0; i < counts.size(); ++i) {
	const W& word = this->word(i);
	std::string countstring(std::to_string( counts[i] ));
	ss << word << " " << countstring << " ";
	expected_size += word.size() + countstring.size() + 2;
      }
      return expected_size;
    }
  };

  template <typename W>
  bool operator==(const Vocabulary<W>& v1, const Vocabulary<W>& v2);

  template <typename T>
  class AnnotatedToken {
  private:
    std::string lemma_;
    std::string original_;
    std::string pos_;
    T view_;
    friend std::ostream & operator<<(std::ostream &os, const AnnotatedToken<T>& token) {
      return os << "{view:" << token.view_ << ", lemma:" << token.lemma_ << "}";
    }
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
}

#endif

#include "ferrum/words.tcc"
