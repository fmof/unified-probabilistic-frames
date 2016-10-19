#ifndef FERRUM_WORDS_TCC_
#define FERRUM_WORDS_TCC_

namespace ferrum {
  // I'd like to use the following specialization, but I ran into ODR problems
  // template <>
  // inline minsky::Vocab Vocabulary<std::string>::minskify() const {
  //   minsky::Vocab voc;
  //   voc.__set_words(words);
  //   if(oov_index_ >= 0) {
  //     voc.__set_oov_index( oov_index_ );
  //     voc.__set_oov( words[oov_index_] );
  //   }
  //   return voc;
  // }

  template <typename W>
  bool operator==(const Vocabulary<W>& v1, const Vocabulary<W>& v2) {
    if(&v1 == &v2) return true;
    if(v1.oov_index_ != v2.oov_index_) return false;
    if(! ferrum::vectors_equal(v1.words, v2.words) ) return false;
    if(! ferrum::maps_equal(v1.word_idx, v2.word_idx) ) return false;
    if(v1.allow_new_words_ != v2.allow_new_words_) return false;
    return true;
  }

  template <typename W>
  minsky::Vocab Vocabulary<W>::minskify() const {
    minsky::Vocab voc;
    std::vector<std::string> w;
    for(const W& word : words) {
      w.push_back( word );
    }
    voc.__set_words(w);
    if(oov_index_ >= 0) {
      voc.__set_oov_index( oov_index_ );
      voc.__set_oov( words[oov_index_] );
    }
    return voc;
  }

  template <typename W>
  void Vocabulary<W>::clear() {
    oov_index_ = -1;
    allow_new_words_ = true;
    words.clear();
    word_idx.clear();
  }

  template <typename W>
  size_t Vocabulary<W>::reset_with(const minsky::Vocab& mvoc) {
    if(! allow_new_words() ) {
      ERROR << "To reset a vocabulary, it must allow new words.";
      return num_words();
    }
    clear();
    if(! mvoc.__isset.words) {
      ERROR << "minsky::Vocab.words must be set";
      throw 10;
    }
    size_t expected_num_words = mvoc.words.size();
    __force_set_size(expected_num_words);

    int oov_index_i = mvoc.__isset.oov_index ? mvoc.oov_index : -1;
    if(oov_index_i >= 0) {
      if(! mvoc.__isset.oov) {
	ERROR << "oov must be set when oov_index is >= 0";
	throw 11;
      }
      __force_set_oov( W(mvoc.oov.c_str(), mvoc.oov.size()),
			     oov_index_i);
    }

    for(size_t i = 0; i < expected_num_words; ++i) {
      if((int)i == oov_index_i) continue;
      __force_set_word(mvoc.words[i], i);
    }

    return (size_t)num_words();
  }
}

#endif
