#include "ferrum/logging.hpp"
#include "ferrum/words.hpp"

namespace ferrum {
  void StopWordList::init(const std::string& fpath, const std::string& delim) {
    std::ifstream ifile(fpath);
    std::string line;
    std::vector<std::string> strs;
    if(ifile.is_open()) {
      while( getline(ifile,line) ) {
	boost::split(strs, line , boost::is_any_of(delim));
	std::string w = strs[0];
	words.insert(w);
	strs.clear();
      }
      ifile.close();
    } else {
      ERROR << "StopWordList was unable to open file " << fpath; 
      throw 19;
    }
  }
  StopWordList::StopWordList() {
  }
  StopWordList::StopWordList(const std::string& fpath) {
    init(fpath, "\t");
  };
  const bool StopWordList::contains(const std::string& word) const {
    return words.find(word) != words.end();
  }
  const int StopWordList::num_words() const {
    return words.size();
  }

  void WordMapper::init(const std::string& fpath, const std::string& delim,
			int orig_word_col, int mapped_word_col) {
    std::ifstream ifile(fpath);
    std::string line;
    std::vector<std::string> strs;
    if(ifile.is_open()) {
      while( getline(ifile,line) ) {
	boost::split(strs, line , boost::is_any_of(delim));
	std::string w = strs[orig_word_col];
	std::string m = strs[mapped_word_col];
	if(words.find(w) == words.end()) {
	  words[w] = m;
	} else {
	  ERROR << "Trying to overwrite existing entry: (" << w << " -> " << words[w] << ") with (" << w << " -> " << m << ")";
	}
	strs.clear();
      }
      ifile.close();
    } else {
      ERROR << "WordMapper was unable to open file " << fpath; 
      throw 19;
    }
  }

  WordMapper::WordMapper(const std::string& oov) : OOV(oov) {
  }
  WordMapper::WordMapper
  (
   const std::string& oov,
   const std::string& fpath,
   const std::string& delim,
   int orig_col,
   int mapped_col
   ) : OOV(oov) {
    init(fpath, "\t", orig_col, mapped_col);
  };
  const std::string WordMapper::representation(const std::string& word) {
      return (words.find(word) == words.end()) ? OOV : words[word];
    }
  const std::string WordMapper::representation(const std::string& word) const {
      auto finder = words.find(word);
      return (finder == words.end()) ? OOV : finder->second;
    }
  void WordMapper::add_word(const std::string& orig, const std::string& mapped) {
    if(words.find(orig) == words.end()) {
      words[orig] = mapped;
    } else {
      ERROR << "Trying to overwrite existing entry: (" << orig << " -> " << words[orig] << ") with (" << orig << " -> " << mapped << ")";
    }
      
  }
  const int WordMapper::num_words() const {
    return words.size();
  }
}
