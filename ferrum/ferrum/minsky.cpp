#include "ferrum/logging.hpp"
#include "ferrum/minsky.hpp"

#include <map>

namespace minsky {
  bool operator<(const MinskyAnnotationWrapper& x, const MinskyAnnotationWrapper& y) {
    return std::tie(x.annot_level, x.st_level, x.word_level) < std::tie(y.annot_level, y.st_level, y.word_level);
  }
  minsky::Frame frame_unnormalized(const std::vector<double>& v) {
    minsky::Frame f;
    minsky::Distribution d;
    minsky::Weights w;
    w.__set_unnormalized(v);
    d.__set_weights(w);
    f.__set_distr(d);
    return f;
  }
  size_t num_entities(const EDoc& doc) {
    return doc.entities.size();
  }
  std::map<int, int> predicate_histogram
  (
   const Entity& entity,
   size_t which_structure) {
    std::map<int, int> hist;
    for(const minsky::Mention& mention : entity.mentions) {
      const minsky::PredArg& pa = mention.structures[which_structure];
      hist[pa.predicate.word] += 1;
    }
    return hist;
  }
  std::map<int, int> predicate_histogram_on_level
  (
   const Entity& entity,
   const AnnotationLevel::type& al
   ) {
    std::map<int, int> hist;
    bool unhit = true;
    for(const minsky::Mention& mention : entity.mentions) {
      for(const minsky::PredArg& pa : mention.structures) {
	if(pa.annot_level != al) continue;
	hist[pa.predicate.word] += 1;
	unhit = false;
      }
    }
    if(unhit) {
      WARN << "None of the " << entity.mentions.size() << " mentions of entity " << entity.id << " had a structure level of " << al;
    }
    return hist;
  }
  std::map<int, int> relation_histogram
  (
   const Entity& entity,
   size_t which_structure
   ) {
    std::map<int, int> hist;
    for(const minsky::Mention& mention : entity.mentions) {
      const minsky::PredArg& pa = mention.structures[which_structure];
      hist[pa.relation] += 1;
    }
    return hist;
  }
  std::map<int, int> relation_histogram_on_level
  (
   const Entity& entity,
   const AnnotationLevel::type& al
   ) {
    std::map<int, int> hist;
    bool unhit = true;
    for(const minsky::Mention& mention : entity.mentions) {
      for(const minsky::PredArg& pa : mention.structures) {
	if(pa.annot_level != al) continue;
	hist[pa.relation] += 1;
	unhit = false;
      }
    }
    if(unhit) {
      WARN << "None of the " << entity.mentions.size() << " mentions of entity " << entity.id << " had a relation structure level of " << al;
    }
    return hist;
  }

  std::map<int, int> predicate_histogram
  (
   const EDoc& doc,
   const AnnotationLevel::type& al
   ) {
    std::map<int, int> hist;
    for(const minsky::Entity& entity : doc.entities) {
      for(const minsky::Mention& mention : entity.mentions) {
	for(const minsky::PredArg& pa : mention.structures) {
	  if(pa.annot_level == al) {
	    hist[pa.predicate.word] += 1;
	  }
	}
      }
    }
    return hist;
  }
  std::map<int, int> relation_histogram
  (
   const EDoc& doc,
   const AnnotationLevel::type& al
   ) {
    std::map<int, int> hist;
    for(const minsky::Entity& entity : doc.entities) {
      for(const minsky::Mention& mention : entity.mentions) {
	for(const minsky::PredArg& pa : mention.structures) {
	  if(pa.annot_level == al) {
	    hist[pa.relation] += 1;
	  }
	}
      }
    }
    return hist;
  }

  int find_level(const Mention& ment, const AnnotationLevel::type& al, bool print) {
    const size_t s = ment.structures.size();
    int ret = -1;
    switch(s) {
    case 0:
      if(print)
	ERROR << "No set structures in mention " << ment.id;
      ret = -1;
      break;
    case 2:
      if(ment.structures[1].annot_level == al) {
	ret = 1;
	goto finally;
      }
      // deliberate fallthrough
    case 1:
      if(ment.structures[0].annot_level == al) {
	ret = 0;
	goto finally;
      }
      else goto notfound;
      break;
    default:
      ret = 0;
      for(const auto& pa : ment.structures) {
	if(pa.annot_level == al) {
	  goto finally;
	}
	++ret;
      }
      goto notfound;
    }
  notfound:
    if(print)
      ERROR << "Did not find any annotation level " << al << " in mention " << ment.id;
    ret = -1;
  finally:
    return ret;
  }

  minsky::CountList& operator+=(minsky::CountList& recept, const minsky::CountList& other) {
    for(const auto& opair : other.icounts) {
      recept.icounts[opair.first] += opair.second;
    }
    recept.__isset.icounts = true;
    return recept;
  }
  bool empty(const minsky::WordsClause& wc) {
    if(wc.__isset.words) {
      return wc.words.empty();
    }
    if(wc.__isset.counts) {
      if(wc.counts.__isset.icounts) {
	return wc.counts.icounts.empty();
      }
      return true;
    }
    return true;
  }

  minsky::CountList bow_counts(const SimpleDoc& doc) {
    minsky::CountList cl;
    for(const minsky::WordsClause& sentence : doc.sentences) {
      if(sentence.__isset.counts) {
	cl += sentence.counts;
      } else if(sentence.__isset.words) {
	for(int word : sentence.words) {
	  cl.icounts[word] += 1;
	}
      }
    }
    return cl;
  }

  size_t num_words(const SimpleDoc& doc) {
    size_t n{};
    for(const minsky::WordsClause& sentence : doc.sentences) {
      n += num_words(sentence.counts);
    }
    return n;
  }
  size_t num_words(const CountList& counts) {
    size_t n{};
    for(const auto& pair : counts.icounts) {
      n += pair.second;
    }
    return n;
  }

  SyntacticEDocVocabUpdater::SyntacticEDocVocabUpdater() :
    MinskyDocVocabUpdater<SyntacticEDocVocabUpdater, minsky::EDoc>() {
  }
  void SyntacticEDocVocabUpdater::operator()(minsky::EDoc& doc, const std::vector<std::vector<int> >& vocab_mappers) {
    if(vocab_mappers.size() != 2) {
      ERROR << __PRETTY_FUNCTION__ << " can only have 2 vocab mappers passed in, not " << vocab_mappers.size();
      return;
    }
    const std::vector<int>& pred_map = vocab_mappers[0];
    INFO << "Remapping observables in " << doc.id;
    INFO << "Pred mapper knows of " << pred_map.size() << " predicates";
    const std::vector<int>& rel_map = vocab_mappers[1];
    INFO << "Relation mapper knows of " << rel_map.size() << " relations";
    for(Entity& ent : doc.entities) {
      for(Mention& mention : ent.mentions) {
	for(PredArg& pa : mention.structures) {
	  if(pa.annot_level != minsky::AnnotationLevel::SYNTAX) {
	    continue;
	  }
	  if(pa.__isset.predicate && pa.predicate.__isset.word) {
	    int orig = pa.predicate.word;
	    INFO << "Original: " << orig << ", out of " << pred_map.size();
	    int other = pred_map.at(orig);
	    INFO << "Transformed p: " << orig << " --> " << other;
	    pa.predicate.word = other;
	  } else {
	    WARN << "Predicate is not set for mention " << mention.id << " in document " << doc.id;
	  }
	  if(pa.__isset.relation) {
	    int orig = pa.relation;
	    INFO << "Original: " << orig << ", out of " << rel_map.size();
	    int other = rel_map.at(orig);
	    INFO << "Transformed r: " << orig << " --> " << other;
	    pa.relation = other;
	  } else {
	    WARN << "Relation is not set for mention " << mention.id << " in document " << doc.id;
	  }
	}
      }
    }
  }

  PredArgEDocVocabUpdater::PredArgEDocVocabUpdater(const std::vector<minsky::MinskyAnnotationWrapper>& order) :
    MinskyDocVocabUpdater<PredArgEDocVocabUpdater, minsky::EDoc>() {
    int i = 0;
    for(const auto& p : order) {
      orders_[p] = (i++);
    }
  }
  void PredArgEDocVocabUpdater::operator()(minsky::EDoc& doc, const std::vector<std::vector<int> >& vocab_mappers) {
    if(vocab_mappers.size() != orders_.size()) {
      ERROR << __PRETTY_FUNCTION__ << " must have " << orders_.size() << " vocab mappers passed in";
      return;
    }
    INFO << "Remapping observables in " << doc.id;
    for(const auto& pk : orders_) {
      minsky::MinskyAnnotationWrapper p = pk.first;
      INFO << minsky::_AnnotationLevel_VALUES_TO_NAMES.at(p.annot_level) << " " << minsky::_StructureType_VALUES_TO_NAMES.at(p.st_level) << " mapper knows of " << vocab_mappers[pk.second].size()  << minsky::_StructureType_VALUES_TO_NAMES.at(p.st_level) << "s";
    }
    for(Entity& ent : doc.entities) {
      for(Mention& mention : ent.mentions) {
	for(PredArg& pa : mention.structures) {
	  minsky::MinskyAnnotationWrapper maw;
	  maw.annot_level = pa.annot_level;
	  if(pa.__isset.predicate && pa.predicate.__isset.word) {
	    maw.st_level = minsky::StructureType::PREDICATE;
	    auto it = orders_.find(maw);
	    if(it == orders_.end()) continue;
	    const std::vector<int>& w_map = vocab_mappers[it->second];
	    int orig = pa.predicate.word;
	    INFO << "Original: " << orig << ", out of " << w_map.size();
	    int other = w_map.at(orig);
	    INFO << "Transformed p: " << orig << " --> " << other;
	    pa.predicate.word = other;
	  }
	  if(pa.__isset.relation) {
	    maw.st_level = minsky::StructureType::RELATION;
	    auto it = orders_.find(maw);
	    if(it == orders_.end()) continue;
	    const std::vector<int>& w_map = vocab_mappers[it->second];
	    ////
	    int orig = pa.relation;
	    INFO << "Original: " << orig << ", out of " << w_map.size();
	    int other = w_map.at(orig);
	    INFO << "Transformed r: " << orig << " --> " << other;
	    pa.relation = other;
	  }
	}
      }
    }
  }

  SimpleDocVocabUpdater::SimpleDocVocabUpdater() :
    MinskyDocVocabUpdater<SimpleDocVocabUpdater, minsky::SimpleDoc>(){
  }
  void SimpleDocVocabUpdater::operator()(minsky::SimpleDoc& doc, const std::vector<std::vector<int> >& vocab_mappers) {
    if(vocab_mappers.size() != 1) {
      ERROR << __PRETTY_FUNCTION__ << " can only have 1 vocab mapper passed in, not " << vocab_mappers.size();
      return;
    }
    const std::vector<int>& vmap = vocab_mappers[0];
    INFO << "Word mapper knows of " << vmap.size() << " predicates";
    for(WordsClause& wc : doc.sentences) {
      if(wc.__isset.counts) {
	std::map<int, int> nm;
	for(const std::pair<int, int> pc : wc.counts.icounts) {
	  nm[ vmap.at(pc.first) ] = pc.second;
	}
	wc.counts.__set_icounts(nm);
      } else if(wc.__isset.words) {
	const size_t s = wc.words.size();
	for(size_t i = 0; i < s; ++i) {
	  const int old = wc.words[i];
	  wc.words[i] = vmap.at(old);
	}
      }
    }
  }

  bool LabelCmp::operator()(const minsky::Label& lhs, const minsky::Label& rhs) const {
    if(lhs.label < rhs.label) return true;
    else if(lhs.label > rhs.label) return false;
    else return lhs.weight < rhs.weight;
  }

#ifdef THRIFT091
  std::ostream& operator<<(std::ostream& os, const minsky::Label& label) {
    if(label.__isset.label) {
      os << label.label;
      if(label.__isset.weight) {
  	os << "(" << label.weight << ")";
      }
    }
    return os;
  }
#endif

  std::string get_label(const minsky::SimpleDoc& doc) {
    if(doc.labels.size() == 0) {
      return "No label";
    }
    std::stringstream ss;
    const size_t L = doc.labels.size();
    size_t l = 0;
    for(const auto& label : doc.labels) {
      ss << label;
      ++l;
      if(l < L) {
	ss << ", ";
      }
    }
    return ss.str();
  }
} // ends namespace minsky

namespace ferrum {
  std::map<int, int> doc_multinomial(const minsky::EDoc& doc,
				     const minsky::MinskyAnnotationWrapper& maw) {
				     //const minsky::AnnotationLevel::type& annot_level,
				     //const minsky::StructureType::type& st_level) {
    switch(maw.st_level) {
    case minsky::StructureType::PREDICATE:
      return minsky::predicate_histogram(doc, maw.annot_level);
    case minsky::StructureType::RELATION:
      return minsky::relation_histogram(doc, maw.annot_level);
    default:
      ERROR << "Unable to handle StructureLevel " << maw.st_level << " for document " << doc.id;
      throw 10;
    }
  }

  std::map<int, int> doc_multinomial(const minsky::SimpleDoc& doc,
				     const minsky::MinskyAnnotationWrapper& maw) {
				     //const minsky::AnnotationLevel::type& annot_level,
				     //const minsky::StructureType::type& st_level) {
    switch(maw.word_level) {
    case minsky::WordAnnotation::ORTHOGRAPHIC:
    case minsky::WordAnnotation::LEMMA:
      {
      minsky::CountList cl = minsky::bow_counts(doc);
      return cl.icounts;
      }
      break;
    default:
      ERROR << "Unable to handle StructureLevel " << maw.word_level << " for document " << doc.id;
      throw 10;
    }
  }

  const std::vector<minsky::Entity>& get_entities(const minsky::EDoc& doc) {
    return doc.entities;
  }
  const std::vector<minsky::Mention>& get_mentions(const minsky::Entity& entity) {
    return entity.mentions;
  }
  int num_entities(const minsky::EDoc& doc) {
    return (int)doc.entities.size();
  }
  int num_entities(const minsky::Entity& entity) {
    return (int)entity.mentions.size();
  }
  int get_num_mentions(const minsky::EDoc& doc) {
    int nm = 0;
    for(const auto& ent : doc.entities) {
      nm += ent.mentions.size();
    }
    return nm;
  }
}
