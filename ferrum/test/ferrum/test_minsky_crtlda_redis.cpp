#include "gtest/gtest.h"

#include "concrete_util/io.h"
#include "ferrum/minsky.hpp"
#include "ferrum/crtlda.hpp"
#include "ferrum/crtlda_concrete.hpp"
#include "ferrum/crtlda_minsky.hpp"
#include "ferrum/crtlda_pruner_minsky.hpp"
#include "ferrum/data_util.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/redis_corpus.hpp"

#include <iostream>

TEST(MinskySituationGovernedPruner, from_compact_NYT_ENG_19980113_0597_comm) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/NYT_ENG_19980113.0597.tcompact.with_semafor.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("NYT_ENG_19980113.0597", comm.id);

  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  minsky::EDoc edoc_full =
    ferrum::MinskySituationGovernedPruner::make(comm, &gov_voc, &rel_voc);
  ASSERT_EQ(40, minsky::num_entities(edoc_full));

  ferrum::RedisCorpus<minsky::EDoc> rc({"localhost", 4532});
  //ferrum::RedisCorpus<minsky::EDoc> rc;
  ASSERT_EQ(0, rc.num_docs());
  rc.add_document(edoc_full);
  ASSERT_EQ(1, rc.num_docs());
  int num_from_iterator = 0;
  for(ferrum::RedisCorpus<minsky::EDoc>::const_iterator it = rc.begin();
      it != rc.end();
      ++it) {
    const minsky::EDoc& doc = *(it->document);
    INFO << (it->iteration_idx) << " :: " << doc.id;
    ++num_from_iterator;
  }
  ASSERT_GT(ferrum::num_mentions(&rc), 0);
}

TEST(MinskySituationGovernedPruner, from_tgz) {
  
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  std::string name("resources/twenty_nyt_semafor_tcompact-v4.tar.gz");
  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  typedef minsky::EDoc Doc;
  ferrum::RedisCorpus<Doc>* rc =
    ferrum::get_archived_corpus<ferrum::RedisCorpus<Doc>,
				ferrum::MinskySituationGovernedPruner,
				concrete::util::TCompactProtocol,
				std::string,
				ferrum::db::Address>
    (
     name,
     gov_voc,
     rel_voc,
     "train_corpus",
    {"localhost", 4532}
     );
  ASSERT_TRUE(rc != NULL);
  ASSERT_EQ(20, rc->num_docs());
  int num_from_iterator = 0;
  for(ferrum::RedisCorpus<minsky::EDoc>::const_iterator it = rc->begin();
      it != rc->end();
      ++it) {
    const minsky::EDoc& doc = *(it->document);
    INFO << (it->iteration_idx) << " :: " << doc.id;
    ++num_from_iterator;
  }
  delete rc;
}

TEST(RedisCorpus, save_load_vocab) {
  typedef std::string string;
  ferrum::Vocabulary<string> vocab("my_oov");
  vocab.make_word("foo");
  vocab.make_word("bar");
  ferrum::RedisCorpus<minsky::EDoc> rc({"localhost", 4532});
  ASSERT_NO_THROW( ( rc.save_vocab("sample_corpus", __func__, vocab) ) );
  ferrum::Vocabulary<string> recon;
  ASSERT_NO_THROW( ( rc.load_vocab("sample_corpus", __func__, recon) ) );
  ASSERT_TRUE( ( vocab == recon ) );
}

TEST(RedisCorpus, unify_vocabs) {
  typedef std::string string;
  ferrum::Vocabulary<string> base_vocab("my_oov");
  base_vocab.make_word("foo");
  base_vocab.make_word("bar");
  ferrum::Vocabulary<string> base_rel_vocab("rel_oov");
  base_rel_vocab.make_word("rel1");
  std::string pf(__PRETTY_FUNCTION__);
  std::string corpus_name(pf + ":0");
  ferrum::db::Address addr("localhost", 4532);
  std::shared_ptr<ferrum::db::Redis> redis(new ferrum::db::Redis(addr));
  ferrum::RedisCorpus<minsky::EDoc> rc(corpus_name, redis);
  {
    minsky::EDoc doc;
    doc.__set_id("first_doc");
    minsky::Entity entity;
    entity.__set_id("first_entity");
    minsky::Mention mention;
    mention.__set_id("first_mention");
    minsky::PredArg syntactic;
    minsky::RelationFiller pred_rf;
    pred_rf.__set_word(base_vocab("foo"));
    syntactic.__set_predicate(pred_rf);
    syntactic.__set_relation(base_rel_vocab("rel1"));
    syntactic.__set_annot_level(minsky::AnnotationLevel::SYNTAX);
    mention.structures.push_back(syntactic);
    mention.__isset.structures = true;
    entity.mentions.push_back(mention);
    entity.__isset.mentions = true;
    doc.entities.push_back(entity);
    doc.__isset.entities = true;
    // push the doc into the corpus
    rc.add_document(doc);
  }
  ASSERT_EQ(1, rc.num_docs());
  {
    const minsky::EDoc& stored_doc = rc[0];
    ASSERT_EQ(1, stored_doc.entities[0].mentions[0].structures[0].predicate.word);
    ASSERT_EQ(1, stored_doc.entities[0].mentions[0].structures[0].relation);
  }
  std::vector<std::string> vocab_names = {"gov_vocab", "rel_vocab"};
  std::vector<ferrum::Vocabulary<std::string>* > vocabs = {&base_vocab, &base_rel_vocab};
  minsky::SyntacticEDocVocabUpdater dm; // document mapper
  ASSERT_NO_THROW(rc.unify_vocabs("unified", vocab_names, vocabs, dm));
  ASSERT_EQ(1, rc.num_docs());
  ferrum::db::RedisQuery num_unified_g_words("gov_vocab:unified", "num_words");
  num_unified_g_words.hget();
  redis->operator()(num_unified_g_words);
  ASSERT_EQ("3", num_unified_g_words.value());
  {
    const minsky::EDoc& stored_doc = rc[0];
    ASSERT_EQ(1, stored_doc.entities[0].mentions[0].structures[0].predicate.word);
    ASSERT_EQ(1, stored_doc.entities[0].mentions[0].structures[0].relation);
  }

  ferrum::Vocabulary<string> base_vocab1("my_oov1");
  base_vocab1.make_word("bar");
  base_vocab1.make_word("foo");
  ferrum::Vocabulary<string> base_rel_vocab1("rel_oov1");
  base_rel_vocab1.make_word("rel2");
  base_rel_vocab1.make_word("rel3");
  base_rel_vocab1.make_word("rel1");

  std::string corpus_name1(pf + ":1");
  std::shared_ptr<ferrum::db::Redis> redis1(new ferrum::db::Redis(addr));
  ferrum::RedisCorpus<minsky::EDoc> rc1(corpus_name1, redis1);
  {
    minsky::EDoc doc;
    doc.__set_id("second_doc");
    minsky::Entity entity;
    entity.__set_id("first_entity");
    minsky::Mention mention;
    mention.__set_id("first_mention");
    minsky::PredArg syntactic;
    minsky::RelationFiller pred_rf;
    pred_rf.__set_word(base_vocab1("foo"));
    syntactic.__set_predicate(pred_rf);
    syntactic.__set_relation(base_rel_vocab1("rel1"));
    syntactic.__set_annot_level(minsky::AnnotationLevel::SYNTAX);
    mention.structures.push_back(syntactic);
    mention.__isset.structures = true;
    entity.mentions.push_back(mention);
    entity.__isset.mentions = true;
    doc.entities.push_back(entity);
    doc.__isset.entities = true;
    // push the doc into the corpus
    rc1.add_document(doc);
  }
  ASSERT_EQ(1, rc1.num_docs());

  vocabs[0] = &base_vocab1;
  vocabs[1] = &base_rel_vocab1;
  // check that rc1 has stored the docs wrt the base*vocab1 values
  {
    const minsky::EDoc& stored_doc = rc1[0];
    ASSERT_EQ(2, stored_doc.entities[0].mentions[0].structures[0].predicate.word);
    ASSERT_EQ(3, stored_doc.entities[0].mentions[0].structures[0].relation);
  }
  ASSERT_NO_THROW(rc1.unify_vocabs("unified", vocab_names, vocabs, dm));
  ASSERT_EQ(1, rc.num_docs());
  ASSERT_EQ(1, rc1.num_docs());
  // now check that rc1 has stored the docs wrt the base*vocab values
  {
    const minsky::EDoc& stored_doc = rc1[0];
    ASSERT_EQ(1, stored_doc.entities[0].mentions[0].structures[0].predicate.word);
    ASSERT_EQ(1, stored_doc.entities[0].mentions[0].structures[0].relation);
  }
  ferrum::db::RedisQuery num_unified_g_words1("gov_vocab:unified", "num_words");
  num_unified_g_words1.hget();
  redis1->operator()(num_unified_g_words1);
  ASSERT_EQ("3", num_unified_g_words1.value());
  ferrum::db::RedisQuery num_unified_r_words1("rel_vocab:unified", "num_words");
  num_unified_r_words1.hget();
  redis1->operator()(num_unified_r_words1);
  ASSERT_EQ("4", num_unified_r_words1.value());
}

TEST(RedisCorpus, subset_0) {
  std::string name("resources/twenty_nyt_semafor_tcompact-v4.tar.gz");
  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  typedef minsky::EDoc Doc;
  ferrum::RedisCorpus<Doc>* rc =
    ferrum::get_archived_corpus<ferrum::RedisCorpus<Doc>,
				ferrum::MinskySituationGovernedPruner,
				concrete::util::TCompactProtocol,
				std::string,
				ferrum::db::Address>
    (
     name,
     gov_voc,
     rel_voc,
     "train_corpus",
    {"localhost", 4532}
     );
  ASSERT_TRUE(rc != NULL);
  ASSERT_EQ(20, rc->num_docs());
  int num_from_iterator = 0;
  {
    ferrum::RedisCorpus<Doc> subset5_10 =
      rc->subset(5, 10);
    ASSERT_EQ(5, subset5_10.num_docs());
    const std::vector<std::string>& dnames = rc->doc_names();
    ASSERT_EQ(20, dnames.size());
    const std::vector<std::string>& sub_dnames = subset5_10.doc_names();
    ASSERT_EQ(5, sub_dnames.size());
    for(size_t i = 5; i < 10; ++i) {
      ASSERT_EQ(dnames[i], sub_dnames[i - 5]);
    }
    for(ferrum::RedisCorpus<minsky::EDoc>::const_iterator it = subset5_10.begin();
	it != subset5_10.end();
	++it) {
      const minsky::EDoc& doc = *(it->document);
      ASSERT_EQ(dnames[ it->iteration_idx + 5], doc.id);
      ++num_from_iterator;
    }
    ASSERT_EQ(5, num_from_iterator);
  }
  {
    ferrum::RedisCorpus<Doc> subset17_20 =
      rc->subset(17, 24); // 24 intentional here
    const size_t expected_size = 3;
    ASSERT_EQ(expected_size, subset17_20.num_docs());
    const std::vector<std::string>& dnames = rc->doc_names();
    ASSERT_EQ(20, dnames.size());
    const std::vector<std::string>& sub_dnames = subset17_20.doc_names();
    ASSERT_EQ(expected_size, sub_dnames.size());
    for(size_t i = 17; i < 20; ++i) {
      ASSERT_EQ(dnames[i], sub_dnames[i - 17]);
    }
    num_from_iterator = 0;
    for(ferrum::RedisCorpus<minsky::EDoc>::const_iterator it = subset17_20.begin();
	it != subset17_20.end();
	++it) {
      const minsky::EDoc& doc = *(it->document);
      ASSERT_EQ(dnames[ it->iteration_idx + 17], doc.id);
      ++num_from_iterator;
    }
    ASSERT_EQ(expected_size, num_from_iterator);
  }
  delete rc;
}

TEST(RedisCorpus, subset_0_mt) {
  std::string name("resources/twenty_nyt_semafor_tcompact-v4.tar.gz");
  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  typedef minsky::EDoc Doc;
  ferrum::RedisCorpus<Doc>* rc =
    ferrum::get_archived_corpus<ferrum::RedisCorpus<Doc>,
				ferrum::MinskySituationGovernedPruner,
				concrete::util::TCompactProtocol,
				std::string,
				ferrum::db::Address>
    (
     name,
     gov_voc,
     rel_voc,
     "train_corpus",
    {"localhost", 4532}
     );
  ASSERT_TRUE(rc != NULL);
  ASSERT_EQ(20, rc->num_docs());
  rc->multithreaded(true);
#pragma omp parallel for num_threads(2)
  for(int thread_i = 0; thread_i < 4; ++thread_i) {
    int num_from_iterator = 0;
    {
      ferrum::RedisCorpus<Doc> subset5_10 =
	rc->subset(5, 10);
      ASSERT_EQ(5, subset5_10.num_docs());
      const std::vector<std::string>& dnames = rc->doc_names();
      ASSERT_EQ(20, dnames.size());
      const std::vector<std::string>& sub_dnames = subset5_10.doc_names();
      ASSERT_EQ(5, sub_dnames.size());
      for(size_t i = 5; i < 10; ++i) {
	ASSERT_EQ(dnames[i], sub_dnames[i - 5]);
      }
      for(ferrum::RedisCorpus<minsky::EDoc>::const_iterator it = subset5_10.begin();
	  it != subset5_10.end();
	  ++it) {
	const minsky::EDoc& doc = *(it->document);
	ASSERT_EQ(dnames[ it->iteration_idx + 5], doc.id);
	++num_from_iterator;
      }
      ASSERT_EQ(5, num_from_iterator);
    }
    {
      ferrum::RedisCorpus<Doc> subset17_20 =
	rc->subset(17, 24); // 24 intentional here
      const size_t expected_size = 3;
      ASSERT_EQ(expected_size, subset17_20.num_docs());
      const std::vector<std::string>& dnames = rc->doc_names();
      ASSERT_EQ(20, dnames.size());
      const std::vector<std::string>& sub_dnames = subset17_20.doc_names();
      ASSERT_EQ(expected_size, sub_dnames.size());
      for(size_t i = 17; i < 20; ++i) {
	ASSERT_EQ(dnames[i], sub_dnames[i - 17]);
      }
      num_from_iterator = 0;
      for(ferrum::RedisCorpus<minsky::EDoc>::const_iterator it = subset17_20.begin();
	  it != subset17_20.end();
	  ++it) {
	const minsky::EDoc& doc = *(it->document);
	ASSERT_EQ(dnames[ it->iteration_idx + 17], doc.id);
	++num_from_iterator;
      }
      ASSERT_EQ(expected_size, num_from_iterator);
    }
  }
  delete rc;
}

TEST(MinskySituationGovernedPruner, background_computation) {
  std::string name("resources/twenty_nyt_semafor_tcompact-v4.tar.gz");
  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  typedef minsky::EDoc Doc;
  ferrum::RedisCorpus<Doc>* rc =
    ferrum::get_archived_corpus<ferrum::RedisCorpus<Doc>,
				ferrum::MinskySituationGovernedPruner,
				concrete::util::TCompactProtocol,
				std::string,
				ferrum::db::Address>
    (
     name,
     gov_voc,
     rel_voc,
     "train_corpus",
    {"localhost", 4532}
     );
  ASSERT_TRUE(rc != NULL);
  ASSERT_EQ(20, rc->num_docs());

  ferrum::MinskyEntityCounter mec(minsky::AnnotationLevel::SYNTAX);
  ferrum::RedisBackgroundComputer<Doc, ferrum::Vocabulary<std::string> > rbc;
  ferrum::EDocBackgroundCounter<ferrum::RedisCorpus<Doc>, ferrum::RedisBackgroundComputer<Doc, ferrum::Vocabulary<std::string> > > bec;

  std::shared_ptr<std::vector<double> > gov_background(new std::vector<double>(gov_voc.num_words(), 0.0));

  bec.in_memory_as(gov_background)
    .defined_by(&gov_voc)
    .how(&mec)
    .with(&rbc)
    .over(rc)
    .compute_predicate_background<std::string>();

  delete rc;
}

TEST(MinskyDocBOWPruner, background_computation) {
  std::string name("resources/twenty_nyt_semafor_tcompact-v4.tar.gz");
  typedef std::string string;
  ferrum::Vocabulary<string> voc;
  typedef minsky::SimpleDoc Doc;
  ferrum::Toolnames tools;
  // ferrum::RedisCorpus<Doc>* rc = new ferrum::RedisCorpus<Doc>("train_simpledoc_corpus_ortho",
  //   {"localhost", 4532});
  std::shared_ptr<ferrum::db::Redis> db(new ferrum::db::Redis({"localhost", 4532}));
  ferrum::RedisCorpus<Doc>* rc =
    ferrum::get_archived_corpus<ferrum::RedisCorpus<Doc>,
  				ferrum::MinskyDocBOWPruner, // the pruner type
  				concrete::util::TCompactProtocol, // how to read the Communications
  				std::string, // RedisCorpus ctor args...
  				std::shared_ptr<ferrum::db::Redis> >
    (
     name,
     voc,
     tools,
     minsky::WordAnnotation::ORTHOGRAPHIC,
     "train_simpledoc_corpus_ortho",
     db
     );
  ASSERT_TRUE(rc != NULL);
  ASSERT_EQ(20, rc->num_docs());

  ferrum::RedisBackgroundComputer<Doc, ferrum::Vocabulary<std::string> > rbc;
  ferrum::BDocBackgroundCounter<ferrum::RedisCorpus<Doc>, ferrum::RedisBackgroundComputer<Doc, ferrum::Vocabulary<std::string> > > bec;

  std::shared_ptr<std::vector<double> > background(new std::vector<double>(voc.num_words(), 0.0));

  bec.in_memory_as(background)
    .defined_by(&voc)
    .with(&rbc)
    .over(rc)
    .compute_background<std::string>();

  delete rc;
}

TEST(MinskyVerbBOWPruner, from_tgz) {
  std::string name("resources/twenty_nyt_semafor_tcompact-v4.tar.gz");
  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  typedef minsky::SimpleDoc Doc;
  ferrum::RedisCorpus<Doc>* rc = new ferrum::RedisCorpus<Doc>("train_corpus", ferrum::db::Address("localhost", 4532));
  ferrum::set_archived_corpus<ferrum::RedisCorpus<Doc>,
			      ferrum::MinskyVerbBOWPruner,
			      concrete::util::TCompactProtocol>
    (
     name,
     rc,
     &gov_voc
     );
  ASSERT_TRUE(rc != NULL);
  ASSERT_EQ(20, rc->num_docs());
  int num_from_iterator = 0;
  for(ferrum::RedisCorpus<minsky::SimpleDoc>::const_iterator it = rc->begin();
      it != rc->end();
      ++it) {
    const minsky::SimpleDoc& doc = *(it->document);
    INFO << (it->iteration_idx) << " :: " << doc.id;
    ++num_from_iterator;
  }
  delete rc;
}
