#ifndef PUBLIC_UPF
#include "gtest/gtest.h"

#include "concrete_util/io.h"
#include "ferrum/minsky.hpp"
#include "ferrum/crtlda_minsky.hpp"
#include "ferrum/crtlda.hpp"
#include "ferrum/crtlda_pruner_minsky.hpp"
#include "ferrum/sage_tm_variational.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/data_util.hpp"
#include "ferrum/redis_corpus.hpp"
#include "ferrum/crtlda_concrete.hpp"

TEST(SAGETopicModel, from_twenty_nyt_semafor_tcompact_v4_tar_gz) {
  const char *name = "resources/twenty_nyt_semafor_tcompact-v4.tar.gz";
  typedef std::string string;
  ferrum::Vocabulary<string> vocab("__WORD_OOV__");
  typedef minsky::SimpleDoc Doc;

  ferrum::db::Address addr("localhost", 4532);
  std::shared_ptr<ferrum::db::Redis> main_redis_connection(new ferrum::db::Redis(addr));
  std::string corpus_name(__FUNCTION__);
  corpus_name += std::string(name);
  ferrum::RedisCorpus<Doc>* corpus =
    new ferrum::RedisCorpus<Doc>(corpus_name, main_redis_connection);
  ferrum::Toolnames tools;
  minsky::WordAnnotation::type word_form = minsky::WordAnnotation::ORTHOGRAPHIC;
  ferrum::set_archived_corpus<ferrum::RedisCorpus<Doc>,
			      ferrum::MinskyDocBOWPruner, // the pruner type
			      concrete::util::TCompactProtocol // how to read the Communications
			      >
    (
     name,
     corpus,
     &vocab,
     tools,
     word_form,
     0 // number of dependency hops
     );
  for(auto it = corpus->begin(); it != corpus->end(); ++it) {
    const Doc& doc = *(it->document);
    EXPECT_EQ(1, doc.sentences.size());
    EXPECT_TRUE(doc.sentences[0].counts.__isset.icounts);
    EXPECT_GT(doc.sentences[0].counts.icounts.size(), 0);
    for(const auto& pair : doc.sentences[0].counts.icounts) {
      DEBUG << vocab.word(pair.first) << " occurs " << pair.second << " times";
    }
  }

  std::shared_ptr<std::vector<double> > background(new std::vector<double>(vocab.num_words(), 0.0));

  INFO << "Computing word background model...";      
  ferrum::RedisBackgroundComputer<Doc, ferrum::Vocabulary<std::string> > rbc;
  ferrum::BDocBackgroundCounter<ferrum::RedisCorpus<Doc>, ferrum::RedisBackgroundComputer<Doc, ferrum::Vocabulary<std::string> > > bec;
  bec.in_memory_as(background)
    .defined_by(&vocab)
    .with(&rbc)
    .over(corpus)
    .compute_background<std::string>();
  INFO << "... done computing word background model...";

  typedef ferrum::SageTMVariational Inferencer;
  const int num_topics = 2;
  Inferencer var_inf(num_topics, vocab.num_words(), 0.1, 0.1);
  ferrum::StringDiscreteKindPrinter inf_printer;
  inf_printer.g_vocab = &vocab;

  INFO << "Initializing training model";
  const int num_train_tokens = ferrum::num_words(corpus);
  ferrum::SageInitializer initializer(num_train_tokens, 1);
  initializer.sparse(true, true, 1E-7);
  var_inf.init(initializer, vocab, corpus, background);
  INFO << "... done initializing training model";

  ferrum::DKVWriters sw_wrapper;
  ferrum::VStrategy strategy;
  strategy.batch_size = 5;
  for(int epoch = 0; epoch < 20; ++epoch) {
  var_inf.learn< ferrum::RedisCorpus<Doc>,
		 ::ferrum::db::RedisThriftSmartWriter,
		 decltype(main_redis_connection) >
    (
     corpus,
     strategy,
     epoch,
     inf_printer,
     &sw_wrapper,
     main_redis_connection
     );
  }

  minsky::residual::ResidualTopicModel rgs =
    var_inf.create_minsky_view(vocab);
  ASSERT_EQ(num_topics, rgs.topics.size());
  for(int i = 0; i < num_topics; ++i) {
    ASSERT_TRUE(rgs.topics[i].__isset.frame);
    ASSERT_TRUE(rgs.topics[i].frame.__isset.distr);
    ASSERT_TRUE(rgs.topics[i].frame.distr.weights.__isset.residual);
    std::vector<size_t> srt = ferrum::sort_indices(rgs.topics[i].frame.distr.weights.residual, false);
    for(int j = 0; j < 20; ++j) {
      EXPECT_GT(rgs.topics[i].frame.distr.weights.residual[srt[j]], 1E-6);
      INFO << "Topic " << i << ", word #" << j << ": " << vocab.word(srt[j]) << " : " << rgs.topics[i].frame.distr.weights.residual[srt[j]];
    }
  }
}

TEST(SAGETopicModel, reset_buffers) {
  const char *name = "resources/twenty_nyt_semafor_tcompact-v4.tar.gz";
  typedef std::string string;
  ferrum::Vocabulary<string> vocab("__WORD_OOV__");
  typedef minsky::SimpleDoc Doc;

  ferrum::db::Address addr("localhost", 4532);
  std::shared_ptr<ferrum::db::Redis> main_redis_connection(new ferrum::db::Redis(addr));
  std::string corpus_name(__FUNCTION__);
  corpus_name += std::string(name);
  ferrum::RedisCorpus<Doc>* corpus =
    new ferrum::RedisCorpus<Doc>(corpus_name, main_redis_connection);
  ferrum::Toolnames tools;
  minsky::WordAnnotation::type word_form = minsky::WordAnnotation::ORTHOGRAPHIC;
  ferrum::set_archived_corpus<ferrum::RedisCorpus<Doc>,
			      ferrum::MinskyDocBOWPruner, // the pruner type
			      concrete::util::TCompactProtocol // how to read the Communications
			      >
    (
     name,
     corpus,
     &vocab,
     tools,
     word_form,
     0 // number of dependency hops
     );
  for(auto it = corpus->begin(); it != corpus->end(); ++it) {
    const Doc& doc = *(it->document);
    EXPECT_EQ(1, doc.sentences.size());
    EXPECT_TRUE(doc.sentences[0].counts.__isset.icounts);
    EXPECT_GT(doc.sentences[0].counts.icounts.size(), 0);
    for(const auto& pair : doc.sentences[0].counts.icounts) {
      DEBUG << vocab.word(pair.first) << " occurs " << pair.second << " times";
    }
  }

  std::shared_ptr<std::vector<double> > background(new std::vector<double>(vocab.num_words(), 0.0));

  INFO << "Computing word background model...";      
  ferrum::RedisBackgroundComputer<Doc, ferrum::Vocabulary<std::string> > rbc;
  ferrum::BDocBackgroundCounter<ferrum::RedisCorpus<Doc>, ferrum::RedisBackgroundComputer<Doc, ferrum::Vocabulary<std::string> > > bec;
  bec.in_memory_as(background)
    .defined_by(&vocab)
    .with(&rbc)
    .over(corpus)
    .compute_background<std::string>();
  INFO << "... done computing word background model...";

  typedef ferrum::SageTMVariational Inferencer;
  const int num_topics = 2;
  Inferencer var_inf(num_topics, vocab.num_words(), 0.1, 0.1);
  ferrum::StringDiscreteKindPrinter inf_printer;
  inf_printer.g_vocab = &vocab;

  const int num_train_tokens = ferrum::num_words(corpus);
  INFO << "Initializing training model with " << num_train_tokens << " tokens";
  ferrum::SageInitializer initializer(num_train_tokens, 1);
  initializer.sparse(true, true, 1E-7);
  var_inf.init(initializer, vocab, corpus, background);
  INFO << "... done initializing training model";

  ferrum::DKVWriters sw_wrapper;
  ferrum::VStrategy strategy;
  strategy.num_e_iters = 1;
  for(int epoch = 0; epoch < 20; ++epoch) {
    for(size_t t = 0; t < (size_t)num_topics; ++t) {
      double t_sum = ferrum::sum(var_inf.buffer_topic_word_params(t));
      ASSERT_LT(t_sum, 2.0*num_train_tokens) << "Failed for topic " << t << " on epoch " << epoch;
    }
    var_inf.learn< ferrum::RedisCorpus<Doc>,
		   ::ferrum::db::RedisThriftSmartWriter,
		   decltype(main_redis_connection) >
      (
       corpus,
       strategy,
       epoch,
       inf_printer,
       &sw_wrapper,
       main_redis_connection
       );
  }
}
#endif
