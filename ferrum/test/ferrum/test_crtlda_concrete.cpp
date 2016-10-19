#include "gtest/gtest.h"

#include "ferrum/minsky.hpp"
#include "ferrum/crtlda_concrete.hpp"
#include "ferrum/crtlda_pruner_minsky.hpp"
#include "ferrum/logging.hpp"


#include <iostream>

TEST(ConcreteCRTLDAIntegration, get_archived_corpus_compact) {
  typedef ferrum::InMemoryCorpus<minsky::EDoc> CorpusT;
  typedef ferrum::MinskySituationGovernedPruner PrunerT;
  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  CorpusT* corpus = NULL;
  ASSERT_NO_THROW({
      // NOTE: GoogleTest is sensitive to commas in templates
      corpus = (ferrum::get_archived_corpus<CorpusT, PrunerT>
	(
	 "test_corpus",
	 "resources/AFP_ENG_19940531.0390.tcompact.2copies.tar.gz",
	 gov_voc,
	 rel_voc 
	 ));
	}
    );
  ASSERT_TRUE(corpus != NULL);
  ASSERT_EQ(2, corpus->num_docs());
  delete corpus;
}

TEST(ConcreteCRTLDAIntegration, get_archived_corpus_compact2) {
  typedef ferrum::InMemoryCorpus<minsky::EDoc> CorpusT;
  typedef ferrum::MinskySituationGovernedPruner PrunerT;
  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  CorpusT* corpus = NULL;
  ASSERT_NO_THROW({
      // NOTE: GoogleTest is sensitive to commas in templates
      corpus = (ferrum::get_archived_corpus<CorpusT, PrunerT>
	(
	 "test_corpus",
	 "resources/twenty_nyt_semafor_tcompact-v4.tar.gz",
	 gov_voc,
	 rel_voc 
	 ));
	}
    );
  ASSERT_TRUE(corpus != NULL);
  ASSERT_EQ(20, corpus->num_docs());
  std::vector<std::string> expected_ids = {
    "NYT_ENG_19950614.0461",
    "NYT_ENG_19950913.0259",
    "NYT_ENG_19951003.0428",
    "NYT_ENG_19951219.0207",
    "NYT_ENG_19980515.0355",
    "NYT_ENG_19990217.0032",
    "NYT_ENG_19990419.0195",
    "NYT_ENG_19991214.0407",
    "NYT_ENG_20000305.0062",
    "NYT_ENG_20010202.0043",
    "NYT_ENG_20010806.0019",
    "NYT_ENG_20041105.0019",
    "NYT_ENG_20050121.0174",
    "NYT_ENG_20050828.0022",
    "NYT_ENG_20060803.0060",
    "NYT_ENG_20061019.0280",
    "NYT_ENG_20070629.0031",
    "NYT_ENG_20071108.0195",
    "NYT_ENG_20080531.0057",
    "NYT_ENG_20080813.0177"
  };
  for(size_t i = 0; i < 20; ++i) {
    EXPECT_EQ(expected_ids[i], corpus->operator[](i).id) << "document " << i << " (0-indexed) doesn't match";
  }
  delete corpus;
}

