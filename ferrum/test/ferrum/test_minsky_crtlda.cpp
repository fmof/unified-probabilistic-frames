#include "gtest/gtest.h"

#include "concrete_util/io.h"
#include "ferrum/minsky.hpp"
#include "ferrum/crtlda_minsky.hpp"
#include "ferrum/crtlda.hpp"
#include "ferrum/crtlda_pruner_minsky.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/data_util.hpp"

#include <iostream>

TEST(MinskyVerbGovernedPruner, from_compact_NYT_ENG_19980113_0597_comm) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/NYT_ENG_19980113.0597.tcompact.with_semafor.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("NYT_ENG_19980113.0597", comm.id);

  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  const ferrum::MinskyVerbGovernedPruner sgp(comm, &gov_voc, &rel_voc);
  concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention =
    concrete_util::mention_id_to_mention(comm, "Stanford");
  const concrete::UUID& mention_id = comm.entityMentionSetList[0].mentionList[0].uuid;
  ASSERT_NE(mention_id_to_mention.find(mention_id), mention_id_to_mention.end());
  const concrete::EntityMention& em = mention_id_to_mention[mention_id];
  concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization =
    concrete_util::mention_id_to_tokenization(comm, "Stanford");
  ASSERT_EQ(1, mention_id_to_tokenization.count(mention_id));
  ASSERT_NE(mention_id_to_tokenization.find(mention_id), mention_id_to_tokenization.end());
  const concrete::Tokenization& tokenization = mention_id_to_tokenization.at(mention_id);
  DEBUG << "Mention id: " << mention_id.uuidString;
  ASSERT_EQ("beae7914-9205-420c-b762-7a4bfa8d388c", em.uuid.uuidString);
  ASSERT_EQ("2fc3e7e9-e103-4791-a82e-1adb525f647d", em.tokens.tokenizationId.uuidString);
  DEBUG << "Trying to inspect tokenization " << tokenization.uuid.uuidString ;
  DEBUG << "Tokenization is at " << &tokenization;
  ASSERT_TRUE(tokenization.__isset.dependencyParseList) << "Tokenization " << tokenization.uuid.uuidString << " does not have a dependency parse list set";
  const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(tokenization, "col-ccproc-deps");
  ASSERT_TRUE(dep_parse != NULL);
  std::vector< minsky::Mention > created_mentions = sgp.prune(em);
  ASSERT_NO_THROW(created_mentions.size());
  //////////////////////////////////////////////////
  minsky::EDoc edoc_full =
    ferrum::MinskyVerbGovernedPruner::make(comm, &gov_voc, &rel_voc);
  ASSERT_EQ(70, minsky::num_entities(edoc_full));
  //////////////////////////////////////////////////
  minsky::EDoc copy = minsky::EDoc(edoc_full);
  ASSERT_TRUE( edoc_full == copy );
  //////////////////////////////////////////////////
  std::string estring = ferrum::thrift::thrift_struct_to_string<ferrum::thrift::TCompactProtocol>( edoc_full );
  minsky::EDoc reconstructed;
  ferrum::thrift::thrift_struct_from_string<ferrum::thrift::TCompactProtocol>
    (
     estring,
     &reconstructed
     );
  ASSERT_TRUE( edoc_full == reconstructed );
  ////////////////////////////////////////////////
  std::vector<double> counts(gov_voc.num_words(), 0.0);
  for(const auto& entity : ferrum::get_entities(edoc_full) ) {
    for(const auto& mention : entity.mentions) {
      for(const minsky::PredArg& pa : mention.structures) {
	ASSERT_TRUE( const_cast<minsky::PredArg*>(&pa)->__isset.annot_level );
	ASSERT_NE(pa.annot_level, minsky::AnnotationLevel::UNSPECIFIED);
	if(pa.annot_level == minsky::AnnotationLevel::SYNTAX) {
	  counts[pa.predicate.word] += 1;
	  DEBUG << "Seeing predicate " << pa.predicate.word << " (count = " << counts[pa.predicate.word] << ")";
	}
      }
    }
  }

  std::vector<double> counts0(gov_voc.num_words(), 0.0);
  std::vector<double> counts1(gov_voc.num_words(), 0.0);
  ferrum::MinskyEntityCounter mec(minsky::AnnotationLevel::SYNTAX);
  int num_hit = 0;
  for(const auto& entity : ferrum::get_entities(edoc_full) ) {
    ferrum::TemplatedEntityInterface< decltype(entity) > tei = mec.make_entity_interface< decltype(entity) >(entity);
    //mec.update_predicate_counts(&gov_voc, &tei, counts0);
    tei.update_p_count(&mec, &gov_voc, counts0);
    for(const auto& mention : tei.entity.mentions) {
      ++num_hit;
      ASSERT_GT(num_hit, 0);
      for(const minsky::PredArg& pa : mention.structures) {
	ASSERT_TRUE( const_cast<minsky::PredArg*>(&pa)->__isset.annot_level );
	ASSERT_NE(pa.annot_level, minsky::AnnotationLevel::UNSPECIFIED);
	if(pa.annot_level == minsky::AnnotationLevel::SYNTAX) {
	  counts1[pa.predicate.word] += 1;
	}
      }
    }
  }

  for(size_t i = 0; i < (size_t)gov_voc.num_words(); ++i) {
    ASSERT_EQ(counts[i], counts0[i]);
    ASSERT_EQ(counts[i], counts1[i]);
  }
}

TEST(MinskySituationGovernedPruner, from_compact_NYT_ENG_19980113_0597_comm) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/NYT_ENG_19980113.0597.tcompact.with_semafor.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("NYT_ENG_19980113.0597", comm.id);

  typedef std::string string;
  ferrum::Vocabulary<string> gov_voc;
  ferrum::Vocabulary<string> rel_voc;
  const ferrum::MinskySituationGovernedPruner sgp(comm, &gov_voc, &rel_voc);
  concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention =
    concrete_util::mention_id_to_mention(comm, "Stanford");
  const concrete::UUID& mention_id = comm.entityMentionSetList[0].mentionList[0].uuid;
  ASSERT_NE(mention_id_to_mention.find(mention_id), mention_id_to_mention.end());
  const concrete::EntityMention& em = mention_id_to_mention[mention_id];
  concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization =
    concrete_util::mention_id_to_tokenization(comm, "Stanford");
  ASSERT_EQ(1, mention_id_to_tokenization.count(mention_id));
  ASSERT_NE(mention_id_to_tokenization.find(mention_id), mention_id_to_tokenization.end());
  const concrete::Tokenization& tokenization = mention_id_to_tokenization.at(mention_id);
  DEBUG << "Mention id: " << mention_id.uuidString;
  ASSERT_EQ("beae7914-9205-420c-b762-7a4bfa8d388c", em.uuid.uuidString);
  ASSERT_EQ("2fc3e7e9-e103-4791-a82e-1adb525f647d", em.tokens.tokenizationId.uuidString);
  DEBUG << "Trying to inspect tokenization " << tokenization.uuid.uuidString ;
  DEBUG << "Tokenization is at " << &tokenization;
  ASSERT_TRUE(tokenization.__isset.dependencyParseList) << "Tokenization " << tokenization.uuid.uuidString << " does not have a dependency parse list set";
  const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(tokenization, "col-ccproc-deps");
  ASSERT_TRUE(dep_parse != NULL);
  std::vector< minsky::Mention > created_mentions = sgp.prune(em);
  ASSERT_NO_THROW(created_mentions.size());
  //////////////////////////////////////////////////
  minsky::EDoc edoc_full =
    ferrum::MinskySituationGovernedPruner::make(comm, &gov_voc, &rel_voc);
  ASSERT_EQ(40, minsky::num_entities(edoc_full));
  //////////////////////////////////////////////////
  minsky::EDoc copy = minsky::EDoc(edoc_full);
  ASSERT_TRUE( edoc_full == copy );
  //////////////////////////////////////////////////
  std::string estring = ferrum::thrift::thrift_struct_to_string<ferrum::thrift::TCompactProtocol>( edoc_full );
  minsky::EDoc reconstructed;
  ferrum::thrift::thrift_struct_from_string<ferrum::thrift::TCompactProtocol>
    (
     estring,
     &reconstructed
     );
  ASSERT_TRUE( edoc_full == reconstructed );
  ////////////////////////////////////////////////
  std::vector<double> counts(gov_voc.num_words(), 0.0);
  for(const auto& entity : ferrum::get_entities(edoc_full) ) {
    for(const auto& mention : entity.mentions) {
      for(const minsky::PredArg& pa : mention.structures) {
	ASSERT_TRUE( const_cast<minsky::PredArg*>(&pa)->__isset.annot_level );
	ASSERT_NE(pa.annot_level, minsky::AnnotationLevel::UNSPECIFIED);
	if(pa.annot_level == minsky::AnnotationLevel::SYNTAX) {
	  counts[pa.predicate.word] += 1;
	  DEBUG << "Seeing predicate " << pa.predicate.word << ", " << gov_voc.word(pa.predicate.word)  << ", (count = " << counts[pa.predicate.word] << ")";
	}
      }
    }
  }

  std::vector<double> counts0(gov_voc.num_words(), 0.0);
  std::vector<double> counts1(gov_voc.num_words(), 0.0);
  ferrum::MinskyEntityCounter mec(minsky::AnnotationLevel::SYNTAX);
  int num_hit = 0;
  for(const auto& entity : ferrum::get_entities(edoc_full) ) {
    ferrum::TemplatedEntityInterface< decltype(entity) > tei = mec.make_entity_interface< decltype(entity) >(entity);
    //mec.update_predicate_counts(&gov_voc, &tei, counts0);
    tei.update_p_count(&mec, &gov_voc, counts0);
    for(const auto& mention : tei.entity.mentions) {
      ++num_hit;
      ASSERT_GT(num_hit, 0);
      for(const minsky::PredArg& pa : mention.structures) {
	ASSERT_TRUE( const_cast<minsky::PredArg*>(&pa)->__isset.annot_level );
	ASSERT_NE(pa.annot_level, minsky::AnnotationLevel::UNSPECIFIED);
	if(pa.annot_level == minsky::AnnotationLevel::SYNTAX) {
	  counts1[pa.predicate.word] += 1;
	}
      }
    }
  }

  for(size_t i = 0; i < (size_t)gov_voc.num_words(); ++i) {
    ASSERT_EQ(counts[i], counts0[i]);
    ASSERT_EQ(counts[i], counts1[i]);
  }
}

void _MinskyDocBOWPruner_from_compact_NYT_ENG_19980113_0597_comm_token_check(const std::map<int,int>&, const ferrum::Vocabulary<std::string>&);
void _MinskyDocBOWPruner_from_compact_NYT_ENG_19980113_0597_comm_token_count(const std::map<int,int>&, const ferrum::Vocabulary<std::string>&);

TEST(MinskyDocBOWPruner, from_compact_NYT_ENG_19980113_0597_comm) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/NYT_ENG_19980113.0597.tcompact.with_semafor.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("NYT_ENG_19980113.0597", comm.id);

  typedef std::string string;
  std::vector<size_t> sizes = {279};
  for(auto form : {minsky::WordAnnotation::LEMMA, minsky::WordAnnotation::ORTHOGRAPHIC} ) {
    ferrum::Vocabulary<string> voc("__OOV__");
    minsky::SimpleDoc doc = ferrum::MinskySentencePruner<ferrum::MinskyDocBOWPruner, concrete::Communication>::make(comm, &voc, form);
    ASSERT_EQ("NYT_ENG_19980113.0597", doc.id);
    ASSERT_TRUE(doc.__isset.sentences);
    ASSERT_EQ(1, doc.sentences.size());
  }
  int i = 0;
  for(auto form : {minsky::WordAnnotation::ORTHOGRAPHIC} ) {
    ferrum::Vocabulary<string> voc("__OOV__");
    minsky::SimpleDoc doc = ferrum::MinskySentencePruner<ferrum::MinskyDocBOWPruner, concrete::Communication>::make(comm, &voc, form);
    ASSERT_EQ("NYT_ENG_19980113.0597", doc.id);
    ASSERT_TRUE(doc.__isset.sentences);
    ASSERT_EQ(1, doc.sentences.size());
    ASSERT_EQ(voc.num_words(), sizes[i] + 1);
    ASSERT_EQ(doc.sentences[0].counts.icounts.size(), sizes[i]);
    _MinskyDocBOWPruner_from_compact_NYT_ENG_19980113_0597_comm_token_check(doc.sentences[0].counts.icounts, voc);
    _MinskyDocBOWPruner_from_compact_NYT_ENG_19980113_0597_comm_token_count(doc.sentences[0].counts.icounts, voc);
    ++i;
  }
  // const ferrum::MinskyDocBOWPruner sgp(comm, &voc);

}

/*
 * $ concrete_inspect.py --dependency test/resources/NYT_ENG_19980113.0597.tcompact.with_semafor.concrete| cut -f2 | tail -n+3 | perl -nle 'print lc' | grep -vP '^$' | sort | uniq | perl -ne 'chomp; print("EXPECT_TRUE(m.find(voc.index(\"$_\")) != m.end());\n")'
 */
void _MinskyDocBOWPruner_from_compact_NYT_ENG_19980113_0597_comm_token_check(const std::map<int,int>& m, const ferrum::Vocabulary<std::string>& voc) {
  EXPECT_TRUE(m.find(voc.index("``")) != m.end());
  EXPECT_TRUE(m.find(voc.index(",")) != m.end());
  EXPECT_TRUE(m.find(voc.index(".")) != m.end());
  EXPECT_TRUE(m.find(voc.index("''")) != m.end());
  EXPECT_TRUE(m.find(voc.index("$")) != m.end());
  EXPECT_TRUE(m.find(voc.index("10")) != m.end());
  EXPECT_TRUE(m.find(voc.index("1.19")) != m.end());
  EXPECT_TRUE(m.find(voc.index("12")) != m.end());
  EXPECT_TRUE(m.find(voc.index("13")) != m.end());
  EXPECT_TRUE(m.find(voc.index("14")) != m.end());
  EXPECT_TRUE(m.find(voc.index("1.5")) != m.end());
  EXPECT_TRUE(m.find(voc.index("16,600")) != m.end());
  EXPECT_TRUE(m.find(voc.index("17")) != m.end());
  EXPECT_TRUE(m.find(voc.index("18")) != m.end());
  EXPECT_TRUE(m.find(voc.index("1996")) != m.end());
  EXPECT_TRUE(m.find(voc.index("1997")) != m.end());
  EXPECT_TRUE(m.find(voc.index("1 3/8")) != m.end());
  EXPECT_TRUE(m.find(voc.index("25")) != m.end());
  EXPECT_TRUE(m.find(voc.index("3")) != m.end());
  EXPECT_TRUE(m.find(voc.index("34")) != m.end());
  EXPECT_TRUE(m.find(voc.index("3.93")) != m.end());
  EXPECT_TRUE(m.find(voc.index("400")) != m.end());
  EXPECT_TRUE(m.find(voc.index("4th-qtr")) != m.end());
  EXPECT_TRUE(m.find(voc.index("61 5/16")) != m.end());
  EXPECT_TRUE(m.find(voc.index("65")) != m.end());
  EXPECT_TRUE(m.find(voc.index("7:30")) != m.end());
  EXPECT_TRUE(m.find(voc.index("75")) != m.end());
  EXPECT_TRUE(m.find(voc.index("8.8")) != m.end());
  EXPECT_TRUE(m.find(voc.index("a")) != m.end());
  EXPECT_TRUE(m.find(voc.index("about")) != m.end());
  EXPECT_TRUE(m.find(voc.index("addition")) != m.end());
  EXPECT_TRUE(m.find(voc.index("adjust")) != m.end());
  EXPECT_TRUE(m.find(voc.index("against")) != m.end());
  EXPECT_TRUE(m.find(voc.index("alex")) != m.end());
  EXPECT_TRUE(m.find(voc.index("already")) != m.end());
  EXPECT_TRUE(m.find(voc.index("a.m.")) != m.end());
  EXPECT_TRUE(m.find(voc.index("amid")) != m.end());
  EXPECT_TRUE(m.find(voc.index("an")) != m.end());
  EXPECT_TRUE(m.find(voc.index("analyst")) != m.end());
  EXPECT_TRUE(m.find(voc.index("analysts")) != m.end());
  EXPECT_TRUE(m.find(voc.index("and")) != m.end());
  EXPECT_TRUE(m.find(voc.index("are")) != m.end());
  EXPECT_TRUE(m.find(voc.index("as")) != m.end());
  EXPECT_TRUE(m.find(voc.index("assault")) != m.end());
  EXPECT_TRUE(m.find(voc.index("at")) != m.end());
  EXPECT_TRUE(m.find(voc.index("attached")) != m.end());
  EXPECT_TRUE(m.find(voc.index("average")) != m.end());
  EXPECT_TRUE(m.find(voc.index("b.")) != m.end());
  EXPECT_TRUE(m.find(voc.index("back")) != m.end());
  EXPECT_TRUE(m.find(voc.index("barney")) != m.end());
  EXPECT_TRUE(m.find(voc.index("based")) != m.end());
  EXPECT_TRUE(m.find(voc.index("bc-kodak-forecast-update1-bloom")) != m.end());
  EXPECT_TRUE(m.find(voc.index("be")) != m.end());
  EXPECT_TRUE(m.find(voc.index("beat")) != m.end());
  EXPECT_TRUE(m.find(voc.index("billion")) != m.end());
  EXPECT_TRUE(m.find(voc.index("bloomberg")) != m.end());
  EXPECT_TRUE(m.find(voc.index("both")) != m.end());
  EXPECT_TRUE(m.find(voc.index("business")) != m.end());
  EXPECT_TRUE(m.find(voc.index("by")) != m.end());
  EXPECT_TRUE(m.find(voc.index("category")) != m.end());
  EXPECT_TRUE(m.find(voc.index("cause")) != m.end());
  EXPECT_TRUE(m.find(voc.index("cents")) != m.end());
  EXPECT_TRUE(m.find(voc.index("change")) != m.end());
  EXPECT_TRUE(m.find(voc.index("charge")) != m.end());
  EXPECT_TRUE(m.find(voc.index("charles")) != m.end());
  EXPECT_TRUE(m.find(voc.index("co.")) != m.end());
  EXPECT_TRUE(m.find(voc.index("color")) != m.end());
  EXPECT_TRUE(m.find(voc.index("comment")) != m.end());
  EXPECT_TRUE(m.find(voc.index("compact")) != m.end());
  EXPECT_TRUE(m.find(voc.index("companies")) != m.end());
  EXPECT_TRUE(m.find(voc.index("company")) != m.end());
  EXPECT_TRUE(m.find(voc.index("competition")) != m.end());
  EXPECT_TRUE(m.find(voc.index("computers")) != m.end());
  EXPECT_TRUE(m.find(voc.index("consumer")) != m.end());
  EXPECT_TRUE(m.find(voc.index("continued")) != m.end());
  EXPECT_TRUE(m.find(voc.index("currency")) != m.end());
  EXPECT_TRUE(m.find(voc.index("cut")) != m.end());
  EXPECT_TRUE(m.find(voc.index("decline")) != m.end());
  EXPECT_TRUE(m.find(voc.index("declined")) != m.end());
  EXPECT_TRUE(m.find(voc.index("degree")) != m.end());
  EXPECT_TRUE(m.find(voc.index("deutsche")) != m.end());
  EXPECT_TRUE(m.find(voc.index("digital")) != m.end());
  EXPECT_TRUE(m.find(voc.index("diluted")) != m.end());
  EXPECT_TRUE(m.find(voc.index("disclosed")) != m.end());
  EXPECT_TRUE(m.find(voc.index("discount")) != m.end());
  EXPECT_TRUE(m.find(voc.index("discs")) != m.end());
  EXPECT_TRUE(m.find(voc.index("does")) != m.end());
  EXPECT_TRUE(m.find(voc.index("dollar")) != m.end());
  EXPECT_TRUE(m.find(voc.index("dollars")) != m.end());
  EXPECT_TRUE(m.find(voc.index("dominates")) != m.end());
  EXPECT_TRUE(m.find(voc.index("double-digit")) != m.end());
  EXPECT_TRUE(m.find(voc.index("down")) != m.end());
  EXPECT_TRUE(m.find(voc.index("earned")) != m.end());
  EXPECT_TRUE(m.find(voc.index("earnings")) != m.end());
  EXPECT_TRUE(m.find(voc.index("eastman")) != m.end());
  EXPECT_TRUE(m.find(voc.index("effort")) != m.end());
  EXPECT_TRUE(m.find(voc.index("established")) != m.end());
  EXPECT_TRUE(m.find(voc.index("estimate")) != m.end());
  EXPECT_TRUE(m.find(voc.index("evidence")) != m.end());
  EXPECT_TRUE(m.find(voc.index("expect")) != m.end());
  EXPECT_TRUE(m.find(voc.index("expected")) != m.end());
  EXPECT_TRUE(m.find(voc.index("extent")) != m.end());
  EXPECT_TRUE(m.find(voc.index("fall")) != m.end());
  EXPECT_TRUE(m.find(voc.index("fell")) != m.end());
  EXPECT_TRUE(m.find(voc.index("fewer")) != m.end());
  EXPECT_TRUE(m.find(voc.index("film")) != m.end());
  EXPECT_TRUE(m.find(voc.index("financial")) != m.end());
  EXPECT_TRUE(m.find(voc.index("first")) != m.end());
  EXPECT_TRUE(m.find(voc.index("flexibility")) != m.end());
  EXPECT_TRUE(m.find(voc.index("focus")) != m.end());
  EXPECT_TRUE(m.find(voc.index("for")) != m.end());
  EXPECT_TRUE(m.find(voc.index("fourth-quarter")) != m.end());
  EXPECT_TRUE(m.find(voc.index("from")) != m.end());
  EXPECT_TRUE(m.find(voc.index("fuji")) != m.end());
  EXPECT_TRUE(m.find(voc.index("gain")) != m.end());
  EXPECT_TRUE(m.find(voc.index("getting")) != m.end());
  EXPECT_TRUE(m.find(voc.index("going")) != m.end());
  EXPECT_TRUE(m.find(voc.index("gold")) != m.end());
  EXPECT_TRUE(m.find(voc.index("growing")) != m.end());
  EXPECT_TRUE(m.find(voc.index("had")) != m.end());
  EXPECT_TRUE(m.find(voc.index("has")) != m.end());
  EXPECT_TRUE(m.find(voc.index("have")) != m.end());
  EXPECT_TRUE(m.find(voc.index("he")) != m.end());
  EXPECT_TRUE(m.find(voc.index("henderson")) != m.end());
  EXPECT_TRUE(m.find(voc.index("high-profit")) != m.end());
  EXPECT_TRUE(m.find(voc.index("his")) != m.end());
  EXPECT_TRUE(m.find(voc.index("hit")) != m.end());
  EXPECT_TRUE(m.find(voc.index("hold")) != m.end());
  EXPECT_TRUE(m.find(voc.index("holds")) != m.end());
  EXPECT_TRUE(m.find(voc.index("holiday")) != m.end());
  EXPECT_TRUE(m.find(voc.index("holidays")) != m.end());
  EXPECT_TRUE(m.find(voc.index("hurts")) != m.end());
  EXPECT_TRUE(m.find(voc.index("ibes")) != m.end());
  EXPECT_TRUE(m.find(voc.index("imaging")) != m.end());
  EXPECT_TRUE(m.find(voc.index("in")) != m.end());
  EXPECT_TRUE(m.find(voc.index("inc.")) != m.end());
  EXPECT_TRUE(m.find(voc.index("include")) != m.end());
  EXPECT_TRUE(m.find(voc.index("includes")) != m.end());
  EXPECT_TRUE(m.find(voc.index("international")) != m.end());
  EXPECT_TRUE(m.find(voc.index("is")) != m.end());
  EXPECT_TRUE(m.find(voc.index("issued")) != m.end());
  EXPECT_TRUE(m.find(voc.index("it")) != m.end());
  EXPECT_TRUE(m.find(voc.index("its")) != m.end());
  EXPECT_TRUE(m.find(voc.index("jan.")) != m.end());
  EXPECT_TRUE(m.find(voc.index("japan")) != m.end());
  EXPECT_TRUE(m.find(voc.index("jobs")) != m.end());
  EXPECT_TRUE(m.find(voc.index("jonathan")) != m.end());
  EXPECT_TRUE(m.find(voc.index("just")) != m.end());
  EXPECT_TRUE(m.find(voc.index("key")) != m.end());
  EXPECT_TRUE(m.find(voc.index("kodak")) != m.end());
  EXPECT_TRUE(m.find(voc.index("largest")) != m.end());
  EXPECT_TRUE(m.find(voc.index("last")) != m.end());
  EXPECT_TRUE(m.find(voc.index("local")) != m.end());
  EXPECT_TRUE(m.find(voc.index("looking")) != m.end());
  EXPECT_TRUE(m.find(voc.index("loss")) != m.end());
  EXPECT_TRUE(m.find(voc.index("losses")) != m.end());
  EXPECT_TRUE(m.find(voc.index("lost")) != m.end());
  EXPECT_TRUE(m.find(voc.index("lower")) != m.end());
  EXPECT_TRUE(m.find(voc.index("-lrb-")) != m.end());
  EXPECT_TRUE(m.find(voc.index("made")) != m.end());
  EXPECT_TRUE(m.find(voc.index("margins")) != m.end());
  EXPECT_TRUE(m.find(voc.index("mark")) != m.end());
  EXPECT_TRUE(m.find(voc.index("market")) != m.end());
  EXPECT_TRUE(m.find(voc.index("market-share")) != m.end());
  EXPECT_TRUE(m.find(voc.index("microfilm")) != m.end());
  EXPECT_TRUE(m.find(voc.index("million")) != m.end());
  EXPECT_TRUE(m.find(voc.index("month")) != m.end());
  EXPECT_TRUE(m.find(voc.index("more")) != m.end());
  EXPECT_TRUE(m.find(voc.index("mounting")) != m.end());
  EXPECT_TRUE(m.find(voc.index("much")) != m.end());
  EXPECT_TRUE(m.find(voc.index("need")) != m.end());
  EXPECT_TRUE(m.find(voc.index("negative")) != m.end());
  EXPECT_TRUE(m.find(voc.index("neutral")) != m.end());
  EXPECT_TRUE(m.find(voc.index("new")) != m.end());
  EXPECT_TRUE(m.find(voc.index("n't")) != m.end());
  EXPECT_TRUE(m.find(voc.index("of")) != m.end());
  EXPECT_TRUE(m.find(voc.index("older")) != m.end());
  EXPECT_TRUE(m.find(voc.index("on")) != m.end());
  EXPECT_TRUE(m.find(voc.index("operating")) != m.end());
  EXPECT_TRUE(m.find(voc.index("or")) != m.end());
  EXPECT_TRUE(m.find(voc.index("others")) != m.end());
  EXPECT_TRUE(m.find(voc.index("overseas")) != m.end());
  EXPECT_TRUE(m.find(voc.index("percent")) != m.end());
  EXPECT_TRUE(m.find(voc.index("percentage")) != m.end());
  EXPECT_TRUE(m.find(voc.index("photo")) != m.end());
  EXPECT_TRUE(m.find(voc.index("photography")) != m.end());
  EXPECT_TRUE(m.find(voc.index("plans")) != m.end());
  EXPECT_TRUE(m.find(voc.index("point")) != m.end());
  EXPECT_TRUE(m.find(voc.index("pressure")) != m.end());
  EXPECT_TRUE(m.find(voc.index("price")) != m.end());
  EXPECT_TRUE(m.find(voc.index("prices")) != m.end());
  EXPECT_TRUE(m.find(voc.index("products")) != m.end());
  EXPECT_TRUE(m.find(voc.index("profit")) != m.end());
  EXPECT_TRUE(m.find(voc.index("promotion")) != m.end());
  EXPECT_TRUE(m.find(voc.index("prudential")) != m.end());
  EXPECT_TRUE(m.find(voc.index("quarter")) != m.end());
  EXPECT_TRUE(m.find(voc.index("quite")) != m.end());
  EXPECT_TRUE(m.find(voc.index("rates")) != m.end());
  EXPECT_TRUE(m.find(voc.index("rating")) != m.end());
  EXPECT_TRUE(m.find(voc.index("'re")) != m.end());
  EXPECT_TRUE(m.find(voc.index("release")) != m.end());
  EXPECT_TRUE(m.find(voc.index("replace")) != m.end());
  EXPECT_TRUE(m.find(voc.index("report")) != m.end());
  EXPECT_TRUE(m.find(voc.index("respond")) != m.end());
  EXPECT_TRUE(m.find(voc.index("result")) != m.end());
  EXPECT_TRUE(m.find(voc.index("rises")) != m.end());
  EXPECT_TRUE(m.find(voc.index("rising")) != m.end());
  EXPECT_TRUE(m.find(voc.index("rochester")) != m.end());
  EXPECT_TRUE(m.find(voc.index("rose")) != m.end());
  EXPECT_TRUE(m.find(voc.index("rosenzweig")) != m.end());
  EXPECT_TRUE(m.find(voc.index("-rrb-")) != m.end());
  EXPECT_TRUE(m.find(voc.index("'s")) != m.end());
  EXPECT_TRUE(m.find(voc.index("said")) != m.end());
  EXPECT_TRUE(m.find(voc.index("sales")) != m.end());
  EXPECT_TRUE(m.find(voc.index("salomon")) != m.end());
  EXPECT_TRUE(m.find(voc.index("season")) != m.end());
  EXPECT_TRUE(m.find(voc.index("second")) != m.end());
  EXPECT_TRUE(m.find(voc.index("securities")) != m.end());
  EXPECT_TRUE(m.find(voc.index("see")) != m.end());
  EXPECT_TRUE(m.find(voc.index("seen")) != m.end());
  EXPECT_TRUE(m.find(voc.index("september")) != m.end());
  EXPECT_TRUE(m.find(voc.index("share")) != m.end());
  EXPECT_TRUE(m.find(voc.index("shares")) != m.end());
  EXPECT_TRUE(m.find(voc.index("shopping")) != m.end());
  EXPECT_TRUE(m.find(voc.index("since")) != m.end());
  EXPECT_TRUE(m.find(voc.index("slowing")) != m.end());
  EXPECT_TRUE(m.find(voc.index("slumping")) != m.end());
  EXPECT_TRUE(m.find(voc.index("smith")) != m.end());
  EXPECT_TRUE(m.find(voc.index("socked")) != m.end());
  EXPECT_TRUE(m.find(voc.index("sold")) != m.end());
  EXPECT_TRUE(m.find(voc.index("some")) != m.end());
  EXPECT_TRUE(m.find(voc.index("spends")) != m.end());
  EXPECT_TRUE(m.find(voc.index("spokesman")) != m.end());
  EXPECT_TRUE(m.find(voc.index("steal")) != m.end());
  EXPECT_TRUE(m.find(voc.index("stem")) != m.end());
  EXPECT_TRUE(m.find(voc.index("stock")) != m.end());
  EXPECT_TRUE(m.find(voc.index("storage")) != m.end());
  EXPECT_TRUE(m.find(voc.index("such")) != m.end());
  EXPECT_TRUE(m.find(voc.index("surveyed")) != m.end());
  EXPECT_TRUE(m.find(voc.index("systems")) != m.end());
  EXPECT_TRUE(m.find(voc.index("take")) != m.end());
  EXPECT_TRUE(m.find(voc.index("than")) != m.end());
  EXPECT_TRUE(m.find(voc.index("that")) != m.end());
  EXPECT_TRUE(m.find(voc.index("the")) != m.end());
  EXPECT_TRUE(m.find(voc.index("they")) != m.end());
  EXPECT_TRUE(m.find(voc.index("thing")) != m.end());
  EXPECT_TRUE(m.find(voc.index("third")) != m.end());
  EXPECT_TRUE(m.find(voc.index("through")) != m.end());
  EXPECT_TRUE(m.find(voc.index("thursday")) != m.end());
  EXPECT_TRUE(m.find(voc.index("time")) != m.end());
  EXPECT_TRUE(m.find(voc.index("to")) != m.end());
  EXPECT_TRUE(m.find(voc.index("today")) != m.end());
  EXPECT_TRUE(m.find(voc.index("triple")) != m.end());
  EXPECT_TRUE(m.find(voc.index("tumbling")) != m.end());
  EXPECT_TRUE(m.find(voc.index("two")) != m.end());
  EXPECT_TRUE(m.find(voc.index("unchanged")) != m.end());
  EXPECT_TRUE(m.find(voc.index("update1")) != m.end());
  EXPECT_TRUE(m.find(voc.index("u.s.")) != m.end());
  EXPECT_TRUE(m.find(voc.index("u.s.-based")) != m.end());
  EXPECT_TRUE(m.find(voc.index("using")) != m.end());
  EXPECT_TRUE(m.find(voc.index("version")) != m.end());
  EXPECT_TRUE(m.find(voc.index("we")) != m.end());
  EXPECT_TRUE(m.find(voc.index("weak")) != m.end());
  EXPECT_TRUE(m.find(voc.index("were")) != m.end());
  EXPECT_TRUE(m.find(voc.index("what")) != m.end());
  EXPECT_TRUE(m.find(voc.index("which")) != m.end());
  EXPECT_TRUE(m.find(voc.index("while")) != m.end());
  EXPECT_TRUE(m.find(voc.index("who")) != m.end());
  EXPECT_TRUE(m.find(voc.index("widening")) != m.end());
  EXPECT_TRUE(m.find(voc.index("wider")) != m.end());
  EXPECT_TRUE(m.find(voc.index("will")) != m.end());
  EXPECT_TRUE(m.find(voc.index("workforce")) != m.end());
  EXPECT_TRUE(m.find(voc.index("world")) != m.end());
  EXPECT_TRUE(m.find(voc.index("would")) != m.end());
  EXPECT_TRUE(m.find(voc.index("wrote")) != m.end());
  EXPECT_TRUE(m.find(voc.index("year")) != m.end());
  EXPECT_TRUE(m.find(voc.index("year-earlier")) != m.end());
  EXPECT_TRUE(m.find(voc.index("york")) != m.end());
  EXPECT_TRUE(m.find(voc.index("york-based")) != m.end());
}

/*
 * $ concrete_inspect.py --dependency test/resources/NYT_ENG_19980113.0597.tcompact.with_semafor.concrete| cut -f2 | tail -n+3 | perl -nle 'print lc' | grep -vP '^$' | sort | uniq -c | perl -ne 'chomp; /^\s*(\d+)\s+(\S+)$/; print("EXPECT_EQ(m.at(voc.index(\"$2\")), $1);\n")'
 */
void _MinskyDocBOWPruner_from_compact_NYT_ENG_19980113_0597_comm_token_count(const std::map<int,int>& m, const ferrum::Vocabulary<std::string>& voc) {
  EXPECT_EQ(m.at(voc.index("``")), 7);
  EXPECT_EQ(m.at(voc.index(",")), 33);
  EXPECT_EQ(m.at(voc.index(".")), 24);
  EXPECT_EQ(m.at(voc.index("''")), 9);
  EXPECT_EQ(m.at(voc.index("$")), 4);
  EXPECT_EQ(m.at(voc.index("10")), 1);
  EXPECT_EQ(m.at(voc.index("1.19")), 1);
  EXPECT_EQ(m.at(voc.index("12")), 1);
  EXPECT_EQ(m.at(voc.index("13")), 1);
  EXPECT_EQ(m.at(voc.index("14")), 1);
  EXPECT_EQ(m.at(voc.index("1.5")), 1);
  EXPECT_EQ(m.at(voc.index("16,600")), 1);
  EXPECT_EQ(m.at(voc.index("17")), 1);
  EXPECT_EQ(m.at(voc.index("18")), 1);
  EXPECT_EQ(m.at(voc.index("1996")), 1);
  EXPECT_EQ(m.at(voc.index("1997")), 2);
  EXPECT_EQ(m.at(voc.index("1 3/8")), 1);
  EXPECT_EQ(m.at(voc.index("25")), 2);
  EXPECT_EQ(m.at(voc.index("3")), 1);
  EXPECT_EQ(m.at(voc.index("34")), 1);
  EXPECT_EQ(m.at(voc.index("3.93")), 1);
  EXPECT_EQ(m.at(voc.index("400")), 1);
  EXPECT_EQ(m.at(voc.index("4th-qtr")), 1);
  EXPECT_EQ(m.at(voc.index("61 5/16")), 1);
  EXPECT_EQ(m.at(voc.index("65")), 1);
  EXPECT_EQ(m.at(voc.index("7:30")), 1);
  EXPECT_EQ(m.at(voc.index("75")), 1);
  EXPECT_EQ(m.at(voc.index("8.8")), 1);
  EXPECT_EQ(m.at(voc.index("a")), 11);
  EXPECT_EQ(m.at(voc.index("about")), 4);
  EXPECT_EQ(m.at(voc.index("addition")), 1);
  EXPECT_EQ(m.at(voc.index("adjust")), 1);
  EXPECT_EQ(m.at(voc.index("against")), 1);
  EXPECT_EQ(m.at(voc.index("alex")), 1);
  EXPECT_EQ(m.at(voc.index("already")), 1);
  EXPECT_EQ(m.at(voc.index("a.m.")), 1);
  EXPECT_EQ(m.at(voc.index("amid")), 2);
  EXPECT_EQ(m.at(voc.index("an")), 2);
  EXPECT_EQ(m.at(voc.index("analyst")), 2);
  EXPECT_EQ(m.at(voc.index("analysts")), 6);
  EXPECT_EQ(m.at(voc.index("and")), 4);
  EXPECT_EQ(m.at(voc.index("are")), 2);
  EXPECT_EQ(m.at(voc.index("as")), 10);
  EXPECT_EQ(m.at(voc.index("assault")), 1);
  EXPECT_EQ(m.at(voc.index("at")), 3);
  EXPECT_EQ(m.at(voc.index("attached")), 1);
  EXPECT_EQ(m.at(voc.index("average")), 1);
  EXPECT_EQ(m.at(voc.index("b.")), 1);
  EXPECT_EQ(m.at(voc.index("back")), 1);
  EXPECT_EQ(m.at(voc.index("barney")), 1);
  EXPECT_EQ(m.at(voc.index("based")), 1);
  EXPECT_EQ(m.at(voc.index("bc-kodak-forecast-update1-bloom")), 1);
  EXPECT_EQ(m.at(voc.index("be")), 4);
  EXPECT_EQ(m.at(voc.index("beat")), 1);
  EXPECT_EQ(m.at(voc.index("billion")), 2);
  EXPECT_EQ(m.at(voc.index("bloomberg")), 1);
  EXPECT_EQ(m.at(voc.index("both")), 1);
  EXPECT_EQ(m.at(voc.index("business")), 2);
  EXPECT_EQ(m.at(voc.index("by")), 3);
  EXPECT_EQ(m.at(voc.index("category")), 1);
  EXPECT_EQ(m.at(voc.index("cause")), 1);
  EXPECT_EQ(m.at(voc.index("cents")), 1);
  EXPECT_EQ(m.at(voc.index("change")), 1);
  EXPECT_EQ(m.at(voc.index("charge")), 1);
  EXPECT_EQ(m.at(voc.index("charles")), 1);
  EXPECT_EQ(m.at(voc.index("co.")), 2);
  EXPECT_EQ(m.at(voc.index("color")), 1);
  EXPECT_EQ(m.at(voc.index("comment")), 1);
  EXPECT_EQ(m.at(voc.index("compact")), 1);
  EXPECT_EQ(m.at(voc.index("companies")), 1);
  EXPECT_EQ(m.at(voc.index("company")), 1);
  EXPECT_EQ(m.at(voc.index("competition")), 1);
  EXPECT_EQ(m.at(voc.index("computers")), 1);
  EXPECT_EQ(m.at(voc.index("consumer")), 2);
  EXPECT_EQ(m.at(voc.index("continued")), 1);
  EXPECT_EQ(m.at(voc.index("currency")), 1);
  EXPECT_EQ(m.at(voc.index("cut")), 2);
  EXPECT_EQ(m.at(voc.index("decline")), 1);
  EXPECT_EQ(m.at(voc.index("declined")), 2);
  EXPECT_EQ(m.at(voc.index("degree")), 1);
  EXPECT_EQ(m.at(voc.index("deutsche")), 1);
  EXPECT_EQ(m.at(voc.index("digital")), 5);
  EXPECT_EQ(m.at(voc.index("diluted")), 1);
  EXPECT_EQ(m.at(voc.index("disclosed")), 1);
  EXPECT_EQ(m.at(voc.index("discount")), 2);
  EXPECT_EQ(m.at(voc.index("discs")), 1);
  EXPECT_EQ(m.at(voc.index("does")), 1);
  EXPECT_EQ(m.at(voc.index("dollar")), 4);
  EXPECT_EQ(m.at(voc.index("dollars")), 1);
  EXPECT_EQ(m.at(voc.index("dominates")), 2);
  EXPECT_EQ(m.at(voc.index("double-digit")), 1);
  EXPECT_EQ(m.at(voc.index("down")), 3);
  EXPECT_EQ(m.at(voc.index("earned")), 1);
  EXPECT_EQ(m.at(voc.index("earnings")), 3);
  EXPECT_EQ(m.at(voc.index("eastman")), 1);
  EXPECT_EQ(m.at(voc.index("effort")), 1);
  EXPECT_EQ(m.at(voc.index("established")), 1);
  EXPECT_EQ(m.at(voc.index("estimate")), 2);
  EXPECT_EQ(m.at(voc.index("evidence")), 1);
  EXPECT_EQ(m.at(voc.index("expect")), 1);
  EXPECT_EQ(m.at(voc.index("expected")), 2);
  EXPECT_EQ(m.at(voc.index("extent")), 1);
  EXPECT_EQ(m.at(voc.index("fall")), 1);
  EXPECT_EQ(m.at(voc.index("fell")), 2);
  EXPECT_EQ(m.at(voc.index("fewer")), 1);
  EXPECT_EQ(m.at(voc.index("film")), 11);
  EXPECT_EQ(m.at(voc.index("financial")), 1);
  EXPECT_EQ(m.at(voc.index("first")), 2);
  EXPECT_EQ(m.at(voc.index("flexibility")), 1);
  EXPECT_EQ(m.at(voc.index("focus")), 1);
  EXPECT_EQ(m.at(voc.index("for")), 5);
  EXPECT_EQ(m.at(voc.index("fourth-quarter")), 3);
  EXPECT_EQ(m.at(voc.index("from")), 5);
  EXPECT_EQ(m.at(voc.index("fuji")), 5);
  EXPECT_EQ(m.at(voc.index("gain")), 1);
  EXPECT_EQ(m.at(voc.index("getting")), 1);
  EXPECT_EQ(m.at(voc.index("going")), 1);
  EXPECT_EQ(m.at(voc.index("gold")), 1);
  EXPECT_EQ(m.at(voc.index("growing")), 1);
  EXPECT_EQ(m.at(voc.index("had")), 1);
  EXPECT_EQ(m.at(voc.index("has")), 3);
  EXPECT_EQ(m.at(voc.index("have")), 1);
  EXPECT_EQ(m.at(voc.index("he")), 1);
  EXPECT_EQ(m.at(voc.index("henderson")), 4);
  EXPECT_EQ(m.at(voc.index("high-profit")), 1);
  EXPECT_EQ(m.at(voc.index("his")), 1);
  EXPECT_EQ(m.at(voc.index("hit")), 1);
  EXPECT_EQ(m.at(voc.index("hold")), 2);
  EXPECT_EQ(m.at(voc.index("holds")), 1);
  EXPECT_EQ(m.at(voc.index("holiday")), 1);
  EXPECT_EQ(m.at(voc.index("holidays")), 1);
  EXPECT_EQ(m.at(voc.index("hurts")), 1);
  EXPECT_EQ(m.at(voc.index("ibes")), 1);
  EXPECT_EQ(m.at(voc.index("imaging")), 5);
  EXPECT_EQ(m.at(voc.index("in")), 13);
  EXPECT_EQ(m.at(voc.index("inc.")), 2);
  EXPECT_EQ(m.at(voc.index("include")), 1);
  EXPECT_EQ(m.at(voc.index("includes")), 1);
  EXPECT_EQ(m.at(voc.index("international")), 1);
  EXPECT_EQ(m.at(voc.index("is")), 3);
  EXPECT_EQ(m.at(voc.index("issued")), 1);
  EXPECT_EQ(m.at(voc.index("it")), 4);
  EXPECT_EQ(m.at(voc.index("its")), 5);
  EXPECT_EQ(m.at(voc.index("jan.")), 1);
  EXPECT_EQ(m.at(voc.index("japan")), 1);
  EXPECT_EQ(m.at(voc.index("jobs")), 1);
  EXPECT_EQ(m.at(voc.index("jonathan")), 1);
  EXPECT_EQ(m.at(voc.index("just")), 2);
  EXPECT_EQ(m.at(voc.index("key")), 1);
  EXPECT_EQ(m.at(voc.index("kodak")), 19);
  EXPECT_EQ(m.at(voc.index("largest")), 1);
  EXPECT_EQ(m.at(voc.index("last")), 5);
  EXPECT_EQ(m.at(voc.index("local")), 1);
  EXPECT_EQ(m.at(voc.index("looking")), 1);
  EXPECT_EQ(m.at(voc.index("loss")), 1);
  EXPECT_EQ(m.at(voc.index("losses")), 4);
  EXPECT_EQ(m.at(voc.index("lost")), 1);
  EXPECT_EQ(m.at(voc.index("lower")), 4);
  EXPECT_EQ(m.at(voc.index("-lrb-")), 3);
  EXPECT_EQ(m.at(voc.index("made")), 1);
  EXPECT_EQ(m.at(voc.index("margins")), 1);
  EXPECT_EQ(m.at(voc.index("mark")), 1);
  EXPECT_EQ(m.at(voc.index("market")), 4);
  EXPECT_EQ(m.at(voc.index("market-share")), 1);
  EXPECT_EQ(m.at(voc.index("microfilm")), 1);
  EXPECT_EQ(m.at(voc.index("million")), 1);
  EXPECT_EQ(m.at(voc.index("month")), 1);
  EXPECT_EQ(m.at(voc.index("more")), 2);
  EXPECT_EQ(m.at(voc.index("mounting")), 1);
  EXPECT_EQ(m.at(voc.index("much")), 2);
  EXPECT_EQ(m.at(voc.index("need")), 1);
  EXPECT_EQ(m.at(voc.index("negative")), 1);
  EXPECT_EQ(m.at(voc.index("neutral")), 1);
  EXPECT_EQ(m.at(voc.index("new")), 3);
  EXPECT_EQ(m.at(voc.index("n't")), 1);
  EXPECT_EQ(m.at(voc.index("of")), 10);
  EXPECT_EQ(m.at(voc.index("older")), 1);
  EXPECT_EQ(m.at(voc.index("on")), 8);
  EXPECT_EQ(m.at(voc.index("operating")), 1);
  EXPECT_EQ(m.at(voc.index("or")), 3);
  EXPECT_EQ(m.at(voc.index("others")), 1);
  EXPECT_EQ(m.at(voc.index("overseas")), 1);
  EXPECT_EQ(m.at(voc.index("percent")), 10);
  EXPECT_EQ(m.at(voc.index("percentage")), 1);
  EXPECT_EQ(m.at(voc.index("photo")), 2);
  EXPECT_EQ(m.at(voc.index("photography")), 1);
  EXPECT_EQ(m.at(voc.index("plans")), 1);
  EXPECT_EQ(m.at(voc.index("point")), 1);
  EXPECT_EQ(m.at(voc.index("pressure")), 1);
  EXPECT_EQ(m.at(voc.index("price")), 1);
  EXPECT_EQ(m.at(voc.index("prices")), 5);
  EXPECT_EQ(m.at(voc.index("products")), 1);
  EXPECT_EQ(m.at(voc.index("profit")), 1);
  EXPECT_EQ(m.at(voc.index("promotion")), 1);
  EXPECT_EQ(m.at(voc.index("prudential")), 1);
  EXPECT_EQ(m.at(voc.index("quarter")), 3);
  EXPECT_EQ(m.at(voc.index("quite")), 1);
  EXPECT_EQ(m.at(voc.index("rates")), 1);
  EXPECT_EQ(m.at(voc.index("rating")), 2);
  EXPECT_EQ(m.at(voc.index("'re")), 1);
  EXPECT_EQ(m.at(voc.index("release")), 1);
  EXPECT_EQ(m.at(voc.index("replace")), 1);
  EXPECT_EQ(m.at(voc.index("report")), 3);
  EXPECT_EQ(m.at(voc.index("respond")), 1);
  EXPECT_EQ(m.at(voc.index("result")), 1);
  EXPECT_EQ(m.at(voc.index("rises")), 1);
  EXPECT_EQ(m.at(voc.index("rising")), 3);
  EXPECT_EQ(m.at(voc.index("rochester")), 2);
  EXPECT_EQ(m.at(voc.index("rose")), 2);
  EXPECT_EQ(m.at(voc.index("rosenzweig")), 1);
  EXPECT_EQ(m.at(voc.index("-rrb-")), 3);
  EXPECT_EQ(m.at(voc.index("'s")), 7);
  EXPECT_EQ(m.at(voc.index("said")), 10);
  EXPECT_EQ(m.at(voc.index("sales")), 10);
  EXPECT_EQ(m.at(voc.index("salomon")), 1);
  EXPECT_EQ(m.at(voc.index("season")), 1);
  EXPECT_EQ(m.at(voc.index("second")), 1);
  EXPECT_EQ(m.at(voc.index("securities")), 1);
  EXPECT_EQ(m.at(voc.index("see")), 2);
  EXPECT_EQ(m.at(voc.index("seen")), 1);
  EXPECT_EQ(m.at(voc.index("september")), 1);
  EXPECT_EQ(m.at(voc.index("share")), 3);
  EXPECT_EQ(m.at(voc.index("shares")), 1);
  EXPECT_EQ(m.at(voc.index("shopping")), 1);
  EXPECT_EQ(m.at(voc.index("since")), 1);
  EXPECT_EQ(m.at(voc.index("slowing")), 1);
  EXPECT_EQ(m.at(voc.index("slumping")), 1);
  EXPECT_EQ(m.at(voc.index("smith")), 2);
  EXPECT_EQ(m.at(voc.index("socked")), 1);
  EXPECT_EQ(m.at(voc.index("sold")), 1);
  EXPECT_EQ(m.at(voc.index("some")), 1);
  EXPECT_EQ(m.at(voc.index("spends")), 1);
  EXPECT_EQ(m.at(voc.index("spokesman")), 1);
  EXPECT_EQ(m.at(voc.index("steal")), 1);
  EXPECT_EQ(m.at(voc.index("stem")), 1);
  EXPECT_EQ(m.at(voc.index("stock")), 3);
  EXPECT_EQ(m.at(voc.index("storage")), 1);
  EXPECT_EQ(m.at(voc.index("such")), 2);
  EXPECT_EQ(m.at(voc.index("surveyed")), 1);
  EXPECT_EQ(m.at(voc.index("systems")), 1);
  EXPECT_EQ(m.at(voc.index("take")), 1);
  EXPECT_EQ(m.at(voc.index("than")), 3);
  EXPECT_EQ(m.at(voc.index("that")), 6);
  EXPECT_EQ(m.at(voc.index("the")), 25);
  EXPECT_EQ(m.at(voc.index("they")), 2);
  EXPECT_EQ(m.at(voc.index("thing")), 1);
  EXPECT_EQ(m.at(voc.index("third")), 1);
  EXPECT_EQ(m.at(voc.index("through")), 1);
  EXPECT_EQ(m.at(voc.index("thursday")), 2);
  EXPECT_EQ(m.at(voc.index("time")), 2);
  EXPECT_EQ(m.at(voc.index("to")), 21);
  EXPECT_EQ(m.at(voc.index("today")), 1);
  EXPECT_EQ(m.at(voc.index("triple")), 1);
  EXPECT_EQ(m.at(voc.index("tumbling")), 1);
  EXPECT_EQ(m.at(voc.index("two")), 1);
  EXPECT_EQ(m.at(voc.index("unchanged")), 1);
  EXPECT_EQ(m.at(voc.index("update1")), 1);
  EXPECT_EQ(m.at(voc.index("u.s.")), 5);
  EXPECT_EQ(m.at(voc.index("u.s.-based")), 1);
  EXPECT_EQ(m.at(voc.index("using")), 1);
  EXPECT_EQ(m.at(voc.index("version")), 1);
  EXPECT_EQ(m.at(voc.index("we")), 1);
  EXPECT_EQ(m.at(voc.index("weak")), 1);
  EXPECT_EQ(m.at(voc.index("were")), 2);
  EXPECT_EQ(m.at(voc.index("what")), 1);
  EXPECT_EQ(m.at(voc.index("which")), 5);
  EXPECT_EQ(m.at(voc.index("while")), 1);
  EXPECT_EQ(m.at(voc.index("who")), 2);
  EXPECT_EQ(m.at(voc.index("widening")), 1);
  EXPECT_EQ(m.at(voc.index("wider")), 1);
  EXPECT_EQ(m.at(voc.index("will")), 2);
  EXPECT_EQ(m.at(voc.index("workforce")), 1);
  EXPECT_EQ(m.at(voc.index("world")), 1);
  EXPECT_EQ(m.at(voc.index("would")), 1);
  EXPECT_EQ(m.at(voc.index("wrote")), 1);
  EXPECT_EQ(m.at(voc.index("year")), 4);
  EXPECT_EQ(m.at(voc.index("year-earlier")), 1);
  EXPECT_EQ(m.at(voc.index("york")), 2);
  EXPECT_EQ(m.at(voc.index("york-based")), 1);
}
