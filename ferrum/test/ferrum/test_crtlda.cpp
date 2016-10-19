#include "gtest/gtest.h"

#include "concrete_util/io.h"
#include "ferrum/crtlda.hpp"
#include "ferrum/logging.hpp"

#include <iostream>

ferrum::AnnotatedToken<std::string> get_verb_token() {
  ferrum::AnnotatedToken<std::string> token;
  token.lemma("run");
  token.original("ran");
  token.pos("VBD");
  token.view("run-1.0");
  return token;
}
ferrum::AnnotatedToken<std::string> get_noun_token() {
  ferrum::AnnotatedToken<std::string> token;
  token.lemma("boy");
  token.original("boy");
  token.pos("NN");
  token.view("boy-1.0");
  return token;
}

ferrum::DocumentGRC<std::string,std::string,std::string> get_str_doc() {
  typedef std::string string;
  ferrum::DocumentGRC<std::string,std::string,std::string> doc("my pretend document");
  for(int ei = 0; ei < 3; ei++) {
    ferrum::Entity<string,string,string> entity;
    entity.canonical_name("Bob");
    for(int mi = 0; mi < 5; mi++) {
      ferrum::Mention<string, string> mention;
      mention.rel("nsubj-run");
      mention.gov(get_verb_token());
      mention.head(get_noun_token());  
      entity.add_mention(mention);
    }
    doc.add_entity(entity);
  }
  return doc;
}

TEST(Vocabulary, opLRBRRB){
  ferrum::Vocabulary<char> vocab;
  ASSERT_EQ(0, vocab('a'));
  ASSERT_EQ(1, vocab('b'));
  ASSERT_EQ(0, vocab('a'));
  ASSERT_EQ(2, vocab.num_words());
}

TEST(AnnotatedToken, create_str) {
  ferrum::AnnotatedToken<std::string> token = get_verb_token();
  ASSERT_EQ(token.original(), "ran");
  ASSERT_EQ(token.lemma(), "run");
  ASSERT_EQ(token.pos(), "VBD");
  ASSERT_EQ(token.view(), "run-1.0");
}

TEST(Mention, create_str) {
  typedef std::string string;
  ferrum::Mention<string, string> mention;
  mention.rel("nsubj-run");
  ferrum::AnnotatedToken<std::string> gov = get_verb_token();
  mention.gov(gov);
  ferrum::AnnotatedToken<std::string> head = get_noun_token();
  mention.head(head);
}

TEST(Entity, create_str) {
  typedef std::string string;
  ferrum::Entity<string,string,string> entity;
  entity.canonical_name("Bob");
  ferrum::Mention<string, string> mention;
  mention.rel("nsubj-run");
  mention.gov(get_verb_token());
  mention.head(get_noun_token());  
  ASSERT_EQ(entity.num_mentions(), 0);
  entity.add_mention(mention);
  ASSERT_EQ(entity.num_mentions(), 1);
  ferrum::Mention<string,string> ret_ment = entity[0];
  ASSERT_EQ(ret_ment.gov().view(), "run-1.0");
  ASSERT_EQ(ret_ment.gov().view(), mention.gov().view());
  ASSERT_NE(&ret_ment, &mention);
  ASSERT_NE(&(ret_ment.gov()), &(mention.gov()));
  entity.id("My_id_1");
  ASSERT_EQ("My_id_1", entity.id());
}

TEST(Document, create) {
  typedef ferrum::DocumentGRC<int, int, int> Doc;
  Doc my_doc("doc id");
  ASSERT_EQ(my_doc.id, "doc id");
}

TEST(Document, from_fake_comm) {
  typedef std::string string;
  concrete::Communication comm;
  comm.__set_id("My communication id");
  const ferrum::VerbGovernedPruner< string, string > vgp(comm);
  ferrum::DocumentGRC<string, string, string> doc(comm, vgp);
  ASSERT_EQ("My communication id", comm.id);
  ASSERT_EQ("My communication id", doc.id);
  ASSERT_EQ(doc.id,comm.id);
}

TEST(VerbGovernedPruner, from_compact_AFP_ENG_19940531_0390_comm) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("AFP_ENG_19940531.0390", comm.id);

  typedef std::string string;
  const ferrum::VerbGovernedPruner< string, string > vgp(comm);
  concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention =
    concrete_util::mention_id_to_mention(comm, "Stanford");
  concrete::UUID mention_id = comm.entityMentionSetList[0].mentionList[0].uuid;
  ASSERT_NE(mention_id_to_mention.find(mention_id), mention_id_to_mention.end());
  concrete::EntityMention em = mention_id_to_mention[mention_id];
  const concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization =
    concrete_util::mention_id_to_tokenization(comm, "Stanford");
  ASSERT_EQ(1, mention_id_to_tokenization.count(mention_id));
  ASSERT_NE(mention_id_to_tokenization.find(mention_id), mention_id_to_tokenization.end());
  const concrete::Tokenization tokenization = mention_id_to_tokenization.at(mention_id);
  DEBUG << "Mention id: " << mention_id.uuidString;
  ASSERT_EQ("6f0253b9-b9a5-4982-9be0-39c3ae79d663", em.uuid.uuidString);
  ASSERT_EQ("e5f8ba4a-a22f-4e7e-8499-f05c4c0065fa", em.tokens.tokenizationId.uuidString);
  DEBUG << "Trying to inspect tokenization " << tokenization.uuid.uuidString ;
  DEBUG << "Tokenization is at " << &tokenization;
  ASSERT_TRUE(tokenization.__isset.dependencyParseList) << "Tokenization " << tokenization.uuid.uuidString << " does not have a dependency parse list set";
  const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(tokenization, "col-ccproc-deps");
  ASSERT_TRUE(dep_parse != NULL);
  vgp.prune(em);
}

TEST(SituationGovernedPruner, from_compact_NYT_ENG_19980113_0597_comm) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/NYT_ENG_19980113.0597.tcompact.with_semafor.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("NYT_ENG_19980113.0597", comm.id);

  typedef std::string string;
  const ferrum::SituationGovernedPruner< string, string > sgp(comm);
  concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention =
    concrete_util::mention_id_to_mention(comm, "Stanford");
  concrete::UUID mention_id = comm.entityMentionSetList[0].mentionList[0].uuid;
  ASSERT_NE(mention_id_to_mention.find(mention_id), mention_id_to_mention.end());
  concrete::EntityMention em = mention_id_to_mention[mention_id];
  const concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization =
    concrete_util::mention_id_to_tokenization(comm, "Stanford");
  ASSERT_EQ(1, mention_id_to_tokenization.count(mention_id));
  ASSERT_NE(mention_id_to_tokenization.find(mention_id), mention_id_to_tokenization.end());
  const concrete::Tokenization tokenization = mention_id_to_tokenization.at(mention_id);
  DEBUG << "Mention id: " << mention_id.uuidString;
  ASSERT_EQ("beae7914-9205-420c-b762-7a4bfa8d388c", em.uuid.uuidString);
  ASSERT_EQ("2fc3e7e9-e103-4791-a82e-1adb525f647d", em.tokens.tokenizationId.uuidString);
  DEBUG << "Trying to inspect tokenization " << tokenization.uuid.uuidString ;
  DEBUG << "Tokenization is at " << &tokenization;
  ASSERT_TRUE(tokenization.__isset.dependencyParseList) << "Tokenization " << tokenization.uuid.uuidString << " does not have a dependency parse list set";
  const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(tokenization, "col-ccproc-deps");
  ASSERT_TRUE(dep_parse != NULL);
  std::vector< ferrum::Mention<string, string> > created_mentions = sgp.prune(em);
}

TEST(Document, from_compact_AFP_EN_19940531_0390_comm) {
  typedef std::string string;
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  DEBUG << "preparing to read " << name;
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  DEBUG << "successfully loaded communication from " << name;
  ASSERT_EQ("AFP_ENG_19940531.0390", comm.id);
  const ferrum::VerbGovernedPruner< string, string > vgp(comm);
  DEBUG << "Made VerbGoverendPruner vgp for comm " << comm.id;
  ferrum::DocumentGRC<string, string, string> doc(comm, vgp);
  //ASSERT_EQ("AFP_ENG_19940531.0390", doc.id);
}

TEST(SituationGovernedPruner, logic_use_lex_true) {
  typedef std::string string;
  concrete::Communication comm;
  ferrum::SituationGovernedPruner< string, string > sgp(comm);
  sgp.use_lexical(true);
  ferrum::AnnotatedToken<string> verb = get_verb_token();
  string gov = sgp.make_gov_view(verb, verb.view());
  ASSERT_EQ("run", gov);
  string rel = sgp.make_relation("FOO", "foo", gov, verb.lemma());
  ASSERT_EQ("foo-run", rel);
}
TEST(SituationGovernedPruner, logic_use_lex_false) {
  typedef std::string string;
  concrete::Communication comm;
  ferrum::SituationGovernedPruner< string, string > sgp(comm);
  sgp.use_lexical(false);
  ferrum::AnnotatedToken<string> verb = get_verb_token();
  string gov = sgp.make_gov_view(verb, verb.view());
  ASSERT_EQ("run-1.0", gov);
  string rel = sgp.make_relation("FOO", "foo", gov, verb.lemma());
  ASSERT_EQ("FOO-run-1.0", rel);
}
TEST(SituationGovernedPrunerMask, logic_use_lex_true) {
  typedef std::string string;
  concrete::Communication comm;
  ferrum::WordMapper wm("OOV");
  wm.add_word("ran","0101");
  ferrum::SituationGovernedPrunerMaskLemma< string, string > sgp(comm, wm);
  sgp.use_lexical(true);
  ferrum::AnnotatedToken<string> verb = get_verb_token();
  INFO << verb.view();
  string gov = sgp.make_gov_view(verb, verb.view());
  ASSERT_EQ("0101", gov);
  string rel = sgp.make_relation("FOO", "foo", gov, verb.lemma());
  ASSERT_EQ("foo-0101", rel);
}
TEST(SituationGovernedPrunerMask, logic_use_lex_false) {
  typedef std::string string;
  concrete::Communication comm;
  ferrum::WordMapper wm("OOV");
  wm.add_word("ran","0101");
  ferrum::SituationGovernedPrunerMaskLemma< string, string > sgp(comm, wm);
  sgp.use_lexical(false);
  ferrum::AnnotatedToken<string> verb = get_verb_token();
  string gov = sgp.make_gov_view(verb, verb.view());
  ASSERT_EQ("0101", gov);
  string rel = sgp.make_relation("FOO", "foo", gov, verb.lemma());
  ASSERT_EQ("foo-0101", rel);
}

TEST(Document, from_compact_NYT_ENG_19980113_0597_comm) {
  typedef std::string string;
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/NYT_ENG_19980113.0597.tcompact.with_semafor.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("NYT_ENG_19980113.0597", comm.id);
  const ferrum::SituationGovernedPruner< string, string > sgp(comm);
  ferrum::DocumentGRC<string, string, string> doc(comm, sgp);
  ASSERT_GT(doc.num_entities(), 0);
  INFO << "Created ferrum::DocumentGRC has " << doc.num_entities() << " entities";
}

TEST(Document, fill) {
  ferrum::DocumentGRC<std::string,std::string,std::string> doc = get_str_doc();
  ASSERT_EQ(doc.num_entities(), 3);
  for(int i = 0; i < 3; i++) {
    ferrum::Entity<std::string, std::string, std::string> entity = doc[i];
    ASSERT_EQ(entity.num_mentions(), 5);
    ASSERT_EQ(entity.canonical_name(), "Bob");
    for(int j = 0; j < 3; j++) {
      ferrum::Mention<std::string, std::string> mention = entity[j];
      ASSERT_EQ(mention.gov().view(), entity[j].gov().view());
      ASSERT_EQ(mention.gov().view(), doc[i][j].gov().view());
      ASSERT_EQ(entity[j].gov().view(), doc[i][j].gov().view());
    }
  }
}

TEST(Document, iostream) {
  ferrum::DocumentGRC<std::string, std::string, std::string> doc("my id");
  std::cout << "Testing doc stream: " << doc << std::endl;
}

TEST(InMemoryCorpus, create) {
  typedef ferrum::DocumentGRC<std::string,std::string,std::string> Doc;
  ferrum::InMemoryCorpus<Doc> corpus("my corpus");
}

TEST(InMemoryCorpus, fill) {
  typedef ferrum::DocumentGRC<std::string,std::string,std::string> Doc;
  ferrum::InMemoryCorpus<Doc> corpus("my corpus");
  for(int di = 0; di < 10; di++) {
    Doc doc = get_str_doc();
    corpus.add_document(doc);
  }
  ASSERT_EQ(corpus.num_docs(), 10);
}

void fillSHP(ferrum::SymmetricHyperparams *shp) {
  shp->h_gov = .1;
  shp->h_rel = .1;
  shp->h_theta = .1;
}

TEST(SymmetricHyperparams, valueSet) {
  ferrum::SymmetricHyperparams shp = ferrum::SymmetricHyperparams();
  fillSHP(&shp);
  ASSERT_EQ(shp.h_gov, .1);
}

TEST(DiscreteModel, create) {
  ferrum::SymmetricHyperparams shp = ferrum::SymmetricHyperparams();
  fillSHP(&shp);
  ferrum::DiscreteModel dm(10, 13, &shp);
}

TEST(InMemoryCorpus, iterate) {
  typedef ferrum::DocumentGRC<std::string,std::string,std::string> Doc;
  ferrum::InMemoryCorpus<Doc> corpus("my corpus");
  for(int di = 0; di < 10; di++) {
    Doc doc = get_str_doc();
    corpus.add_document(doc);
  }
  ASSERT_EQ(corpus.num_docs(), 10);
  int num_from_iterator = 0;
  for(ferrum::InMemoryCorpus<Doc>::const_iterator it = corpus.begin();
      it != corpus.end();
      ++it) {
    ASSERT_EQ(num_from_iterator, it->iteration_idx);
    const Doc& doc = *(it->document);
    INFO << (it->iteration_idx) << " :: " << doc.id;
    ++num_from_iterator;
  }
}

TEST(InMemoryCorpus, random_subset) {
  typedef ferrum::DocumentGRC<std::string,std::string,std::string> Doc;
  ferrum::InMemoryCorpus<Doc> corpus("my corpus");
  for(int di = 0; di < 10; di++) {
    Doc doc = get_str_doc();
    corpus.add_document(doc);
  }
  ASSERT_EQ(corpus.num_docs(), 10);
  int num_from_iterator;
  num_from_iterator = 0;
  ferrum::InMemoryCorpus<Doc> rand_corp = *(corpus.random_subset(3));
  for(ferrum::InMemoryCorpus<Doc>::const_iterator it = rand_corp.begin();
      it != rand_corp.end();
      ++it) {
    const Doc& doc = *(it->document);
    INFO << (it->iteration_idx) << " :: " << doc.id;
    ++num_from_iterator;
  }
  ASSERT_EQ(3, num_from_iterator);
}

TEST(Vocabulary, update_with) {
  typedef std::string string;
  ferrum::Vocabulary<string> vocab0("my_oov");
  vocab0.make_word("foo");
  vocab0.make_word("bar");
  ferrum::Vocabulary<string> vocab1("my_oov2");
  vocab1.make_word("bar");
  vocab1.make_word("foo");
  vocab1.make_word("baz");
  ASSERT_EQ(3, vocab0.num_words());
  std::vector<int> mapper = vocab0.update_with(vocab1);
  ASSERT_EQ(4, mapper.size());
  ASSERT_EQ(0, mapper[0]); // OOV --> OOV
  ASSERT_EQ(2, mapper[1]); // bar --> foo
  ASSERT_EQ(1, mapper[2]); // foo --> bar
  ASSERT_EQ(3, mapper[3]); // baz --> baz (newly created)
  ASSERT_EQ(4, vocab0.num_words());
}
