#include "gtest/gtest.h"

//#include "concrete/communication_types.h"
#include "ferrum/concrete.hpp"
#include "ferrum/mathops.hpp"

#include "concrete_util/uuid_util.h"
#include "concrete_util/io.h"

#include "ferrum/logging.hpp"
#include <cstdio>
#include <unordered_map>

TEST(UUIDString, create) {
  concrete::util::uuid_factory uuid_maker;
  DEBUG << "Generated UUID is " << uuid_maker.get_uuid();
}

TEST(Version, correct_version) {
  ASSERT_EQ("4.8", concrete::CONCRETE_VERSION);
}

TEST(Communication, readBinaryThroughGen) {
  concrete::Communication communication;
  concrete::util::concrete_io concrete_reader;
  concrete_reader.deserialize<concrete::util::TBinaryProtocol, concrete::Communication>(&communication, "resources/AFP_ENG_19940531.0390.tbinary.concrete");
  ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
}

TEST(Communication, readCompactThroughGen) {
  concrete::Communication communication;
  concrete::util::concrete_io concrete_reader;
  concrete_reader.deserialize<concrete::util::TCompactProtocol, concrete::Communication>(&communication, "resources/AFP_ENG_19940531.0390.tcompact.concrete");
  ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
}

TEST(Communication, readTBinaryProtocol) {
  concrete::Communication communication;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tbinary.concrete";
  concrete_reader.deserialize_binary<concrete::Communication>(&communication, name);
  ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
}

TEST(Communication, readTCompactProtocol) {
  concrete::Communication communication;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&communication, name);
  ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
}

TEST(EntitySet, find) {
  concrete::Communication communication;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&communication, name);
  const concrete::EntitySet* es = concrete_util::first_entity_set(communication,
								  "Stanford");
  ASSERT_TRUE(es != NULL);
}

TEST(ConcreteUtil, first_set_with_name_mention_from_compact_AFP_EN_19940531_0390) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  DEBUG << "preparing to read " << name;
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  DEBUG << "successfully loaded communication from " << name;
  ASSERT_EQ("AFP_ENG_19940531.0390", comm.id);
  const concrete::EntityMentionSet* ems =
    concrete_util::first_set_with_name<concrete::EntityMentionSet>(comm.entityMentionSetList,
								   "Stanford");
  DEBUG << "Creating mention id to tokenization mapping";
  ASSERT_TRUE(ems != NULL);
  ASSERT_EQ(50, ems->mentionList.size());
}

TEST(ConcreteUtil, iterate_sentences) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  DEBUG << "preparing to read " << name;
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  DEBUG << "successfully loaded communication from " << name;
  ASSERT_EQ("AFP_ENG_19940531.0390", comm.id);
  int num_toks = 0;
  for(concrete::Section section : comm.sectionList) {
    if(!section.__isset.sentenceList) continue;
    for(concrete::Sentence sentence : section.sentenceList) {
      ASSERT_TRUE(sentence.__isset.tokenization);
      ++num_toks;
    }
  }
  ASSERT_EQ(9, num_toks);
}

TEST(ConcreteUtil, uuid_hash) {
  concrete::util::uuid_factory uf;
  concrete::UUID uuid;
  uuid.__set_uuidString(uf.get_uuid());
  const size_t uuid_hash = concrete_util::uuid_hash()(uuid);
  INFO << "Created hash of UUID[" << uuid.uuidString << "] = " << uuid_hash;
}

TEST(ConcreteUtil, uuid_map_no_alias) {
  std::unordered_map< const concrete::UUID, int, concrete_util::uuid_hash > tutt;
  concrete::util::uuid_factory uf;
  concrete::UUID uuid;
  uuid.__set_uuidString(uf.get_uuid());
  DEBUG << "Attempting to add element with value " << uuid.uuidString << " to the map";
  tutt[uuid] = 101;
  ASSERT_EQ(1, tutt.size());
  ASSERT_EQ(101, tutt[uuid]);
}


TEST(ConcreteUtil, uuid_map) {
  concrete_util::uuid_map<int> tutt;
  concrete::util::uuid_factory uf;
  concrete::UUID uuid;
  uuid.__set_uuidString(uf.get_uuid());
  tutt[uuid] = 101;
  ASSERT_EQ(1, tutt.size());
  ASSERT_EQ(101, tutt[uuid]);
}


TEST(ConcreteUtil, tok_id_to_tok_map_from_compact_AFP_EN_19940531_0390_comm) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("AFP_ENG_19940531.0390", comm.id);
  const concrete_util::uuid_map<concrete::Tokenization> tok_id_to_tptr =
    concrete_util::tokenization_id_to_tokenization(comm);
  // check to make sure all UUID strings are the same
  for(auto& item : tok_id_to_tptr) {
    ASSERT_EQ(item.first.uuidString, item.second.uuid.uuidString);
  };
  ASSERT_EQ(9, tok_id_to_tptr.size());
}


TEST(ConcreteUtil, mention_id_to_tokenization_map_from_compact_AFP_EN_19940531_0390_comm) {
  concrete::Communication comm;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&comm, name);
  ASSERT_EQ("AFP_ENG_19940531.0390", comm.id);
  const concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization =
    concrete_util::mention_id_to_tokenization(comm, "Stanford");
  ASSERT_EQ(50, mention_id_to_tokenization.size());
}

TEST(ConcreteSmartWriter, write_fake_uuid_dev_null) {
  concrete::util::ConcreteSmartWriter<concrete::util::TBinaryProtocol> csw("/dev/null");
  auto writer = csw.get();
  concrete::UUID uuid;
  uuid.__set_uuidString("my id");
  uuid.write(writer);
}

TEST(ConcreteSmartWriter, write_fake_uuid_tmp_file) {
  char buffer [L_tmpnam];
  IGNORE_RESULT(static_cast<void>(std::tmpnam(buffer)));
  std::string sbuf(buffer);
  concrete::util::ConcreteSmartWriter<concrete::util::TBinaryProtocol> csw(sbuf);
  auto writer = csw.get();
  std::string fname = csw.name();
  concrete::UUID uuid;
  uuid.__set_uuidString("my id");
  uuid.write(writer);
  csw.close_csw();
  // check the file exists
  // open
  concrete::util::concrete_io cread;
  concrete::UUID ruuid;
  cread.deserialize<concrete::util::TBinaryProtocol, concrete::UUID>
    (
     &ruuid,
     fname.c_str()
     );
  ASSERT_EQ(uuid.uuidString, ruuid.uuidString);
  ASSERT_EQ("my id", ruuid.uuidString);
  ASSERT_EQ(0, remove(fname.c_str()));
}

TEST(StructIterator, Tokenization_4x2) {
  concrete::Communication comm;
  const int num_sects = 4;
  std::vector<int> num_sents(num_sects);
  int num_toks = 0;
  for(int isect = 0; isect < num_sects; ++isect) {
    int n_sent = 1;
    num_sents[isect] = n_sent;
    concrete::Section section;
    for(int isent = 0; isent < n_sent; ++isent) {
      concrete::Sentence sent;
      concrete::Tokenization tok;
      sent.__set_tokenization(tok);
      num_toks++;
      section.sentenceList.push_back(sent);
      section.__isset.sentenceList = true;
    }
    comm.sectionList.push_back(section);
    comm.__isset.sectionList = true;
  }
  int nt = 0;
  concrete::util::StructIterator<concrete::Tokenization, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator> titerstr(comm);
  auto titer = titerstr.begin();
  auto tend = titerstr.end();
  for(; titer != tend; ++titer) {
    ++nt;
  }
  ASSERT_EQ(num_toks, nt);
}

TEST(TokenizationIterator, Tokenization_4x2) {
  concrete::Communication comm;
  const int num_sects = 4;
  std::vector<int> num_sents(num_sects);
  int num_toks = 0;
  for(int isect = 0; isect < num_sects; ++isect) {
    int n_sent = 1;
    num_sents[isect] = n_sent;
    concrete::Section section;
    for(int isent = 0; isent < n_sent; ++isent) {
      concrete::Sentence sent;
      concrete::Tokenization tok;
      sent.__set_tokenization(tok);
      num_toks++;
      section.sentenceList.push_back(sent);
      section.__isset.sentenceList = true;
    }
    comm.sectionList.push_back(section);
    comm.__isset.sectionList = true;
  }
  int nt = 0;
  concrete::util::TokenizationIterator titerstr(comm);
  auto titer = titerstr.begin();
  auto tend = titerstr.end();
  for(; titer != tend; ++titer) {
    ++nt;
  }
  ASSERT_EQ(num_toks, nt);
}

TEST(StructIteratorTokenization, Tokenization_4x2) {
  concrete::Communication comm;
  const int num_sects = 4;
  std::vector<int> num_sents(num_sects);
  int num_toks = 0;
  for(int isect = 0; isect < num_sects; ++isect) {
    int n_sent = 1;
    num_sents[isect] = n_sent;
    concrete::Section section;
    for(int isent = 0; isent < n_sent; ++isent) {
      concrete::Sentence sent;
      concrete::Tokenization tok;
      sent.__set_tokenization(tok);
      num_toks++;
      section.sentenceList.push_back(sent);
      section.__isset.sentenceList = true;
    }
    comm.sectionList.push_back(section);
    comm.__isset.sectionList = true;
  }
  int nt = 0;
  concrete::util::StructIterator<concrete::Tokenization> titerstr(comm);
  auto titer = titerstr.begin();
  auto tend = titerstr.end();
  for(; titer != tend; ++titer) {
    ++nt;
  }
  ASSERT_EQ(num_toks, nt);
}

TEST(StructIterator, Tokenization) {
  for(int i = 0; i < 10; ++i) {
    concrete::Communication comm;
    const int num_sects = (int)gsl_rng_uniform_int(mathops::rnd_gen, 5);;
    std::vector<int> num_sents(num_sects);
    int num_toks = 0;
    std::vector<std::string> uuids;
    for(int isect = 0; isect < num_sects; ++isect) {
      int n_sent = (int)gsl_rng_uniform_int(mathops::rnd_gen, 7);
      num_sents[isect] = n_sent;
      concrete::Section section;
      for(int isent = 0; isent < n_sent; ++isent) {
	concrete::Sentence sent;
	concrete::Tokenization tok;
	concrete::UUID uuid;
	uuid.__set_uuidString("uuid/" + std::to_string(num_toks));
	uuids.push_back(uuid.uuidString);
	tok.__set_uuid(uuid);
	sent.__set_tokenization(tok);
	num_toks++;
	section.sentenceList.push_back(sent);
	section.__isset.sentenceList = true;
      }
      comm.sectionList.push_back(section);
      comm.__isset.sectionList = true;
    }
    {
      int nt = 0;
      concrete::util::StructIterator<concrete::Tokenization, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator> titerstr(comm);
      auto titer = titerstr.begin();
      auto tend = titerstr.end();
      for(; titer != tend; ++titer) {
	EXPECT_EQ(uuids.at(nt), titer->uuid.uuidString);
	++nt;
	if(nt > num_toks) break;
      }
      ASSERT_EQ(num_toks, nt);
    }
    {
      int nt = 0;
      for(const concrete::Tokenization& tokenization : concrete::util::StructIterator<concrete::Tokenization, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator>(comm)) {
	(void)tokenization;
	++nt;
      }
      ASSERT_EQ(num_toks, nt);
    }
  }
}

TEST(StructIterator, no_Tokenizations) {
  concrete::Communication comm;
  const int num_sects = 4;
  std::vector<int> num_sents(num_sects);
  for(int isect = 0; isect < num_sects; ++isect) {
    int n_sent = isect % 2 ? 1 : 2;
    num_sents[isect] = n_sent;
    concrete::Section section;
    for(int isent = 0; isent < n_sent; ++isent) {
      concrete::Sentence sent;
      section.sentenceList.push_back(sent);
      section.__isset.sentenceList = true;
    }
    comm.sectionList.push_back(section);
    comm.__isset.sectionList = true;
  }

  int nt = 0;
  concrete::util::StructIterator<concrete::Tokenization, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator> titerstr(comm);
  auto titer = titerstr.begin();
  auto tend = titerstr.end();
  for(; titer != tend; ++titer) {
    ++nt;
  }
  ASSERT_EQ(0, nt);
}

TEST(StructIterator, Communication) {
  concrete::Communication comm;
  concrete::util::StructIterator<concrete::Communication, int> titerstr(comm);
  int nc = 0;
  for(const concrete::Communication& c : titerstr) {
    (void)c;
    ++nc;
  }
  ASSERT_EQ(1, nc);
}

TEST(StructIterator, Token) {
  for(int i = 0; i < 10; ++i) {
    concrete::Communication comm;
    const int num_sects = (int)gsl_rng_uniform_int(mathops::rnd_gen, 5);
    std::vector<int> num_sents(num_sects);
    std::vector< std::vector<int> > num_tokens(num_sects);
    int num_toks = 0;
    std::vector<std::string> expected_text;
    for(int isect = 0; isect < num_sects; ++isect) {
      int n_sent = (int)gsl_rng_uniform_int(mathops::rnd_gen, 7);
      num_sents[isect] = n_sent;
      num_tokens[isect] = std::vector<int>(n_sent);
      concrete::Section section;
      for(int isent = 0; isent < n_sent; ++isent) {
	concrete::Sentence sent;
	concrete::Tokenization tok;
	int n_token = (int)gsl_rng_uniform_int(mathops::rnd_gen, 14);
	num_tokens[isect][isent] = n_token;
	concrete::TokenList tl;
	for(int itok = 0; itok < n_token; ++itok) {
	  concrete::Token token;
	  token.__set_text("token" + std::to_string(isect) + "/" + std::to_string(isent) +"/" + std::to_string(itok));
	  expected_text.push_back(token.text);
	  //INFO << "Initial (#" <<  (num_toks+1) << "): " << token.text;
	  tl.tokenList.push_back(token);
	  num_toks++;
	}
	tok.__set_tokenList(tl);
	sent.__set_tokenization(tok);
	section.sentenceList.push_back(sent);
	section.__isset.sentenceList = true;
      }
      comm.sectionList.push_back(section);
      comm.__isset.sectionList = true;
    }
    int nt = 0;
    concrete::util::StructIterator<concrete::Token, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator, std::vector<concrete::Token>::const_iterator> titerstr(comm);
    auto titer = titerstr.begin();
    auto tend = titerstr.end();
    for(; titer != tend; ++titer) {
      EXPECT_EQ(expected_text.at(nt), titer->text);
      ++nt;
      //INFO << nt << " ::: " << titer->text;
      if(nt > num_toks) break;
    }
    ASSERT_EQ(num_toks, nt) << "iteration " << i << " failed";
  }
}

TEST(StructIterator, Token_deterministic) {
  for(int i = 0; i < 1; ++i) {
    concrete::Communication comm;
    const int num_sects = 3;
    std::vector<int> num_sents(num_sects);
    std::vector< std::vector<int> > num_tokens(num_sects);
    int num_toks = 0;
    for(int isect = 0; isect < num_sects; ++isect) {
      int n_sent = isect == 1 ? 0 : 1;
      num_sents[isect] = n_sent;
      num_tokens[isect] = std::vector<int>(n_sent);
      concrete::Section section;
      for(int isent = 0; isent < n_sent; ++isent) {
	concrete::Sentence sent;
	concrete::Tokenization tok;
	int n_token = 2;
	num_tokens[isect][isent] = n_token;
	concrete::TokenList tl;
	for(int itok = 0; itok < n_token; ++itok) {
	  concrete::Token token;
	  token.__set_text("token" + std::to_string(isect) + "/" + std::to_string(isent) +"/" + std::to_string(itok));
	  //INFO << "Initial: " << token.text;
	  tl.tokenList.push_back(token);
	  num_toks++;
	}
	tok.__set_tokenList(tl);
	sent.__set_tokenization(tok);
	section.sentenceList.push_back(sent);
	section.__isset.sentenceList = true;
      }
      comm.sectionList.push_back(section);
      comm.__isset.sectionList = true;
    }
    int nt = 0;
    concrete::util::StructIterator<concrete::Token, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator, std::vector<concrete::Token>::const_iterator> titerstr(comm);
    auto titer = titerstr.begin();
    auto tend = titerstr.end();
    for(; titer != tend; ++titer) {
      ++nt;
      //INFO << nt << " ::: " << titer->text;
      if(nt > num_toks) break;
    }
    ASSERT_EQ(num_toks, nt);
  }
}

TEST(StructIterator, Token_no_Section) {
  concrete::Communication comm;
  const int num_sects = 0;
  std::vector<int> num_sents(num_sects);
  std::vector< std::vector<int> > num_tokens(num_sects);
  int num_toks = 0;
  int nt = 0;
  concrete::util::StructIterator<concrete::Token, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator, std::vector<concrete::Token>::const_iterator> titerstr(comm);
  auto titer = titerstr.begin();
  auto tend = titerstr.end();
  for(; titer != tend; ++titer) {
    ++nt;
    if(nt > num_toks) break;
  }
  ASSERT_EQ(num_toks, nt);
}

TEST(StructIterator, Token_no_Sentence) {
  concrete::Communication comm;
  const int num_sects = 4;
  std::vector<int> num_sents(num_sects);
  std::vector< std::vector<int> > num_tokens(num_sects);
  int num_toks = 0;
  for(int isect = 0; isect < num_sects; ++isect) {
    int n_sent = 0; //(int)gsl_rng_uniform_int(mathops::rnd_gen, 7);
    num_sents[isect] = n_sent;
    num_tokens[isect] = std::vector<int>(n_sent);
    concrete::Section section;
    comm.sectionList.push_back(section);
    comm.__isset.sectionList = true;
  }
  int nt = 0;
  concrete::util::StructIterator<concrete::Token, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator, std::vector<concrete::Token>::const_iterator> titerstr(comm);
  auto titer = titerstr.begin();
  auto tend = titerstr.end();
  for(; titer != tend; ++titer) {
    ++nt;
    if(nt > num_toks) break;
  }
  ASSERT_EQ(num_toks, nt);
}

TEST(StructIterator, Token_no_Token) {
  concrete::Communication comm;
  const int num_sects = 4;
  std::vector<int> num_sents(num_sects);
  std::vector< std::vector<int> > num_tokens(num_sects);
  int num_toks = 0;
  for(int isect = 0; isect < num_sects; ++isect) {
    int n_sent = 3; //(int)gsl_rng_uniform_int(mathops::rnd_gen, 7);
    num_sents[isect] = n_sent;
    num_tokens[isect] = std::vector<int>(n_sent);
    concrete::Section section;
    for(int isent = 0; isent < n_sent; ++isent) {
      concrete::Sentence sent;
      concrete::Tokenization tok;
      int n_token = 0; //(int)gsl_rng_uniform_int(mathops::rnd_gen, 7);
      //int n_token = 1;
      num_tokens[isect][isent] = n_token;
      concrete::TokenList tl;
      for(int itok = 0; itok < n_token; ++itok) {
	concrete::Token token;
	token.__set_text("token" + std::to_string(isect) + "/" + std::to_string(isent) +"/" + std::to_string(itok));
	tl.tokenList.push_back(token);
	num_toks++;
      }
      tok.__set_tokenList(tl);
      sent.__set_tokenization(tok);
      section.sentenceList.push_back(sent);
      section.__isset.sentenceList = true;
    }
    comm.sectionList.push_back(section);
    comm.__isset.sectionList = true;
  }
  int nt = 0;
  concrete::util::StructIterator<concrete::Token, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator, std::vector<concrete::Token>::const_iterator> titerstr(comm);
  auto titer = titerstr.begin();
  auto tend = titerstr.end();
  for(; titer != tend; ++titer) {
    ++nt;
    //INFO << nt << " ::: " << titer->text;
    if(nt > num_toks) break;
  }
  ASSERT_EQ(num_toks, nt);
}
