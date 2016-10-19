#include "gtest/gtest.h"

#include "ferrum/concrete.hpp"
#include "ferrum/data_util.hpp"
#include "ferrum/redis_corpus.hpp"

TEST(RedisCorpus, Communication) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  ferrum::RedisCorpus<concrete::Communication> rc({"localhost", 4532});
}

TEST(RedisCorpus, Communication_num_mentions) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  {
    ferrum::RedisCorpus<concrete::Communication> rc({"localhost", 4532});
    concrete::Communication comm;
    rc.add_document(comm);
    ASSERT_EQ(0, rc.num_mentions());
  }
  {
    ferrum::RedisCorpus<concrete::Communication> rc({"localhost", 4532});
    concrete::Communication comm;
    concrete::EntityMention em;
    concrete::EntityMentionSet ems;
    ems.mentionList.push_back(em);
    comm.entityMentionSetList.push_back(ems);
    comm.__isset.entityMentionSetList = true;
    rc.add_document(comm);
    ASSERT_EQ(1, rc.num_mentions());
  }
}

TEST(RedisCorpus, SimpleDoc) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  ferrum::RedisCorpus<minsky::SimpleDoc> rc({"localhost", 4532});
  {
    ferrum::RedisCorpus<concrete::Communication> rc({"localhost", 4532});
    concrete::Communication comm;
    rc.add_document(comm);
    ASSERT_EQ(0, rc.num_mentions());
  }
}

TEST(RedisCorpus, empty) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  ferrum::RedisCorpus<minsky::SimpleDoc> rc({"localhost", 4532});
  for(const auto it : rc) {
  }
}
