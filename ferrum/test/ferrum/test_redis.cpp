#include "gtest/gtest.h"

#include "ferrum/data_util.hpp"
#include <chrono>
#include <thread>

TEST(DBError, create) {
  ferrum::db::DBError err;
}

TEST(Address, create) {
  ferrum::db::Address addr("localhost", 1447);
  ASSERT_EQ("localhost", addr.host());
  ASSERT_EQ(1447, addr.port());
}

void accept_implicit_address(const ferrum::db::Address& addr) {
  ASSERT_EQ("localhost", addr.host());
  ASSERT_EQ(1447, addr.port());
}
TEST(Address, create_implicit) {
  accept_implicit_address({"localhost", 1447});
}

TEST(Redis, DISABLED_pingpong) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  ferrum::db::RedisQuery ping = ferrum::db::RedisQuery::ping();
  redis(ping);
  ASSERT_TRUE(ping.is_set());
  ASSERT_EQ("PONG", ping.value());
  ASSERT_EQ(1, ping.num_queries());
}

TEST(Redis, DISABLED_multi_pingpong) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  ferrum::db::RedisQuery ping = ferrum::db::RedisQuery::ping();
  ASSERT_EQ(1, ping.num_queries());
  ping += ferrum::db::RedisQuery::ping();
  ASSERT_EQ(2, ping.num_queries());
  ping += ferrum::db::RedisQuery::ping();
  ASSERT_EQ(3, ping.num_queries());
  ping += ferrum::db::RedisQuery::ping();
  ASSERT_EQ(4, ping.num_queries());
  redis(ping);
  ASSERT_TRUE(ping.is_set());
  ASSERT_EQ("PONG", ping.value());
  std::vector<std::string> vals = ping.values();
  ASSERT_EQ(4, vals.size());
  for(const std::string v : vals) {
    ASSERT_EQ("PONG", v);
  }
}

TEST(Redis, DISABLED_store_nullterminator) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  ferrum::db::RedisQuery rq("nulltest", "val", "foo\0bar", 7);
  rq.hset();
  redis(rq);
  ferrum::db::RedisQuery rq_getter("nulltest", "val");
  rq.hget();
  redis(rq);
  ASSERT_EQ(7, rq.value().size());
}

TEST(Redis, DISABLED_insitu_pingpong) {
  ferrum::db::Redis redis;
  ferrum::db::RedisQuery ping = ferrum::db::RedisQuery::ping();
  redis(ping);
  ASSERT_TRUE(ping.is_set());
  ASSERT_EQ("PONG", ping.value());
  ASSERT_EQ(1, ping.num_queries());
}

TEST(Redis, DISABLED_insitu_set_int) {
  ferrum::db::Redis redis;
  ferrum::db::RedisQuery set("x", "value", "1");
  set.hset();
  redis(set);
  ferrum::db::RedisQuery get("x", "value");
  get.hget();
  redis(get);
  ASSERT_TRUE(get.is_set());
  ASSERT_EQ("1", get.value());
  ASSERT_EQ(1, get.num_queries());
}

TEST(Redis, lrange) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  {
    ferrum::db::RedisQuery push("my_list", "0");
    push.rpush();
    for(size_t i = 1; i < 5; ++i) {
      ferrum::db::RedisQuery push1("my_list", std::to_string(i));
      push1.rpush();  
      push += push1;
    }
    ASSERT_EQ(5, push.num_queries());
    redis(push);
  }
  {
    ferrum::db::RedisQuery lrange("my_list");
    lrange.lrange();
    redis(lrange);
    ASSERT_TRUE(lrange.is_set());
    ASSERT_EQ(0, lrange.vec_value().size() % 5);
    for(size_t i = 0; i < 5; ++i) {
      ASSERT_EQ(std::to_string(i), lrange.vec_value()[i]);
    }
  }
  {
    ferrum::db::RedisQuery lrange("my_list");
    ferrum::db::RedisQueryOptions opts;
    opts.lrange.start = 0;
    opts.lrange.end = 3;
    opts.lrange.inclusive = false;
    lrange.options(opts);
    lrange.lrange();
    redis(lrange);
    ASSERT_TRUE(lrange.is_set());
    ASSERT_EQ(opts.lrange.end, lrange.vec_value().size());
    for(size_t i = 0; i < (size_t)opts.lrange.end; ++i) {
      ASSERT_EQ(std::to_string(i), lrange.vec_value()[i]);
    }
  }
  {
    ferrum::db::RedisQuery lrange("my_list");
    ferrum::db::RedisQueryOptions opts;
    opts.lrange.start = 0;
    opts.lrange.end = 3;
    opts.lrange.inclusive = true;
    lrange.options(opts);
    lrange.lrange();
    redis(lrange);
    ASSERT_TRUE(lrange.is_set());
    ASSERT_EQ(opts.lrange.end + 1, lrange.vec_value().size());
    for(size_t i = 0; i <= (size_t)opts.lrange.end; ++i) {
      ASSERT_EQ(std::to_string(i), lrange.vec_value()[i]);
    }
  }
  ASSERT_EQ(1, redis.del("my_list"));
}

TEST(Redis, key_exists) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  std::string fname(__func__);
  std::string key("my_list_key_exists:" + fname);
  {
    ferrum::db::RedisQuery push(key, "0");
    push.rpush();
    redis(push);
  }
  ASSERT_TRUE( redis.has_key(key) );
  ASSERT_FALSE( redis.has_key(key + ":not_there") );
  ASSERT_EQ(1, redis.del(key));
}

TEST(Redis, hexists) {
  ferrum::db::Address addr("localhost", 4532);
  ferrum::db::Redis redis(addr);
  std::string fname(__func__);
  std::string key("my_list_key_exists:" + fname);
  redis.hset(key, "field", "value");
  ferrum::db::RedisQuery rq(key, "field");
  rq.hexists();
  redis(rq);
  ASSERT_TRUE( rq.is_set() );
  ASSERT_TRUE( std::stoi(rq.value()) );
  ASSERT_TRUE( redis.hexists(key, "field") );
  ASSERT_EQ(1, redis.del(key));
}

TEST(Redis, multiclientlock) {
  ferrum::db::Address addr("localhost", 4532);
  std::shared_ptr<ferrum::db::Redis> redis(new ferrum::db::Redis(addr));
  typedef ferrum::db::MultiClientRedisLock MCRL;
  std::string pf(__PRETTY_FUNCTION__);
  std::string val1(pf + "1");
  MCRL mcrl(val1, redis, (unsigned int)3000);
  ASSERT_TRUE(mcrl.locked());
  std::string val2(pf + "2");
  MCRL mcrl2(val2, redis, (unsigned int)3000);
  std::chrono::seconds timespan(4);
  std::this_thread::sleep_for(timespan);
  ASSERT_FALSE(mcrl.locked());
  ASSERT_FALSE(mcrl2.locked());
}

