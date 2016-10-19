#include "ferrum/data_util.hpp"
#include "ferrum/lock.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/redis_corpus.hpp"

#include <boost/make_shared.hpp>
#include <boost/regex.hpp>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>

#include <chrono>
#include <hiredis/hiredis.h>
#ifndef PUBLIC_UPF
#include <hihiredis/core.h>
#include <hihiredis/redis_server_wrapper.h>
#endif
#include <iostream>
#include <memory>
#include <mutex> // for std::lock_guard
#include <stdexcept>
#include <string>
#include <thread>

namespace ferrum {
  namespace compress {
    namespace gzip {
      std::string compress(const std::string& large) {
	boost::iostreams::filtering_ostream os;
	os.push(boost::iostreams::gzip_compressor{});
	std::string small;
	boost::iostreams::back_insert_device<std::string> snk(small);
	os.push(snk);
	os << large << std::flush;
	os.pop();
	return small;
      }

      std::string decompress(const std::string& small) {
	boost::iostreams::filtering_istream is;
	is.push( boost::iostreams::gzip_decompressor{} );
	boost::iostreams::array_source src( small.c_str(), small.size() );
	std::stringstream ss;
	is.push( src );
	boost::iostreams::copy(is, ss);
	std::string large(ss.str());
	//is >> large;
	is.pop();
	return large;
      }
    };
  };

  namespace db {
    std::ostream & operator<<(std::ostream &os, const Address& addr) {
      return os << "ferrum::db::Address(" << addr.h_ << ":" << addr.p_ << ")";
    }
    Address::Address(const std::string& host, unsigned int port) :
      h_(host), p_(port) {
      if(host.size() > HOST_NAME_MAX) {
	throw std::length_error("Host name size too long");
      }
    }
    const std::string& Address::host() const {
      return h_;
    }
    unsigned int Address::port() const {
      return p_;
    }

    DBError::DBError() : std::runtime_error("") {
    }
    DBError::DBError(const std::string& mess) : std::runtime_error(mess) {
    }
    // protected ctor
    RedisQuery::RedisQuery()
      : next_(NULL),
	rr_(NULL),
	c_(NOT_SET),
	is_set_(false),
	nil_okay_(false) {
    }
    //public ctor
    RedisQuery::RedisQuery(const std::string& key) 
      : next_(NULL),
	rr_(NULL),
	c_(NOT_SET),
	k_(key),
	is_set_(false),
	nil_okay_(false) {
    }
    //public ctor
    RedisQuery::RedisQuery(const std::string& key, const std::string& field) 
      : next_(NULL),
	rr_(NULL),
	c_(NOT_SET),
	k_(key),
	f_(field),
	is_set_(false),
	nil_okay_(false)  {
    }
    //public ctor
    RedisQuery::RedisQuery
    (
     const std::string& key,
     const std::string& field,
     const std::string& value
     ) : next_(NULL),
	 rr_(NULL),
	 c_(NOT_SET),
	 k_(key),
	 f_(field),
	 v_(value),
	 is_set_(true),
	 nil_okay_(false)  {
    }
    //public ctor
    RedisQuery::RedisQuery
    (
     const std::string& key,
     const std::string& field,
     const char* value,
     size_t size
     ) :
      RedisQuery(key, field, std::string(value, size)) {
    }
    RedisQuery::RedisQuery
    (
     const std::string& key,
     const char* fvalue,
     size_t fsize
     ) :
      RedisQuery(key, std::string(fvalue, fsize)) {
    }
    RedisQuery::RedisQuery
    (
     const std::string& key,
     const char* fvalue,
     size_t fsize,
     const char* value,
     size_t size
     ) :
      RedisQuery(key, std::string(fvalue, fsize), std::string(value, size)) {
    }
    RedisQuery::~RedisQuery() {
      if(rr_ != NULL) {
	freeReplyObject(rr_);
      }
      if(next_ != NULL) {
	delete next_;
      }
    }
    RedisQuery::RedisQuery(RedisQuery&& rq) :
      next_(rq.next_),
      rr_(rq.rr_),
      c_(std::move(rq.c_)),
      k_(std::move(rq.k_)),
      f_(std::move(rq.f_)),
      v_(std::move(rq.v_)),
      opts_(std::move(rq.opts_)),
      is_set_(rq.is_set_),
      nil_okay_(rq.nil_okay_) {
    }
    RedisQuery RedisQuery::ping() {
      RedisQuery rq;
      rq.c_ = PING;
      return rq;
    }
    RedisQuery RedisQuery::exists(const std::string& k) {
      RedisQuery rq;
      rq.c_ = EXISTS;
      rq.k_ = k;
      return rq;
    }
    RedisQuery& RedisQuery::operator=(RedisQuery&& rq) {
      if(this != &rq) {
	if(rr_ != NULL) {
	  freeReplyObject(rr_);
	}
	this->rr_ = rq.rr_;
	rq.rr_ = NULL;
	c_ = std::move(rq.c_);
	k_ = std::move(rq.k_);
	f_ = std::move(rq.f_);
	v_ = std::move(rq.v_);
	is_set_ = rq.is_set_;
	nil_okay_ = rq.nil_okay_;
	opts_ = std::move(rq.opts_);
      }
      return *this;
    }
    RedisQuery* RedisQuery::copy_query() const {
      if(next_ != NULL) {
	ERROR << "next_ must be null";
	throw DBError();
      }
      if(rr_ != NULL) {
	ERROR << "rr_ must be null";
	throw DBError();
      }
      RedisQuery* rq = new RedisQuery;
      rq->c_ = this->c_;
      rq->k_ = this->k_;
      rq->f_ = this->f_;
      rq->v_ = this->v_;
      rq->is_set_ = this->is_set_;
      rq->nil_okay_ = this->nil_okay_;
      rq->opts_ = this->opts_;
      return rq;
    }
    RedisQuery& RedisQuery::operator+=(const RedisQuery& rq) {
      RedisQuery** rqn = &next_;
      while( (*rqn) != NULL) {
	rqn = &((*rqn)->next_);
      }
      *rqn = rq.copy_query();
      return *this;
    }
    bool RedisQuery::is_mutating() {
      switch(c_) {
      case HGET:
      case HMGET:
      case HEXISTS:
      case HSET:
      case HDEL:
      case HLEN:
      case PING:
      case LRANGE:
      case LLEN:
      case SMEMBERS:
      case EXISTS:
      case RPUSH:
      case SET:
      case TTL:
      case DEL:
      case HINCRBY:
	return true;
      default:
	return false;
      }
    }
    void RedisQuery::nil_okay(bool b) {
      nil_okay_ = b;
    }
    bool RedisQuery::nil_okay() {
      return nil_okay_;
    }
    RedisQueryOptions::RedisQueryOptions() :
      nx(false), xx(false), ex(false), px(false), ttl(0),
      lrange({ 0, -1, true }) {
    }
    bool RedisQueryOptions::valid() {
      if(nx && xx) return false;
      if(ex && px) return false;
      if(ttl > 0 && !(ex ^ px)) return false;
      return true;
    }
    RedisQuery& RedisQuery::options(const RedisQueryOptions& opts) {
      opts_ = opts;
      return *this;
    }
    void RedisQuery::check_reply(const redisContext* rc, int expected_t) {
      if(rr_ == NULL) {
	ERROR << rc->errstr;
	throw DBError();
      } else if(rr_->type == REDIS_REPLY_ERROR) {
	ERROR << rr_->str;
	throw DBError();
      } else if(rr_->type == REDIS_REPLY_NIL) {
	if(!nil_okay_) {
	  ERROR << "unexpected redis reply: nil";
	  throw DBError();
	}
      } else if(rr_->type != expected_t) {
	ERROR << "unexpected redis reply type " << rr_->type;
	throw DBError();
      }

      if(rr_->type == REDIS_REPLY_NIL) {
	is_set_ = false;
	INFO << "Got nil reply";
      } else {
	if(this->is_mutating()) {
	  is_set_ = true;
	  switch(rr_->type) {
	  case REDIS_REPLY_ARRAY:
	    for(size_t i = 0; i < rr_->elements; ++i) {
	      vec_v_.push_back
		(
		 std::string(rr_->element[i]->str,
			     rr_->element[i]->len)
		 );
	      TRACE << "Finding " << rr_->element[i]->str << " of length " << rr_->element[i]->len;
	    }
	    break;
	  case REDIS_REPLY_INTEGER:
	    v_ = std::to_string(rr_->integer);
	    break;
	  default:
	    v_ = std::string(rr_->str, rr_->len);
	    break;
	  }
	}
      }
    }
    bool RedisQuery::is_set() {
      return is_set_;
    }
    const std::string& RedisQuery::value() {
      return v_;
    }
    const std::vector<std::string>& RedisQuery::vec_value() {
      return vec_v_;
    }
    std::vector<std::string> RedisQuery::values() {
      std::vector<std::string> ret;
      RedisQuery* head = this;
      while(head != NULL && head->queryable()) {
	ret.push_back( head->value() );
	head = head->next_;
      }
      return ret;
    }
    RedisQuery& RedisQuery::hget() {
      c_ = HGET;
      return *this;
    }
    RedisQuery& RedisQuery::hset() {
      c_ = HSET;
      return *this;
    }
    RedisQuery& RedisQuery::hlen() {
      c_ = HLEN;
      return *this;
    }
    RedisQuery& RedisQuery::hdel() {
      c_ = HDEL;
      return *this;
    }
    RedisQuery& RedisQuery::hexists() {
      c_ = HEXISTS;
      return *this;
    }
    RedisQuery& RedisQuery::hmget() {
      c_ = HMGET;
      return *this;
    }
    RedisQuery& RedisQuery::hmset() {
      c_ = HMSET;
      return *this;
    }
    RedisQuery& RedisQuery::rpush() {
      c_ = RPUSH;
      return *this;
    }
    RedisQuery& RedisQuery::lpush() {
      c_ = LPUSH;
      return *this;
    }
    RedisQuery& RedisQuery::lrange() {
      c_ = LRANGE;
      return *this;
    }
    RedisQuery& RedisQuery::llen() {
      c_ = LLEN;
      return *this;
    }
    RedisQuery& RedisQuery::smembers() {
      c_ = SMEMBERS;
      return *this;
    }
    RedisQuery& RedisQuery::ttl() {
      c_ = TTL;
      return *this;
    }
    RedisQuery& RedisQuery::get() {
      c_ = GET;
      return *this;
    }
    RedisQuery& RedisQuery::set() {
      c_ = SET;
      return *this;
    }
    RedisQuery& RedisQuery::del() {
      c_ = DEL;
      return *this;
    }
    RedisQuery& RedisQuery::hincrby() {
      c_ = HINCRBY;
      return *this;
    }
    void RedisQuery::operator()(redisContext* rc) {
      int et = -1;
      switch(c_) {
      case HGET: 
	rr_ = static_cast<redisReply*>(redisCommand(rc, "HGET %s %s", k_.c_str(), f_.c_str()));
	et = REDIS_REPLY_STRING;
	break;
      case HSET:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "HSET %s %s %b", k_.c_str(), f_.c_str(), v_.c_str(), v_.size()));
	et = REDIS_REPLY_INTEGER;
	break;
      case HDEL:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "HDEL %s %b", k_.c_str(), f_.c_str(), f_.size()));
	et = REDIS_REPLY_INTEGER;
	break;
      case HLEN:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "HLEN %s", k_.c_str()));
	et = REDIS_REPLY_INTEGER;
	break;
      case HEXISTS:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "HEXISTS %s %b", k_.c_str(), f_.c_str(), f_.size()));
	et = REDIS_REPLY_INTEGER;
	break;
      case HMGET: 
	rr_ = static_cast<redisReply*>(redisCommand(rc, "HMGET %s %s", k_.c_str(), f_.c_str()));
	et = REDIS_REPLY_STRING;
	break;
      case HMSET:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "HMSET %s %s %b", k_.c_str(), f_.c_str(), v_.c_str(), v_.size()));
	et = REDIS_REPLY_STATUS;
	break;
      case SMEMBERS:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "SMEMBERS %s %b", k_.c_str(), k_.size()));
	et = REDIS_REPLY_ARRAY;
	break;
      case LPUSH:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "LPUSH %s %b", k_.c_str(), f_.c_str(), f_.size()));
	et = REDIS_REPLY_INTEGER;
	break;
      case RPUSH:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "RPUSH %s %b", k_.c_str(), f_.c_str(), f_.size()));
	et = REDIS_REPLY_INTEGER;
	break;
      case LLEN:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "LLEN %s", k_.c_str()));
	et = REDIS_REPLY_INTEGER;
	break;
      case LRANGE:
	//rr_ = static_cast<redisReply*>(redisCommand(rc, "LRANGE %s %b", k_.c_str(), f_.c_str(), f_.size()));

	et = REDIS_REPLY_ARRAY;
	{
	  std::string cmd_str="LRANGE %s ";
	  cmd_str += std::to_string(opts_.lrange.start);
	  cmd_str += " ";
	  int end = opts_.lrange.inclusive ? (opts_.lrange.end) : (opts_.lrange.end - 1);
	  cmd_str += std::to_string(end);
	  //INFO << "command [" << cmd_str << "]";
	  rr_ = static_cast<redisReply*>(redisCommand(rc, cmd_str.c_str(), k_.c_str()));
	}
	break;
      case EXISTS:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "EXISTS %s", k_.c_str()));
	et = REDIS_REPLY_INTEGER;
	break;
      case TTL:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "TTL %s", k_.c_str()));
	et = REDIS_REPLY_INTEGER;
	break;	
      case DEL:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "DEL %s", k_.c_str()));
	et = REDIS_REPLY_INTEGER;
	break;	
      case HINCRBY:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "HINCRBY %s %s %s", k_.c_str(), f_.c_str(), v_.c_str()));
	et = REDIS_REPLY_INTEGER;
	break;	
      case SET:
	{
	  std::string cmd_str="SET %s %b";
	  if(!opts_.valid()) {
	    ERROR << "Not valid options given";
	    throw DBError();
	  }
	  if(opts_.nx) cmd_str += " NX";
	  else if(opts_.xx) cmd_str += " XX";
	  if(opts_.ex) cmd_str += " EX %d";
	  else if(opts_.px) cmd_str += " PX %d";
	  if(opts_.ttl > 0) {
	    rr_ = static_cast<redisReply*>(redisCommand(rc, cmd_str.c_str(), k_.c_str(), f_.c_str(), f_.size(), opts_.ttl));
	  } else {
	    rr_ = static_cast<redisReply*>(redisCommand(rc, cmd_str.c_str(), k_.c_str(), f_.c_str(), f_.size()));
	  }
	}
	et = REDIS_REPLY_STATUS;
	break;	
      case GET:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "GET %s", k_.c_str()));
	et = REDIS_REPLY_STRING;
	break;
      case PING:
	rr_ = static_cast<redisReply*>(redisCommand(rc, "PING"));
	et = REDIS_REPLY_STATUS;
	break;	
      default:
	ERROR << "Unknown command " << c_;
	throw DBError();
	break;
      }
      this->check_reply(rc, et);
    }
    bool RedisQuery::queryable() {
      return c_ != NOT_SET;
    }
    void RedisQuery::query_all(redisContext* rc) {
      RedisQuery* head = this;
      while(head != NULL && head->queryable()) {
	head->operator()(rc);
	head = head->next_;
      }
    }
    int RedisQuery::num_queries() {
      RedisQuery* head = this;
      int nq = 0;
      while(head != NULL && head->queryable()) {
	++nq;
	head = head->next_;
      }
      return nq;
    }

#ifndef PUBLIC_UPF
    Redis::RedisStarter::RedisStarter() : host(NULL), my_argv(NULL) {
      INFO << "RedisStarter()";
      const char *argv[2] = {"redis-server", NULL};
      my_argv = hhr_build_args(argv, &host, &rs_port);
      rs_pid = fork();
      if(rs_pid < 0) {
	ERROR << "Could not fork";
      } else {	
	INFO << __PRETTY_FUNCTION__ << " is dealing with " << ( rs_pid == 0 ? "child" : "parent" );
      }
      int ret = hhr_redis_start_post_fork(my_argv, rs_pid, rs_port, 5);
      INFO << "pid of " << rs_pid << " trying to connect";
      if (ret != HHR_OK) {
	ERROR << "failed to start redis " << ret;
	throw ret;
      }
    }
    Redis::RedisStarter::~RedisStarter() {
      hhr_free_args(my_argv, host);
    }
#endif

    MultiClientRedisLock::MultiClientRedisLock(const std::string& value, std::shared_ptr<Redis> db, unsigned int timeout, bool print_ok_status) :
      val_(value), db_(db), acquired_(false), pos_(print_ok_status) {
      RedisQuery lock("mcrl", val_);
      RedisQueryOptions opts;
      opts.nx=true;
      opts.px=true;
      opts.ttl = timeout;
      lock.set().options(opts).nil_okay(true);
      db_->operator()(lock);
      int iter = 2;
      if(pos_) {
	INFO << "db lock acquired for mcrl @ " << val_ << ": ? " << lock.value();
      }
      while(! lock.is_set() ) {
	INFO << "Attempt #" << iter << " to acquire a multi-client lock via key mcrl with value \"" << val_ << "\"";
	std::chrono::seconds timespan(5);
	std::this_thread::sleep_for(timespan);
	db_->operator()(lock);
	++iter;
      }
      if(iter > 2) {
	INFO << "Successfully acquired lock after " << (iter-1) << " attempts (key mcrl with value \"" << val_ << "\")";
      }
      acquired_ = true;
    }
    MultiClientRedisLock::MultiClientRedisLock(const std::string& value, std::shared_ptr<Redis> db) :
      MultiClientRedisLock(value, db, DEFAULT_TIMEOUT, true) {
    }
    MultiClientRedisLock::MultiClientRedisLock(const std::string& value, std::shared_ptr<Redis> db, unsigned int to) :
      MultiClientRedisLock(value, db, to, true) {
    }
    MultiClientRedisLock::MultiClientRedisLock(const std::string& value, std::shared_ptr<Redis> db, bool print_if_okay) :
      MultiClientRedisLock(value, db, DEFAULT_TIMEOUT, print_if_okay) {
    }
    MultiClientRedisLock::~MultiClientRedisLock() {
      // get the ttl
      int ttl = db_->ttl("mcrl");
      switch(ttl) {
      case -2: // key doesn't exist
	if(acquired_) {
	  INFO << "Key mcrl @ " << val_ << " no longer exists; perhaps increase the timeout?";
	}
	break;
      case -1: // key exists, no timeout
	WARN << "Key mcrl @ " << val_ << "exists, but without a timeout. This indicates a lack of uniqueness in multiclient locking.";
	break;
      default:
	{
	  RedisQuery del("mcrl");
	  db_->operator()(del.del());
	  int num_deleted = std::stoi(del.value());
	  switch(num_deleted) {
	  case 1:
	    if(pos_) {
	      INFO << "Successfully unlocked mcrl @ " << val_;
	    }
	    break;
	  case 0:
	    WARN << "Did not remove mcrl lock @ " << val_;
	    break;
	  default:
	    WARN << "Removed unexpected number of keys, " << num_deleted << ", in unlocking mcrl @ " << val_;
	  }
	}
	break;
      }
    }
    bool MultiClientRedisLock::locked() {
      int ttl = db_->ttl("mcrl");
      switch(ttl) {
      case -2: // key doesn't exist
	acquired_ = false;
	break;
      case -1: // key exists, no timeout
	WARN << "Key mcrl @ " << val_ << "exists, but without a timeout. This indicates a lack of uniqueness in multiclient locking.";
	acquired_ = false;
	break;
      default:
	acquired_ = true;
	break;
      }
      return acquired_;
    }

    const unsigned int MultiClientRedisLock::DEFAULT_TIMEOUT = 1000000;

#ifndef PUBLIC_UPF
    Redis::Redis() : context_(NULL), port_(0), pid_(0), in_situ_(false), mt_(false), mut_() {
      RedisStarter starter;
      in_situ_ = true;
      port_ = starter.rs_port;
      pid_ = starter.rs_pid;
      context_ = connect({"localhost", (unsigned int) port_});
    }
#endif
    Redis::Redis(const Address& addr) : context_(NULL), in_situ_(false), mt_(false), mut_() {
      context_ = connect(addr);
    }
    Redis::~Redis() {
      if(context_ != NULL) {
	disconnect(context_);
      }
#ifndef PUBLIC_UPF
      if(in_situ_) {
	INFO << "Stopping the in situ redis server";
	int ret = hhr_redis_stop(port_, pid_);
	if (ret != HHR_OK) {
	  ERROR << "failed to stop redis PID " << pid_ << " with return value " << ret;
	  throw ret;
	}
      }
#endif
    }
    typename Redis::context_t Redis::connect(const Address& addr) {
      Redis::context_t context = redisConnect(addr.host().c_str(), addr.port());

      if (context == NULL) {
	ERROR << "redis::connect can't allocate redis context to address " << addr;
	throw DBError();
      } else if (context->err) {
	ERROR << "redis::connect encountered error " << context->errstr;
	throw DBError();
      }
      return context;
    }
    void Redis::disconnect(typename Redis::context_t context) {
      redisFree(context);
    }
    void Redis::operator()(RedisQuery& rq) {
      if(mt_) {
	//ferrum::Lock lock(mut_);
	std::lock_guard<ferrum::Mutex> lock(mut_);
	rq.query_all(context_);
      } else {
	rq.query_all(context_);
      }
    }
    void Redis::multithreaded(bool b) {
      mt_ = b;
      if(mt_) {
	mut_ = ferrum::Mutex();
      }
    }
    bool Redis::multithreaded() const {
      return mt_;
    }
    bool Redis::multithreaded() {
      return const_cast<const Redis*>(this)->multithreaded();
    }

    bool Redis::has_key(const std::string& key) {
      RedisQuery ex = RedisQuery::exists(key);
      this->operator()(ex);
      bool res = false;
      const std::string& val = ex.value();
      if(val == "1") {
	res = true;
      } else if(val == "0") {
      } else {
	ERROR << "Unknown return value " << val;
	throw DBError();
      }
      return res;
    }

    int Redis::llen(const std::string& key) {
      RedisQuery d(key);
      this->operator()(d.llen());
      return std::stoi(d.value());
    }

    int Redis::del(const std::string& key) {
      RedisQuery d(key);
      this->operator()(d.del());
      return std::stoi(d.value());
    }

    int Redis::ttl(const std::string& key) {
      RedisQuery rq_ttl(key);
      rq_ttl.ttl();
      this->operator()(rq_ttl);
      return std::stoi(rq_ttl.value());
    }

    int Redis::rpush(const std::string& key, const std::string& val) {
      RedisQuery rq(key, val);
      this->operator()(rq.rpush());
      return std::stoi(rq.value());
    }

    int Redis::hexists(const std::string& key, const std::string& field) {
      RedisQuery rq(key, field);
      this->operator()(rq.hexists());
      return std::stoi(rq.value());
    }

    int Redis::hset(const std::string& key, const std::string& field, const std::string& val) {
      RedisQuery rq(key, field.c_str(), field.size(), val.c_str(), val.size());
      this->operator()(rq.hset());
      return std::stoi(rq.value());
    }

    std::string Redis::hget(const std::string& key, const std::string& field) {
      RedisQuery rq(key, field.c_str(), field.size());
      this->operator()(rq.hget());
      return rq.value();
    }

    int Redis::hlen(const std::string& key) {
      RedisQuery rq(key);
      this->operator()(rq.hlen());
      return std::stoi(rq.value());
    }

    std::vector<std::string> Redis::smembers(const std::string& key) {
      RedisQuery rq(key);
      this->operator()(rq.smembers());
      return rq.vec_value();
    }

    std::shared_ptr<Redis> Redis::get_shared_ptr() {
      //return std::enable_shared_from_this<ferrum::db::Redis>::shared_from_this();
      return shared_from_this();
    }

    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////

    template <typename P>
    RedisThriftSmartWriter<P>::RedisThriftSmartWriter(const std::string& name, std::shared_ptr<ferrum::db::Redis> redis_db) :
      ferrum::thrift::ThriftSmartWriter<P>(name),
      redis_(redis_db),
      transport_(new MemoryBuffer()) {
      this->protocol_ = boost::make_shared<P>(transport_);
    }
    template <typename P>
    RedisThriftSmartWriter<P>::RedisThriftSmartWriter(const std::string& name) :
      ferrum::thrift::ThriftSmartWriter<P>(name),
      transport_(new MemoryBuffer()) {
      this->protocol_ = boost::make_shared<P>(transport_);
    }

    template <typename P>
    RedisThriftSmartWriter<P>::~RedisThriftSmartWriter() {
    }

    template <typename P>
    void RedisThriftSmartWriter<P>::redis(std::shared_ptr<ferrum::db::Redis> redis_db) {
      redis_ = redis_db;
    }

    template <typename P>
    P* RedisThriftSmartWriter<P>::get(const std::string& suffix) {
      //this->curr_file = ferrum::SmartWriter::next_file_name(suffix);
      this->curr_file = suffix;
      return this->protocol_.get();
    }

    template <typename P>
    void RedisThriftSmartWriter<P>::_save() {
      std::string saved_struct( transport_->getBufferAsString() );
      RedisQuery rq
	(
	 this->base_,
	 this->curr_file, // field
	 saved_struct.c_str(),
	 saved_struct.size()
	 );
      rq.hset();
      (*redis_)(rq);
      transport_->resetBuffer();
    }

    template <typename P>
    std::string RedisThriftSmartWriter<P>::form_id(const std::string& partial) {
      return ferrum::create_doc_hset_id(partial);
    }

    template <typename P>
    RedisThriftSmartWriter<P>* RedisThriftSmartWriter<P>::clone() const {
      RedisThriftSmartWriter<P>* rtsw = new RedisThriftSmartWriter<P>(this->base_, this->redis_);
      return rtsw;
    }
  } // end namespace db

    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////

  template <typename P>
  TarThriftSmartWriter<P>::TarThriftSmartWriter(const std::string& name, std::shared_ptr<ferrum::CompressedTar> archive) :
    ferrum::thrift::ThriftSmartWriter<P>(name),
    archive_(archive),
    transport_(new MemoryBuffer()) {
    this->protocol_ = boost::make_shared<P>(transport_);
  }
  template <typename P>
  TarThriftSmartWriter<P>::TarThriftSmartWriter(const std::string& name) :
    ferrum::thrift::ThriftSmartWriter<P>(name),
    transport_(new MemoryBuffer()) {
    this->protocol_ = boost::make_shared<P>(transport_);
  }

  template <typename P>
  TarThriftSmartWriter<P>::~TarThriftSmartWriter() {
  }

  // template <typename P>
  // void RedisThriftSmartWriter<P>::redis(std::shared_ptr<ferrum::db::Redis> redis_db) {
  //   redis_ = redis_db;
  // }

  template <typename P>
  P* TarThriftSmartWriter<P>::get(const std::string& suffix) {
    //this->curr_file = ferrum::SmartWriter::next_file_name(suffix);
    this->curr_file = suffix;
    return this->protocol_.get();
  }

  template <typename P>
  std::string TarThriftSmartWriter<P>::safe_path(const std::string& in) {
    boost::regex re("[^a-zA-Z0-9_\\-\\.]");
    std::string out = boost::regex_replace(in, re, "_");
    return out;
  }

  template <typename P>
  void TarThriftSmartWriter<P>::_save() {
    std::string saved_struct( transport_->getBufferAsString() );
    std::stringstream archive_path;
    archive_path << this->base_ << "/";
    archive_path << this->curr_file;
    archive_->write_data
      (
       safe_path(archive_path.str()),
       &(saved_struct[0]),
       saved_struct.size()
       );
    transport_->resetBuffer();
  }

  template <typename P>
  std::string TarThriftSmartWriter<P>::form_id(const std::string& partial) {
    return ferrum::create_doc_hset_id(partial);
  }

  template <typename P>
  TarThriftSmartWriter<P>* TarThriftSmartWriter<P>::clone() const {
    TarThriftSmartWriter<P>* rtsw = new TarThriftSmartWriter<P>(this->base_, this->archive_);
    return rtsw;
  }
} // end namespace ferrum

template class ferrum::db::RedisThriftSmartWriter<ferrum::thrift::TBinaryProtocol>;
template class ferrum::db::RedisThriftSmartWriter<ferrum::thrift::TCompactProtocol>;
template class ferrum::db::RedisThriftSmartWriter<ferrum::thrift::TJSONProtocol>;
template class ferrum::TarThriftSmartWriter<ferrum::thrift::TBinaryProtocol>;
template class ferrum::TarThriftSmartWriter<ferrum::thrift::TCompactProtocol>;
template class ferrum::TarThriftSmartWriter<ferrum::thrift::TJSONProtocol>;
