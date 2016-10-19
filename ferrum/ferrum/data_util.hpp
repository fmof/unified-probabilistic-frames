#ifndef FERRUM_DATA_UTIL_H_
#define FERRUM_DATA_UTIL_H_

#include "ferrum/lock.hpp"
#include "ferrum/tar.hpp"
#include "ferrum/thrift_protocol_defs.hpp"
#include "ferrum/thrift_smart_writer.hpp"

#include <boost/shared_ptr.hpp>
#include <hiredis/hiredis.h>
#ifndef PUBLIC_UPF
#include <hihiredis/core.h> // for port_t
#else
typedef unsigned short port_t;
#endif
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <thrift/transport/TTransportUtils.h>

#ifndef HOST_NAME_MAX
#define HOST_NAME_MAX 255
#endif

namespace ferrum {
  namespace thrift {

    template <typename P, typename Struct>
    inline std::string thrift_struct_to_string(const Struct& s) {
      typedef apache::thrift::transport::TMemoryBuffer MemoryBuffer;
      boost::shared_ptr<MemoryBuffer> transport(new MemoryBuffer());
      boost::shared_ptr<P> protocol(new P(transport));
      // I *could* call transport->open(), but for TMemoryBuffer, it's a no-op
      s.write(protocol.get());
      return transport->getBufferAsString();
    };
    template <typename P, typename Struct>
    inline Struct thrift_struct_from_string(const std::string& str) {
      typedef apache::thrift::transport::TMemoryBuffer MemoryBuffer;
      boost::shared_ptr<MemoryBuffer> transport
	(
	 new MemoryBuffer
	 (
	  reinterpret_cast<uint8_t*>
	  (
	   const_cast<char*>(&(str[0]))
	   ),
	  str.size() //size of the buffer
	  )
	 );
      boost::shared_ptr<P> protocol(new P(transport));
      // I *could* call transport->open(), but for TMemoryBuffer, it's a no-op
      Struct s;
      s.read(protocol.get());
      return s;
    };
    template <typename P, typename Struct>
    inline void thrift_struct_from_string(const std::string& str, Struct* s) {
      typedef apache::thrift::transport::TMemoryBuffer MemoryBuffer;
      boost::shared_ptr<MemoryBuffer> transport
	(
	 new MemoryBuffer
	 (
	  reinterpret_cast<uint8_t*>
	  (
	   const_cast<char*>(&(str[0]))
	   ),
	  str.size() //size of the buffer
	  )
	 );
      boost::shared_ptr<P> protocol(new P(transport));
      // I *could* call transport->open(), but for TMemoryBuffer, it's a no-op
      s->read(protocol.get());
    };
    template <typename P, typename Struct>
    inline void thrift_struct_from_buffer(const void* buffer, size_t size, Struct* s) {
      typedef apache::thrift::transport::TMemoryBuffer MemoryBuffer;
      boost::shared_ptr<MemoryBuffer> transport
	(
	 new MemoryBuffer
	 (
	  static_cast<uint8_t*>
	  (
	   const_cast<void*>(buffer)
	   ),
	  size // how much to read
	  )
	 );
      boost::shared_ptr<P> protocol(new P(transport));
      // I *could* call transport->open(), but for TMemoryBuffer, it's a no-op
      s->read(protocol.get());
    };
    
    template <typename Protocol>
    struct FromBuffer {
      template <typename CStruct>
      void operator()(void* buffer, size_t size, CStruct* obj) {
	thrift_struct_from_buffer<Protocol>(buffer, size, obj);
      }
    };
  };

  namespace compress {
    namespace gzip {
      std::string compress(const std::string&);
      std::string decompress(const std::string&);
    };
  };

  namespace db {
    class Address {
    private:
      std::string h_;
      unsigned int p_;
      friend std::ostream & operator<<(std::ostream &os, const Address& addr);
    public:
      Address(const std::string& host, unsigned int port);
      const std::string& host() const;
      unsigned int port() const;
    };

    class DBError : std::runtime_error {
    public:
      DBError();
      DBError(const std::string& mess);
    };

    template <typename DBT>
    class Database {
    public:
      typedef DBT context_t;
      virtual context_t connect(const Address& address) = 0;
      virtual void disconnect(context_t context) = 0;
    };

    enum RedisCommand {
      NOT_SET = -1,
      HGET = 0,
      HMGET,
      HSET,
      HMSET,
      HLEN,
      HEXISTS,
      HDEL,
      LPUSH,
      RPUSH,
      LRANGE,
      LLEN,
      SMEMBERS,
      EXISTS,
      SET,
      GET,
      TTL,
      DEL,
      HINCRBY,
      PING = 1000
    };

    class RedisQueryOptions {
    public:
      bool nx; // create only if not created
      bool xx; // ???
      bool ex; // set time-out in seconds
      bool px; // set time-out in milliseconds
      unsigned int ttl; // time-out
      struct {
	int start;
	int end;
	/**
	 * is this range [start, end] or [start, end) 
	 * default: True (inclusive --> [start, end])
	 */
	bool inclusive;
      } lrange;
      RedisQueryOptions();
      bool valid();
    };

    class RedisQuery {
    protected:
      RedisQuery* next_;
      redisReply* rr_;
      RedisCommand c_; // command
      std::string k_; // key
      std::string f_; // field
      std::string v_; // value
      RedisQueryOptions opts_; // any options for (certain) commands
      std::vector<std::string> vec_v_; // reply_array values
      bool is_set_;
      bool nil_okay_;
      RedisQuery();
      RedisQuery* copy_query() const;
      void check_reply(const redisContext* rc, int expected_t);
      bool is_mutating();
    public:
      RedisQuery(const std::string& key);
      RedisQuery(const std::string& key, const std::string& field);
      RedisQuery(const std::string& key, const std::string& field, const std::string& value);
      RedisQuery(const std::string& key, const std::string& field, const char* value, size_t size);
      RedisQuery(const std::string& key, const char* fvalue, size_t fsize);
      RedisQuery(const std::string& key, const char* fvalue, size_t fsize, const char* value, size_t vsize);
      ~RedisQuery();
      RedisQuery(const RedisQuery& rq) = delete;
      RedisQuery& operator=(const RedisQuery& rq) = delete;
      RedisQuery(RedisQuery&& rq);
      RedisQuery& operator=(RedisQuery&& rq);
      // for pipelining
      RedisQuery& operator+=(const RedisQuery& rq);
      static RedisQuery ping();
      static RedisQuery exists(const std::string& key);
      RedisQuery& options(const RedisQueryOptions& rqo);
      RedisQuery& hget();
      RedisQuery& hset();
      RedisQuery& hexists();
      RedisQuery& hlen();
      RedisQuery& hdel();
      RedisQuery& hmget();
      RedisQuery& hmset();
      RedisQuery& lpush();
      RedisQuery& rpush();
      RedisQuery& lrange();
      RedisQuery& llen();
      RedisQuery& smembers();
      RedisQuery& ttl();
      RedisQuery& get();
      RedisQuery& set();
      RedisQuery& del();
      RedisQuery& hincrby();
      bool is_set();
      bool nil_okay();
      void nil_okay(bool b);
      const std::string& value();
      const std::vector<std::string>& vec_value();
      std::vector<std::string> values();
      void operator()(redisContext* rc);
      bool queryable();
      void query_all(redisContext* rc);
      int num_queries();
    };

    template <typename V>
    class PrimFromStr {
    public:
      V operator()(const std::string&);      
    };

    class Redis : public Database< redisContext* >, public std::enable_shared_from_this< ::ferrum::db::Redis> {
    public:
#ifndef PUBLIC_UPF
      Redis();
#endif
      Redis(const Address& address);
      Redis(const Redis& other) = delete;
      Redis& operator=(const Redis& other) = delete;
      virtual ~Redis();
      virtual context_t connect(const Address& address);
      virtual void disconnect(context_t context);
      void operator()(RedisQuery& rq);
      void multithreaded(bool mt);
      bool multithreaded() const;
      bool multithreaded();
      bool has_key(const std::string& key);
      int del(const std::string& key);
      int ttl(const std::string& key);
      int llen(const std::string& key);
      template <typename V, typename Conv> V list_item(const std::string& key, int idx, Conv converter);
      template <typename V, typename Conv> void list_item(const std::string& key, int idx, V* ret_val, Conv converter);
      int rpush(const std::string&, const std::string&);
      int hexists(const std::string&, const std::string&);
      int hset(const std::string&, const std::string&, const std::string&);
      std::string hget(const std::string&, const std::string&);
      int hlen(const std::string&);
      std::vector<std::string> smembers(const std::string&);
      //template <typename Prim> using PrimToStr = std::string (*)(Prim);
      template <typename V, typename Converter = PrimFromStr<V> > V hget(const std::string&, const std::string&, Converter converter = Converter());
      // template <typename V, typename Converter = PrimToStr<V> > int hset(const std::string&, const std::string&, const V&);
      std::shared_ptr<Redis> get_shared_ptr();

    protected:
      context_t context_;
      port_t port_;
      pid_t pid_;
      bool in_situ_;
      bool mt_; // multithreaded or not
      ferrum::Mutex mut_;

#ifndef PUBLIC_UPF
      class RedisStarter {
      private:
	char *host;
	char **my_argv;
      public:
	RedisStarter();
	~RedisStarter();
	RedisStarter(const RedisStarter&) = delete;
	RedisStarter(RedisStarter&&) = delete;
	RedisStarter& operator=(const RedisStarter&) = delete;
	RedisStarter& operator=(RedisStarter&&) = delete;

	pid_t rs_pid;
	port_t rs_port;
      };
#endif
    }; // end class Redis

    class MultiClientRedisLock {
    public:
      MultiClientRedisLock(const std::string& value, std::shared_ptr<Redis> db);
      explicit MultiClientRedisLock(const std::string& value, std::shared_ptr<Redis> db, unsigned int to);
      explicit MultiClientRedisLock(const std::string& value, std::shared_ptr<Redis> db, bool print_if_okay);
      MultiClientRedisLock(const std::string& value, std::shared_ptr<Redis> db, unsigned int to, bool print_if_okay);
      ~MultiClientRedisLock();
      bool locked();
      const static unsigned int DEFAULT_TIMEOUT;
    private:
      std::string val_;
      std::shared_ptr<Redis> db_;
      bool acquired_;
      bool pos_; // print "OK"/normal statuses? default = true
    };

    template <typename P>
    class RedisThriftSmartWriter : public ferrum::thrift::ThriftSmartWriter<P> {
    private:
      std::shared_ptr<ferrum::db::Redis> redis_;
      typedef apache::thrift::transport::TMemoryBuffer MemoryBuffer;
      boost::shared_ptr<MemoryBuffer> transport_;
    public:
      RedisThriftSmartWriter<P>(const std::string& name, std::shared_ptr<ferrum::db::Redis> redis_db);
      RedisThriftSmartWriter<P>(const std::string& name);
      RedisThriftSmartWriter<P>& operator=(const RedisThriftSmartWriter<P>&) = delete;
      RedisThriftSmartWriter<P>(const RedisThriftSmartWriter<P>&) = delete;
      void redis(std::shared_ptr<ferrum::db::Redis> redis_db);
      virtual ~RedisThriftSmartWriter();
      P* get(const std::string& suffix);
      virtual std::string form_id(const std::string& partial);
      virtual RedisThriftSmartWriter<P>* clone() const;
    protected:
      void _save();
    };
  } // ends namespace db

  template <typename P>
  class TarThriftSmartWriter : public ferrum::thrift::ThriftSmartWriter<P> {
  private:
    std::shared_ptr<ferrum::CompressedTar> archive_;
    typedef apache::thrift::transport::TMemoryBuffer MemoryBuffer;
    boost::shared_ptr<MemoryBuffer> transport_;
  public:
    TarThriftSmartWriter<P>(const std::string& name, std::shared_ptr<ferrum::CompressedTar> archive);
    TarThriftSmartWriter<P>(const std::string& name);
    TarThriftSmartWriter<P>& operator=(const TarThriftSmartWriter<P>&) = delete;
    TarThriftSmartWriter<P>(const TarThriftSmartWriter<P>&) = delete;
    virtual ~TarThriftSmartWriter();
    P* get(const std::string& suffix);
    virtual std::string form_id(const std::string& partial);
    virtual TarThriftSmartWriter<P>* clone() const;
    std::string safe_path(const std::string&);
  protected:
    void _save();
  };

} // ends namespace ferrum

#include "ferrum/data_util.tcc"

#endif
