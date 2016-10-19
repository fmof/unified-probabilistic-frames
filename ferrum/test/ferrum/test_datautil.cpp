#include "gtest/gtest.h"

#include "ferrum/concrete.hpp"
#include "ferrum/data_util.hpp"

#include <string>

TEST(Thrift, write_struct_str) {
  concrete::UUID uuid;
  uuid.__set_uuidString("my-id");
  std::string srep =
    ferrum::thrift::thrift_struct_to_string<ferrum::thrift::TCompactProtocol>(uuid);
  concrete::UUID uuid2;
  ferrum::thrift::thrift_struct_from_string<ferrum::thrift::TCompactProtocol>(srep, &uuid2);
  ASSERT_EQ(uuid.uuidString, uuid2.uuidString);
}

TEST(GZip, roundtrip_easy_str) {
  std::string input("Hello, world");
  std::string compressed(ferrum::compress::gzip::compress(input));
  ASSERT_EQ(32, compressed.size());
  std::string decompressed(ferrum::compress::gzip::decompress(compressed));
  ASSERT_EQ(input, decompressed);
}

TEST(GZip, roundtrip_comm) {
  concrete::Communication communication;
  concrete::util::concrete_io concrete_reader;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  concrete_reader.deserialize_compact<concrete::Communication>(&communication, name);
  ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
  std::string comm_str( ferrum::thrift::thrift_struct_to_string<ferrum::thrift::TCompactProtocol>(communication) );
  std::string gzip_comm_str( ferrum::compress::gzip::compress(comm_str) );
  std::string dgzip_comm_str( ferrum::compress::gzip::decompress(gzip_comm_str) );
  ASSERT_LE(gzip_comm_str.size(), dgzip_comm_str.size());
  concrete::Communication ncomm;
  ferrum::thrift::thrift_struct_from_string<ferrum::thrift::TCompactProtocol>(dgzip_comm_str, &ncomm);
  ASSERT_EQ(communication, ncomm);
}

TEST(TarThriftSmartWriter, safe_path) {
  ferrum::TarThriftSmartWriter<ferrum::thrift::TCompactProtocol> ttsw("my_name");
  EXPECT_EQ("my_path", ttsw.safe_path("my_path"));
  EXPECT_EQ("my_path", ttsw.safe_path("my path"));
  EXPECT_EQ("my__path", ttsw.safe_path("my :path"));
  EXPECT_EQ("_my_path__foo__bar_", ttsw.safe_path(":my_path::foo::bar!"));
}
