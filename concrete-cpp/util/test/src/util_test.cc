#include <gtest/gtest.h>

#include "concrete/communication_types.h"

#include "concrete_util/uuid_util.h"
#include "concrete_util/io.h"

#include <iostream>

#include <vector>

TEST(UUIDString, create) {
  concrete::util::uuid_factory uuid_maker;
  uuid_maker.get_uuid();
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

TEST(Communication, readCompressedCompact) {
  concrete::util::concrete_io concrete_reader;
  std::vector<concrete::Communication> comm_list = concrete_reader.deserialize_compressed_communications<concrete::util::TCompactProtocol>("resources/AFP_ENG_19940531.0390.tcompact.concrete.gz");
  ASSERT_EQ(2, comm_list.size());
  for(concrete::Communication communication : comm_list) {
    ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
  }
}

TEST(GZipCommunicationSequence, readCompressedCompact) {
  int num_comms = 0;
  concrete::util::GZipCommunicationSequence concrete_reader("resources/AFP_ENG_19940531.0390.tcompact.concrete.gz");
  for(concrete_reader.begin(); concrete_reader.keep_reading(); ++concrete_reader) {
    concrete::Communication communication = *concrete_reader;
    ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
    ++num_comms;
  }
  ASSERT_EQ(2, num_comms);
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

TEST(Communication, readTCompactProtocolDirectory) {
  int num_comms = 0;
  concrete::util::DirectoryCommunicationSequence concrete_reader("resources/AFP_ENG_19940531.0390.tcompact_dir");
  for(concrete_reader.begin(); concrete_reader.keep_reading(); ++concrete_reader) {
    concrete::Communication communication = *concrete_reader;
    ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
    ++num_comms;
  }
  ASSERT_EQ(2, num_comms);
}

TEST(Communication, readTCompactProtocolDirectoryFromFactoryMethod) {
  int num_comms = 0;
  concrete::util::CommunicationSequence *concrete_reader;
  concrete::util::get_communication_sequence("resources/AFP_ENG_19940531.0390.tcompact_dir",
					   concrete_reader);
  for(concrete_reader->begin(); concrete_reader->keep_reading(); concrete_reader->operator++()) {
    concrete::Communication communication = *(*concrete_reader);
    ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
    ++num_comms;
  }
  ASSERT_EQ(2, num_comms);
  delete concrete_reader;
}

 
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
