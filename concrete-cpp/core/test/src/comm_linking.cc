#include <gtest/gtest.h>

#include <boost/shared_ptr.hpp>

#include "concrete/communication_types.h"
#include "concrete/version.hpp"

#include "thrift/protocol/TBinaryProtocol.h"
#include "thrift/protocol/TProtocol.h"
#include "thrift/protocol/TProtocolException.h"
#include "thrift/protocol/TCompactProtocol.h"
#include "thrift/transport/TTransportUtils.h"
#include "thrift/transport/TFileTransport.h"
#include "thrift/transport/TFDTransport.h"

#include <iostream>
#include <fcntl.h>

double square(const double a) {
  double b = a*a;
  if(b != b) { // nan check
    return -1.0;
  }else{
    return b;
  }
}


TEST(Version, correct) {
  ASSERT_GT(concrete::CONCRETE_VERSION.size(), 0);
}

// Note that the following two methods do repeat code,
// but since this is _just_ to test that existing communications
// can be read in, and since there's no dependence on util,
// which _does_ have more reusable code, having this repeat
// code is okay.
void open_comm_tbinary(concrete::Communication *communication,
		       const char *file_name) {
  typedef apache::thrift::transport::TFDTransport TFDTransport;
  typedef apache::thrift::transport::TBufferedTransport TBufferedTransport;
  typedef apache::thrift::protocol::TBinaryProtocol TBinaryProtocol;

  int fd = open(file_name, O_RDONLY);

  boost::shared_ptr<TFDTransport> innerTransport(new TFDTransport(fd));
  boost::shared_ptr<TBufferedTransport> transport(new TBufferedTransport(innerTransport));
  boost::shared_ptr<TBinaryProtocol> protocol(new TBinaryProtocol(transport));
  transport->open();

  try {
    communication->read(protocol.get());
  } catch(const apache::thrift::protocol::TProtocolException &tpe) {
    std::cout << "An exception occurred. Exception: " << tpe.what() << std::endl;
    throw tpe;
  }
  transport->close();
  close(fd);
}
void open_comm_tcompact(concrete::Communication *communication,
			const char *file_name) {
  typedef apache::thrift::transport::TFDTransport TFDTransport;
  typedef apache::thrift::transport::TBufferedTransport TBufferedTransport;
  typedef apache::thrift::protocol::TCompactProtocol TCompactProtocol;

  int fd = open(file_name, O_RDONLY);

  boost::shared_ptr<TFDTransport> innerTransport(new TFDTransport(fd));
  boost::shared_ptr<TBufferedTransport> transport(new TBufferedTransport(innerTransport));
  boost::shared_ptr<TCompactProtocol> protocol(new TCompactProtocol(transport));
  transport->open();

  try {
    communication->read(protocol.get());
  } catch(const apache::thrift::protocol::TProtocolException &tpe) {
    std::cout << "An exception occurred. Exception: " << tpe.what() << std::endl;
    throw tpe;
  }
  transport->close();
  close(fd);
}
 
TEST(SquareTest, PositiveNos) { 
  ASSERT_EQ(36, square(6.0));
  ASSERT_EQ(0, square(0.0));
}

TEST(Communication, create) {
  concrete::Communication communication;
  communication.__set_id("my communication");
  ASSERT_EQ("my communication", communication.id);
}

TEST(Communication, readTBinaryProtocol) {
  concrete::Communication communication;
  const char *name = "resources/AFP_ENG_19940531.0390.tbinary.concrete";
  open_comm_tbinary(&communication, name);
  ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
}

TEST(Communication, readTCompactProtocol) {
  concrete::Communication communication;
  const char *name = "resources/AFP_ENG_19940531.0390.tcompact.concrete";
  open_comm_tcompact(&communication, name);
  ASSERT_EQ("AFP_ENG_19940531.0390", communication.id);
}
 
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
