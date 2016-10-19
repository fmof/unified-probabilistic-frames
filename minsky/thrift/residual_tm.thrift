namespace cpp minsky.residual
namespace py minsky.residual.tm

include "frames.thrift"
include "vocab.thrift"
include "model_metadata.thrift"

struct Topic {
  1: optional frames.Frame frame
  2: optional string name
  //3: optional i32 id
}

struct ResidualTopicModel {
  1: required list<Topic> topics
  100: required list<vocab.Vocab> vocabularies

  200: required list<double> background

  300: optional list<double> word_hyper

  1000: optional model_metadata.Metadata metadata
}

