namespace cpp minsky
namespace py minsky.vocab

struct Vocab {
  1: optional list<string> words
  2: optional string oov
  3: optional i32 oov_index = -1
}

struct VocabCollection {
  1: optional list< Vocab > vocabs
}
