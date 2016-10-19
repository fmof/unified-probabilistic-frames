namespace cpp minsky
namespace py minsky.bdoc

include "vocab.thrift"
include "labels.thrift"

enum WordAnnotation {
  UNSPECIFIED = 1000000000,
  ORTHOGRAPHIC = 0,
  LEMMA = 1
}

struct CountList {
  /*1: optional list<i32> items
  2: optional list<i32> icounts*/
   1: optional map<i32, i32> icounts
  /*3: optional list<double> dcounts*/
}

/**
 * A container for words. This can represent either a
 * dense, sequential view, or a sparse (BOW) view.
 */
struct WordsClause {
  /**
   * The following should only be set for BOW views
   */
  1: optional list<i32> words
  /**
   * The following should only be set for BOW views
   */
  2: optional CountList counts
  /*2: optional i32 which_word_vocab*/
  /*100: optional WordAnnotationLevel annot_level = WordAnnotationLevel.UNSPECIFIED*/
}

/**
 * Simple document
 */
struct SimpleDoc {
  1: optional string id
  2: optional list<WordsClause> sentences
  10: optional list<labels.Label> labels
}
