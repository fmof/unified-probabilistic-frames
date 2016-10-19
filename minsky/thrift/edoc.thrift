namespace cpp minsky
namespace py minsky.edoc

include "vocab.thrift"

enum AnnotationLevel {
  UNSPECIFIED = 1000000000,
  SYNTAX = 0,
  SEMANTIC = 1
}

enum StructureType {
  UNSPECIFIED = 1000000000,
  PREDICATE = 0,
  RELATION = 1
}


struct RelationFiller {
  1: optional i32 word
  2: optional i32 position
}

struct PredArg {
  1: optional RelationFiller predicate
  3: optional i32 relation
  10: optional i32 which_word_vocab
  11: optional i32 which_relation_vocab
  100: optional AnnotationLevel annot_level = AnnotationLevel.UNSPECIFIED
}

struct Mention {
  1: optional string id
  2: optional list<PredArg> structures
  3: optional RelationFiller argument
  4: optional string location_id
}

struct Entity {
  1: optional string id
  2: optional list<Mention> mentions
}

/**
 * Entity-based documents
 */
struct EDoc {
  1: optional string id
  2: optional list<Entity> entities
}
