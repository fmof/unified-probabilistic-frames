namespace cpp minsky.residual
namespace py minsky.residual.rgs

include "frames.thrift"
include "vocab.thrift"
include "model_metadata.thrift"

struct SituationSlot {
  1: optional frames.Frame role
  2: optional string name
  //3: optional i32 id
}

struct SituationSlotFrame {
  1: optional frames.Distribution usage
  2: optional list<SituationSlot> slots
  10: optional string name
}

struct SituationTemplate {
  1: optional frames.Frame predicate_frame
  2: optional frames.Frame slot
  3: optional string name
  10: optional SituationSlotFrame unique_slots
  //4: optional i32 id
}

struct ResidualGlobalSlots {
  1: required list<SituationTemplate> thematics
  2: required list<SituationSlot> slots
  100: required list<vocab.Vocab> vocabularies

  200: required list<double> predicate_background
  201: required list<double> relation_background

  300: optional list<double> predicate_hyper
  301: optional list<double> relation_hyper
  302: optional list<double> slot_hyper

  1000: optional model_metadata.Metadata metadata
}

struct Semantics {
  1: optional list<SituationTemplate> preds
  2: optional list<SituationSlot> roles
}

struct ResidualUniqueSlots {
  1: required list<SituationTemplate> thematics
  2: optional Semantics semantics
  100: required list<vocab.Vocab> vocabularies

  200: optional list<double> predicate_background
  201: optional list<double> relation_background

  300: optional list<double> predicate_hyper
  301: optional list<double> relation_hyper
  302: optional list<double> slot_hyper
  303: optional list<double> sem_predicate_hyper
  304: optional list<double> sem_relation_hyper

  1000: optional model_metadata.Metadata metadata
}

