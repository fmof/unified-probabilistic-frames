namespace cpp minsky
namespace py minsky

struct Label {
  1: optional i32 label
  10: optional double weight
}

struct ClassificationWeights {
  1: required list<double> weights
  2: optional i32 label
}
