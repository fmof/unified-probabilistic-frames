namespace cpp minsky
namespace py minsky.frames

struct Weights {
  1: optional list<double> normalized
  2: optional list<double> unnormalized
  3: optional list<double> residual
  4: optional list<double> biases
}

struct Distribution {
  1: required Weights weights
  2: optional list<i32> items
  3: optional i32 vocab_idx = -1
  4: optional i32 support_size
  5: optional double sparsity_threshold = 0.0
  6: optional bool sparse = false
}

struct Frame {
  1: optional Distribution distr
  2: optional string name
}

struct FrameCollection {
  1: optional list<Frame> frames
  2: optional list<double> weights
}

