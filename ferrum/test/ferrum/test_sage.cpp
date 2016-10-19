#ifndef PUBLIC_UPF
#include "gtest/gtest.h"

#include "ferrum/logging.hpp"
#include "ferrum/sage_defs.hpp"

#include <cstdio>
#include <unordered_map>
#include <utility>
#include <vector>

typedef ferrum::SageTopic<std::vector<double> > Topic;

Topic create() {
  Topic topic(2, 1.0, ferrum::SageTopicRegularization::L2);
  std::shared_ptr<std::vector<double> > background(new std::vector<double>(2));
  background->operator[](0) = 5;
  background->operator[](1) = 1;
  topic.background(background);
  return topic;
}

TEST(SageTopic, create_function) {
  create();
}

TEST(SageTopic, create_copy_ctor) {
  Topic topic = create();
  // under g++ -g -O0, the following actually manifests
  // in the copy constructor
  Topic t2 = topic;
  std::vector<double>* bkg = topic.background();
  EXPECT_TRUE(bkg != NULL);
  std::vector<double>* bkg2 = t2.background();
  EXPECT_TRUE(bkg2 != NULL);
  EXPECT_TRUE(bkg == bkg2);
  ASSERT_EQ(2, bkg->size());
}
TEST(SageTopic, create_move_ctor) {
  Topic topic = create();
  // under g++ -g -O0, the following actually manifests
  // in the move constructor
  Topic t2 = std::move(topic);
  std::vector<double>* bkg = topic.background();
  EXPECT_TRUE(bkg == NULL);
  std::vector<double>* bkg2 = t2.background();
  ASSERT_TRUE(bkg2 != NULL);
  ASSERT_EQ(2, bkg2->size());
}
#endif
