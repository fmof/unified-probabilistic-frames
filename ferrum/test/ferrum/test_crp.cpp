#include "gtest/gtest.h"

//#include "concrete/communication_types.h"
//#include "ferrum/concrete.hpp"
#include "ferrum/crp.hpp"
#include "ferrum/logging.hpp"

#include <unordered_map>

TEST(OccupancyHistogram, create) {
  crp::OccupancyHistogram oh;
}

TEST(OccupancyHistogram, createTable) {
  crp::OccupancyHistogram oh;
  oh.create_table();
  ASSERT_EQ(1, oh.tables());
  auto ohh = oh.histogram();
  ASSERT_EQ(1, ohh.size());
  ASSERT_TRUE(ohh.find(1) != ohh.end());
  ASSERT_EQ(1, oh[1]);
}

TEST(OccupancyHistogram, createTwoTables) {
  crp::OccupancyHistogram oh;
  oh.create_table();
  oh.create_table();
  ASSERT_EQ(2, oh.tables());
  auto ohh = oh.histogram();
  ASSERT_EQ(1, ohh.size());
  ASSERT_TRUE(ohh.find(1) != ohh.end());
  ASSERT_EQ(2, oh[1]);
}

TEST(OccupancyHistogram, accessInnerHist) {
  crp::OccupancyHistogram oh;
  oh.create_table();
  oh.create_table();
  ASSERT_EQ(2, oh.tables());
  oh[1] = 4;
  auto ohh = oh.histogram();
  ASSERT_EQ(1, ohh.size());
  ASSERT_TRUE(ohh.find(1) != ohh.end());
  ASSERT_EQ(4, oh[1]);
  ASSERT_EQ(4, ohh[1]);
}

TEST(CRP, create_conc1_discount0) {
  crp::CRP<> base_crp(1.0, 0.0);
}

TEST(CRP, seating) {
  crp::CRP<> base_crp(1.0, 0.0);
  base_crp.seat(1, 1.0/3.0);
  base_crp.seat(2, 1.0/3.0);
  base_crp.seat(1, 1.0/3.0);
  base_crp.seat(3, 1.0/3.0);
}

TEST(CRP, seat_unseat) {
  crp::CRP<> base_crp(1.0, 0.0);
  base_crp.seat(1, 1.0/3.0);
  base_crp.unseat(1);
}
