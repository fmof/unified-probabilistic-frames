#include "ferrum/util.hpp"
#include "ferrum/util.tcc"

#include "gtest/gtest.h"
#include "ferrum/logging.hpp"

#include <iostream>
#include <sstream>

#include <vector>

TEST(Util, sort_indices_ascending) {
  std::vector<int> vec;
  vec.push_back(4);
  vec.push_back(10);
  vec.push_back(0);
  auto s_indices = ferrum::sort_indices(vec, true);
  EXPECT_EQ(2, s_indices[0]);
  EXPECT_EQ(0, s_indices[1]);
  EXPECT_EQ(1, s_indices[2]);
}

TEST(Util, sort_indices_descending) {
  std::vector<int> vec;
  vec.push_back(4);
  vec.push_back(10);
  vec.push_back(0);
  auto s_indices = ferrum::sort_indices(vec, false);
  std::vector<size_t> expected = { 1, 0, 2};
  EXPECT_EQ(1, expected[0]);
  EXPECT_EQ(0, expected[1]);
  EXPECT_EQ(2, expected[2]);

  EXPECT_EQ(1, s_indices[0]);
  EXPECT_EQ(0, s_indices[1]);
  EXPECT_EQ(2, s_indices[2]);
}

TEST(Util, sum_1d_vecs) {
  std::vector<int> vec;
  vec.push_back(4);
  vec.push_back(10);
  vec.push_back(0);
  std::vector<int> v = ferrum::sum(vec,vec);
  ASSERT_EQ(3, v.size());
  EXPECT_EQ(8, v[0]);
  EXPECT_EQ(20, v[1]);
  EXPECT_EQ(0, v[2]);
}

TEST(Util, max_1d) {
  std::vector<int> vec;
  vec.push_back(4);
  vec.push_back(10);
  vec.push_back(0);
  int m = ferrum::max(vec);
  ASSERT_EQ(10, m);
}

TEST(Util, max_1d_all_negative) {
  std::vector<int> vec;
  vec.push_back(-4);
  vec.push_back(-10);
  vec.push_back(-1);
  int m = ferrum::max(vec);
  ASSERT_EQ(-1, m);
}

TEST(Util, max_2d) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 0, 3, 2});
  vec.push_back(std::vector<int> { -1, 192, -10});
  int m = ferrum::max<int>(vec);
  ASSERT_EQ(192, m);
}

TEST(Util, column_max) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 0, 3, 2});
  vec.push_back(std::vector<int> { -1, 192, -10});
  std::vector<int> c_maxes = ferrum::column_max<int>(vec);
  ASSERT_EQ(3, c_maxes.size());
  EXPECT_EQ(0, c_maxes[0]);
  EXPECT_EQ(192, c_maxes[1]);
  EXPECT_EQ(2, c_maxes[2]);
}
TEST(Util, row_max) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 0, 3, 2});
  vec.push_back(std::vector<int> { -1, 192, -10});
  std::vector<int> r_maxes = ferrum::row_max<int>(vec);
  ASSERT_EQ(2, r_maxes.size());
  EXPECT_EQ(3, r_maxes[0]);
  EXPECT_EQ(192, r_maxes[1]);
}

TEST(Util, histogram) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 0, 3, 2});
  vec.push_back(std::vector<int> { 6, 1, 2});
  vec.push_back(std::vector<int> { 4, 0, 0});
  std::vector< std::vector<int > > histogram = ferrum::histogram(vec);
  ASSERT_EQ(3, histogram.size());
  ASSERT_EQ(7, histogram[0].size());
  ASSERT_EQ(7, histogram[1].size());
  ASSERT_EQ(7, histogram[2].size());
  std::vector< std::vector<int> > expected;
  expected.push_back({1, 0, 0, 0, 1, 0, 1});
  expected.push_back({1, 1, 0, 1, 0, 0, 0});
  expected.push_back({1, 0, 2, 0, 0, 0, 0});
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 7; ++j) {
      EXPECT_EQ(expected[i][j], histogram[i][j]) << "(" << i << ", " << j << ") not equal";
    }
  }
}

TEST(Util, histogram1) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 3, 2});
  vec.push_back(std::vector<int> { 3, 2});
  vec.push_back(std::vector<int> { 2, 1});
  std::vector< std::vector<int > > histogram = ferrum::histogram(vec);
  ASSERT_EQ(2, histogram.size());
  ASSERT_EQ(4, histogram[0].size());
  ASSERT_EQ(4, histogram[1].size());
  std::vector< std::vector<int> > expected;
  expected.push_back({0,0,1,2});
  expected.push_back({0, 1, 2, 0});
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 4; ++j) {
      EXPECT_EQ(expected[i][j], histogram[i][j]) << "(" << i << ", " << j << ") not equal";
    }
  }
}

TEST(Util, sum_1d) {
  std::vector<int> vec;
  vec.push_back(4);
  vec.push_back(10);
  vec.push_back(0);
  int m = ferrum::sum(vec);
  ASSERT_EQ(14, m);
}

TEST(Util, sum_in_first) {
  std::vector<int> vec = {4, 10, 0};
  std::vector<int> other = {7, -3, 1};
  ferrum::sum_in_first(&vec, other);
  std::vector<int> expected = {11, 7, 1};
  for(size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(expected[i], vec[i]);
  }
}
TEST(Util, sum_in_first_map) {
  std::vector<int> vec = {4, 10, 0};
  std::map<int, int> other;
  other[1] = -3;
  ferrum::sum_in_first(&vec, other);
  std::vector<int> expected = {4, 7, 0};
  for(size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(expected[i], vec[i]);
  }
}
TEST(Util, sum_in_first_double) {
  std::vector<double> vec = {4, 10, 0};
  std::vector<double> other = {7, -3, 1};
  ferrum::sum_in_first(&vec, other);
  std::vector<double> expected = {11, 7, 1};
  for(size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(expected[i], vec[i]);
  }
}
TEST(Util, lc_in_first) {
  std::vector<int> vec = {4, 10, 0};
  std::vector<int> other = {7, -3, 1};
  int a = 5, b = 3;
  ferrum::linear_combination_in_first(&vec, other, a, b);
  std::vector<int> expected = {41, 41, 3};
  for(size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(expected[i], vec[i]);
  }
}

TEST(Util, lc_in_first_double) {
  std::vector<double> vec = {4, 10, 0};
  std::vector<double> other = {7, -3, 1};
  double a = 5, b = 3;
  ferrum::linear_combination_in_first(&vec, other, a, b);
  std::vector<double> expected = {41, 41, 3};
  for(size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(expected[i], vec[i]);
  }
}

TEST(Util, scalar_product_int_int_1d) {
  std::vector<int> vec;
  vec.push_back(4);
  vec.push_back(10);
  vec.push_back(0);
  const int scale = 5;
  std::vector<int> ret = ferrum::scalar_product(scale, vec);
  ASSERT_EQ(20, ret[0]);
  ASSERT_EQ(50, ret[1]);
  ASSERT_EQ(0,  ret[2]);
}

TEST(Util, scalar_product_int_double_1d) {
  std::vector<double> vec;
  vec.push_back(4.0);
  vec.push_back(10.0);
  vec.push_back(0.0);
  const int scale = 5;
  std::vector<double> ret = ferrum::scalar_product(scale, vec);
  ASSERT_EQ(20.0, ret[0]);
  ASSERT_EQ(50.0, ret[1]);
  ASSERT_EQ(0.0,  ret[2]);
}

TEST(Util, marginals) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 0, 3, 2});
  vec.push_back(std::vector<int> { -1, 192, -10});
  std::vector<int> marginals = ferrum::marginals(vec);
  ASSERT_EQ(2, marginals.size());
  ASSERT_EQ(5, marginals[0]);
  ASSERT_EQ(181, marginals[1]);
}

TEST(Util, sparse_histogram) {
  std::vector<int> marginals({5, 181, 182, 181});
  std::map<int, int> s_hist = ferrum::sparse_histogram(marginals);
  ASSERT_EQ(3, s_hist.size());
  EXPECT_EQ(1, s_hist[5]);
  EXPECT_EQ(2, s_hist[181]);
  EXPECT_EQ(1, s_hist[182]);
}

TEST(Util, marginal_histogram) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 0, 3, 2});
  vec.push_back(std::vector<int> { 6, 1, 2});
  vec.push_back(std::vector<int> { 4, 0, 0});
  std::vector< int > m_histogram = ferrum::marginal_histogram(vec);
  ASSERT_EQ(10, m_histogram.size());
  std::vector< int > expected({0, 0, 0, 0, 1, 1, 0, 0, 0, 1});
  for(int i = 0; i < 10; ++i) {
    EXPECT_EQ(expected[i], m_histogram[i]) << "(" << i << ") not equal";
  }
}

TEST(Util, marginal_histogram1) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 3, 2});
  vec.push_back(std::vector<int> { 3, 2});
  vec.push_back(std::vector<int> { 2, 1});
  std::vector<int > m_histogram = ferrum::marginal_histogram(vec);
  ASSERT_EQ(6, m_histogram.size());
  std::vector<int> expected({0,0,0,1,0,2});
  for(int i = 0; i < 6; ++i) {
    EXPECT_EQ(expected[i], m_histogram[i]) << "(" << i << ") not equal";
  }
}

TEST(Util, square_ref) {
  std::vector<int> vec = {5, 4};
  std::vector<int> expected = {25, 16};
  std::vector<int> res = ferrum::square(vec);
  for(int i = 0; i < 2; ++i) {
    ASSERT_NEAR(expected[i], res[i], 1E-6) << "(" << i << ") not equal";
  }
}
TEST(Util, cube_ref) {
  std::vector<int> vec = {5, 4};
  std::vector<int> expected = {125, 64};
  std::vector<int> res = ferrum::cube(vec);
  for(int i = 0; i < 2; ++i) {
    ASSERT_NEAR(expected[i], res[i], 1E-6) << "(" << i << ") not equal";
  }
}
TEST(Util, quartic_ref) {
  std::vector<int> vec = {5, 4};
  std::vector<int> expected = {625, 256};
  std::vector<int> res = ferrum::quartic(vec);
  for(int i = 0; i < 2; ++i) {
    ASSERT_NEAR(expected[i], res[i], 1E-6) << "(" << i << ") not equal";
  }
}

TEST(Util, slice) {
  std::vector< std::vector<int> > vec;
  vec.push_back(std::vector<int> { 0, 3, 2});
  vec.push_back(std::vector<int> { -1, 192, -10});
  std::vector<int> slice;
  for(size_t col = 0; col < 3; ++col) {
    slice = ferrum::column(vec, col);
    ASSERT_EQ(2, slice.size());
    ASSERT_EQ(vec[0][col], slice[0]);
    ASSERT_EQ(vec[1][col], slice[1]);
  }
}

TEST(Util, frobenius_norm) {
  std::vector<std::vector<int> > matrix = {{2, -1}, {-1, 2}};
  double expected = sqrt(10);
  double val = ferrum::frobenius_norm(matrix);
  ASSERT_NEAR(expected, val, 1E-5);
}

TEST(Util, smart_writer_to_dev_null) {
  ferrum::SmartWriter sw("/dev/null");
  std::stringstream buffer;
  std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
  std::ostream& writer = sw.get();
  writer << "This is a test\n";
  std::string text = buffer.str();
  ASSERT_EQ(0, text.size());
  std::cout.rdbuf(old);
}

TEST(Util, mat_vec_product1) {
  std::vector< std::vector<int> > mat;
  mat.push_back(std::vector<int> { 3, 2});
  mat.push_back(std::vector<int> { 9, -1});
  mat.push_back(std::vector<int> { 2, 1});
  std::vector< int > vec = {2, -1};
  std::vector< int > prod(3);
  ferrum::product(mat, vec, &prod);
  ASSERT_EQ(3, prod.size());
  std::vector<int> expected = { 4, 19, 3};
  for(size_t i = 0; i < 3; ++i) {
    EXPECT_EQ(expected[i], prod[i]) << "(" << i << ") not equal";
  }
}

TEST(Util, Eigen_mat_vec_product1) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat(3, 2);
  mat << 3, 2, 9, -1, 2, 1;
  std::vector< double > vec = {2, -1};
  std::vector< double > prod(3);
  ferrum::product(mat, vec, &prod);
  ASSERT_EQ(3, prod.size());
  std::vector<double> expected = { 4, 19, 3};
  for(size_t i = 0; i < 3; ++i) {
    EXPECT_EQ(expected[i], prod[i]) << "(" << i << ") not equal";
  }
}

TEST(Util, Eigen_mat_vec_product_aliasing) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat(2, 2);
  mat << 3, 2, -1, 2 ;
  std::vector< double > vec = {2, -1};
  ferrum::product(mat, vec, &vec);
  ASSERT_EQ(2, vec.size());
  std::vector<double> expected = { 4, -4 };
  for(size_t i = 0; i < 2; ++i) {
    EXPECT_EQ(expected[i], vec[i]) << "(" << i << ") not equal";
  }
}

TEST(Util, mat_vec_transpose_product1) {
  std::vector< std::vector<int> > mat;
  mat.push_back(std::vector<int> { 3, 2});
  mat.push_back(std::vector<int> { 9, -1});
  mat.push_back(std::vector<int> { 2, 1});
  std::vector< int > vec = {7, 2, -1};
  std::vector< int > prod(2);
  ferrum::transpose_product(mat, vec, &prod);
  ASSERT_EQ(2, prod.size());
  std::vector<int> expected = { 37, 11};
  for(size_t i = 0; i < 2; ++i) {
    EXPECT_EQ(expected[i], prod[i]) << "(" << i << ") not equal";
  }
}

TEST(Util, vectors_equal) {
  CAPTURE_STDERR;
  {
    std::vector<int> v1 = {1, 2, 6};
    std::vector<int> v2 = {1, 2, 6};
    ASSERT_TRUE(ferrum::vectors_equal(v1, v2));
  }
  {
    std::vector<int> v1 = {1, 6, 2};
    std::vector<int> v2 = {1, 2, 6};
    ASSERT_FALSE(ferrum::vectors_equal(v1, v2));
  }
  {
    std::vector<int> v1 = {1, 2, 6};
    std::vector<int> v2 = {1, 2, 4};
    ASSERT_FALSE(ferrum::vectors_equal(v1, v2));
  }
  {
    std::vector<int> v1 = {1, 2, 4, 6};
    std::vector<int> v2 = {1, 2, 6};
    ASSERT_FALSE(ferrum::vectors_equal(v1, v2));
  }
  {
    std::vector<int> v1 = {1, 2, 6};
    std::vector<int> v2 = {1, 2, 4, 6};
    ASSERT_FALSE(ferrum::vectors_equal(v1, v2));
  }
  {
    std::vector<int> v1 = {1, 2, 6};
    std::vector<int> v2 = {-3, -2, -10};
    ASSERT_FALSE(ferrum::vectors_equal(v1, v2));
  }
}

TEST(Util, maps_equal) {
  CAPTURE_STDERR;
  typedef std::unordered_map<int, int> map;
  {
    map m1 = { {1, 1}, {0, 2}, {6, 1} };
    map m2 = { {1, 1}, {0, 2}, {6, 1} };
    bool e = ferrum::maps_equal(m1, m2);
    ASSERT_TRUE( ( e ) );
  }
  {
    map m1 = { {1, 1}, {0, 2}, {6, 1} };
    map m2 = { {1, 5}, {0, 2}, {6, 1} };
    bool e = ferrum::maps_equal(m1, m2);
    ASSERT_FALSE( ( e ) );
  }
}

TEST(Util, dot) {
  std::vector<int> v1 = {1, 2, 6};
  std::vector<int> v2 = {1, 2, 4};
  ASSERT_EQ(29, ferrum::dot(v1, v2));
  ASSERT_EQ(0, ferrum::dot(std::vector<int>(), std::vector<int>()));
}

TEST(Util, median) {
  {
    std::vector<int> v1 = {1, 2, 6};
    ASSERT_EQ(2, ferrum::median(v1));
  }
  {
    std::vector<int> v2 = {9, 1, 2, 4};
    ASSERT_EQ(3, ferrum::median(v2));
  }
  {
    std::vector<int> v2 = {9, 1, 2, 4, 9};
    ASSERT_EQ(4, ferrum::median(v2));
  }
  {
    std::vector<int> v2 = {9, 1, 5, 4};
    ASSERT_EQ(4, ferrum::median(v2));
  }
  {
    std::vector<double> v2 = {9.0, 1.0, 5.0, 4.0};
    ASSERT_EQ(4.5, ferrum::median(v2));
  }
}

TEST(Util, ensure_min_ptr) {
  std::vector<double> vec = {5, 4, 2};
  std::vector<double> expected = {5, 4, 2.1};
  ferrum::ensure_min(2.1, vec.data(), vec.size());
  for(int i = 0; i < 3; ++i) {
    ASSERT_NEAR(expected[i], vec[i], 1E-6) << "(" << i << ") not equal";
  }
}

TEST(Util, sum_in_first_eigen_rvalue) {
  Eigen::VectorXd v(4);
  v << 1, 2, 3, 4;
  std::vector<double> x = {1, 2, 3, 4};
  ferrum::sum_in_first(&x, -1*v);
  for(auto y : x) {
    ASSERT_EQ(y, 0.0);
  }
}
