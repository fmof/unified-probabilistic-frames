#include "gtest/gtest.h"

#include "ferrum/minsky.hpp"

#include <iostream>

TEST(MinskyMention, find_level) {
  minsky::PredArg syn;
  syn.annot_level = minsky::AnnotationLevel::SYNTAX;
  minsky::PredArg sem;
  sem.annot_level = minsky::AnnotationLevel::SEMANTIC;
  // 1
  {
    minsky::Mention m;
    m.structures.push_back(syn);
    m.structures.push_back(sem);
    ASSERT_EQ(0, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(1, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
  {
    minsky::Mention m;
    m.structures.push_back(sem);
    m.structures.push_back(syn);
    ASSERT_EQ(1, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(0, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
  {
    minsky::Mention m;
    m.structures.push_back(syn);
    ASSERT_EQ(0, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(-1, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
  {
    minsky::Mention m;
    m.structures.push_back(sem);
    ASSERT_EQ(-1, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(0, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
  {
    minsky::Mention m;
    ASSERT_EQ(-1, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(-1, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
  {
    minsky::Mention m;
    m.structures.push_back(sem);
    m.structures.push_back(sem);
    m.structures.push_back(sem);
    ASSERT_EQ(-1, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(0, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
  {
    minsky::Mention m;
    m.structures.push_back(sem);
    m.structures.push_back(sem);
    m.structures.push_back(syn);
    ASSERT_EQ(2, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(0, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
  {
    minsky::Mention m;
    m.structures.push_back(minsky::PredArg());
    m.structures.push_back(sem);
    m.structures.push_back(syn);
    ASSERT_EQ(2, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(1, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
  {
    minsky::Mention m;
    m.structures.push_back(minsky::PredArg());
    ASSERT_EQ(-1, minsky::find_level(m, minsky::AnnotationLevel::SYNTAX, false));
    ASSERT_EQ(-1, minsky::find_level(m, minsky::AnnotationLevel::SEMANTIC, false));
  }
}
