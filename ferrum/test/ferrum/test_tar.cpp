#include "gtest/gtest.h"

#include "ferrum/crtlda_defs.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/tar.hpp"

#include <boost/filesystem.hpp>
#include <cstdio>
#include <unordered_map>
#include <utility>
#include <vector>


TEST(Tar, list) {
  ferrum::CompressedTar ct("resources/cat_jumped.tar.gz");
  ct.read();
  int i = 0;
  for(ferrum::CompressedTar::iterator it = ct.begin();
      it != ct.end();
      ++it) {
    //std::string actual(archive_entry_pathname( it() ));
    std::string actual(ct.name(*it));
    switch(i) {
    case 0:
      ASSERT_EQ("cat_jumped.1.txt", actual);
      break;
    case 1:
      ASSERT_EQ("cat_jumped.2.txt", actual);
      break;
    default:
      ASSERT_TRUE(false);
    }
    const void* buffer = 0;
    size_t s;
    off_t o;
    int r = ct(&buffer, &s, &o);
    ASSERT_EQ(r, ARCHIVE_OK);
    std::string read_text(static_cast<const char*>(buffer), s);
    ASSERT_EQ("The cat jumped.\n", read_text);
    ++i;
  }
}
TEST(Tar, list_forrange) {
  ferrum::CompressedTar ct("resources/cat_jumped.tar.gz");
  ct.read();
  int i = 0;
  for(const auto& it : ct) {
    //std::string actual(archive_entry_pathname( it() ));
    std::string actual(ct.name(it));
    switch(i) {
    case 0:
      ASSERT_EQ("cat_jumped.1.txt", actual);
      break;
    case 1:
      ASSERT_EQ("cat_jumped.2.txt", actual);
      break;
    default:
      ASSERT_TRUE(false);
    }
    const void* buffer = 0;
    size_t s;
    off_t o;
    int r = ct(&buffer, &s, &o);
    ASSERT_EQ(r, ARCHIVE_OK);
    std::string read_text(static_cast<const char*>(buffer), s);
    ASSERT_EQ("The cat jumped.\n", read_text);
    ++i;
  }
}


TEST(Tar, list2) {
  ferrum::CompressedTar ct("resources/animals.tar.gz");
  ct.read();
  int i = 0;
  for(ferrum::CompressedTar::iterator it = ct.begin();
      it != ct.end();
      ++it) {
    std::string actual(ct.name(*it));
    const void* buffer = 0;
    size_t s;
    off_t o;
    int r = ct(&buffer, &s, &o);
    ASSERT_EQ(r, ARCHIVE_OK);
    std::string read_text(static_cast<const char*>(buffer), s);
    switch(i) {
    case 0:
      ASSERT_EQ("cat_jumped.1.txt", actual);
      ASSERT_EQ("The cat jumped.\n", read_text);
      break;
    case 1:
      ASSERT_EQ("dog_ran.1.txt", actual);
      ASSERT_EQ("A dog ran fast.\n", read_text);
      break;
    default:
      ASSERT_TRUE(false);
    }
    ++i;
  }
}

TEST(Tar, get_archived_corpus_compact1) {
  ferrum::Vocabulary<std::string> gov_voc;
  ferrum::Vocabulary<std::string> rel_voc;
  int num_comms = 0;
  std::string archive_name("resources/AFP_ENG_19940531.0390.tcompact.2copies.tar");
  const size_t expected_size = 49454;

  struct archive* arc_;
  arc_ = archive_read_new();
  struct archive_entry* tar_entry = NULL;
  archive_read_support_filter_none(arc_);
  archive_read_support_format_tar(arc_);
  ASSERT_EQ(ARCHIVE_OK, archive_read_open_filename(arc_, archive_name.c_str(), 10240));
  void* buffer = 0;
  size_t size = 0;
  int r;
  while((r = archive_read_next_header(arc_, &tar_entry)) == ARCHIVE_OK ) {
    std::string entry_name( archive_entry_pathname(tar_entry) );
    size = archive_entry_size(tar_entry);
    EXPECT_EQ(expected_size, size) << "expected size from TAR ENTRY failed on comm " << (num_comms+1);
    buffer = malloc(size);
    EXPECT_EQ(size, archive_read_data(arc_, buffer, size));
    free(buffer);
    ++num_comms;
  }

  EXPECT_EQ(ARCHIVE_EOF, r);
  ASSERT_EQ(2, num_comms);
  ASSERT_EQ(ARCHIVE_OK, archive_read_free(arc_));
}

TEST(Tar, write) {
  std::string fname("resources/test_cat_jumped.WRITE..tar.gz");
  if(boost::filesystem::exists(fname)) {
    boost::filesystem::remove(fname);
  }
  {
    ferrum::CompressedTar ct(fname);
    ct.write();
    std::string cat("The cat jumped high.\n");
    ASSERT_EQ(ARCHIVE_OK, ct.write_data("cat_jumped_memory.1.txt", &(cat[0]), cat.size()));
  }
  {
    ferrum::CompressedTar ctr(fname);
    ctr.read();
    int i = 0;
    for(ferrum::CompressedTar::iterator it = ctr.begin();
	it != ctr.end();
	++it) {
      std::string actual(ctr.name(*it));
      switch(i) {
      case 0:
	ASSERT_EQ("cat_jumped_memory.1.txt", actual);
	break;
      default:
	ASSERT_TRUE(false);
      }
      const void* buffer = 0;
      size_t s;
      off_t o;
      int r = ctr(&buffer, &s, &o);
      ASSERT_EQ(r, ARCHIVE_OK);
      std::string read_text(static_cast<const char*>(buffer), s);
      ASSERT_EQ("The cat jumped high.\n", read_text);
      ++i;
    }
  }
  if(boost::filesystem::exists(fname)) {
    boost::filesystem::remove(fname);
  } else {
    ASSERT_TRUE(false) << fname << " does not exist (to remove)";
  }
}
