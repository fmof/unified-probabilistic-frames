TEST_SRC_DIRS = test/upf
TEST_CPP_FILES_ = $(foreach subdir,$(TEST_SRC_DIRS),$(wildcard $(subdir)/*.cpp))
SRC += $(TEST_CPP_FILES_)
REPLACED_SRC += $(subst $(TEST_SRC_DIRS)/,$(SHARED_BUILD_WHERE)/$(TEST_SRC_DIRS)/,$(TEST_CPP_FILES_))
UNIT_TESTS_NO_FULL = $(filter-out $(ALL_TEST_BUILD_DIRS),$(subst test/,,$(TEST_CPP_FILES_:.cpp=)))
