# concrete-cpp/util/test/Makefile

UTIL_INCLUDE ?= .
UTIL_INCLUDE_H := $(UTIL_INCLUDE)/concrete_util/*.h

TEST_EXC_DIR = exc
TEST_BUILD_DIR = build
TEST_LIB_DIR = lib
UNIT_TESTS = $(TEST_EXC_DIR)/main_unittest


GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

CXXFLAGS = $(PASSED_CXX) -isystem $(GTEST_DIR)/include -I$(GTEST_DIR)/include -I$(GTEST_DIR)
CXXFLAGS += $(TEST_INCLUDE) $(TEST_LIB) $(TEST_FLAGS)

AR=ar -rv

# Find the multi-threaded version of the Boost logging library
#
# If you install Boost from source, the multi-threaded library will be:
#   lib/libboost_log.dylib
# If you use Brew to install Boost on OS X, Brew will build Boost
# with the flags:
#   ./b2 --layout=tagged threading=multi,single
# which will change the default name of the multi-threaded version
# of the Boost logging library:
#   lib/libboost_log-mt.dylib  -  multi-threaded version
#   lib/libboost_log.dylib     -  single-threaded version
ifneq ($(wildcard $(BOOST_LIB_DIR)/libboost_log-mt$(SHLIB_SUFFIX)),)
	BOOST_LOG_LIB=boost_log-mt
else
	BOOST_LOG_LIB=boost_log
endif

.PHONY: all setup_env

all: clean setup_env $(UNIT_TESTS)

setup_env:
	@if [ ! -d $(TEST_EXC_DIR) ]; then mkdir $(TEST_EXC_DIR); fi
	@if [ ! -d $(TEST_BUILD_DIR) ]; then mkdir $(TEST_BUILD_DIR); fi
	@if [ ! -d $(TEST_LIB_DIR) ]; then mkdir $(TEST_LIB_DIR); fi

clean:
	rm -fr $(TESTS) $(TEST_BUILD_DIR) $(TEST_EXC_DIR) $(TEST_LIB_DIR)


# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
$(TEST_LIB_DIR)/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CXXFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc && \
	mv gtest-all.o $(TEST_LIB_DIR)/

$(TEST_LIB_DIR)/gtest.a : $(TEST_LIB_DIR)/gtest-all.o
	$(AR) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.

CPP_FILES := $(wildcard src/*.cc)
OBJ_FILES := $(subst src/,$(TEST_BUILD_DIR)/,$(CPP_FILES:.cc=.o))

$(TEST_BUILD_DIR)/%.o: src/%.cc $(GTEST_HEADERS)
	$(CXX) $(CXXFLAGS) -c src/$*.cc $(UTIL_INCLUDE_H) && mv $*.o $(TEST_BUILD_DIR)/.

$(TEST_EXC_DIR)/main_unittest: $(OBJ_FILES) $(TEST_LIB_DIR)/gtest.a 
	@$(ECHO) "$(STATUS_COLOR)Trying to compile and assemble unit test $(OK_COLOR) $@$(NO_COLOR)"
	$(CXX) $(CXXFLAGS) -o $@ -lthrift $^ -lthrift -lthriftz -lz -l$(BOOST_LOG_LIB) -lboost_system -lboost_filesystem $(CORE_LIB)/libconcrete$(SHLIB_SUFFIX) $(UTIL_LIB)/libconcrete_util$(SHLIB_SUFFIX) $(THRIFT_LIB_DIR)/libthrift$(SHLIB_SUFFIX)
	$(call verify_exist,$@)
	@$(ECHO) "Running unit test $(STATUS_COLOR)$@$(NO_COLOR)"
	$(RUNNER) $@
	@$(ECHO) ""
