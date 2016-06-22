# concrete-cpp/core/test/Makefile

TEST_EXC_DIR = exc
TEST_BUILD_DIR = build
TEST_LIB_DIR = lib
UNIT_TESTS = $(TEST_EXC_DIR)/main_unittest

GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

CXXFLAGS = $(PASSED_CXX) -isystem $(GTEST_DIR)/include -I$(GTEST_DIR)/include -I$(GTEST_DIR)
CXXFLAGS += $(TEST_INCLUDE) $(TEST_LIB) $(TEST_FLAGS)
CXXFLAGS += -I$(CONCRETE_INCLUDE_DIR) -L$(CONCRETE_LIB_DIR)

AR=ar -rv

.PHONY: all setup_env

all: clean setup_env $(UNIT_TESTS)

setup_env:
	echo $(CONCRETE_CPP_PROJ_DIR)
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

$(TEST_BUILD_DIR)/comm_linking.o: src/comm_linking.cc 
	$(CXX) $(CXXFLAGS) -c src/comm_linking.cc && mv comm_linking.o $(TEST_BUILD_DIR)/.
$(TEST_EXC_DIR)/main_unittest: $(TEST_BUILD_DIR)/comm_linking.o $(TEST_LIB_DIR)/gtest.a $(CONCRETE_SO)
	@$(ECHO) "$(STATUS_COLOR)Trying to compile and assemble unit test $(OK_COLOR) $@$(NO_COLOR)"
	$(CXX) $(CXXFLAGS) -o $@ -lthrift $^ $(THRIFT_SO) $(CONCRETE_SO)
#	$(CXX) $(CXXFLAGS) -o $@ -lthrift $^ $(THRIFT_SO) ../lib/concrete.so
	$(call verify_exist,$@)
	@$(ECHO) "Running unit test $(STATUS_COLOR)$@$(NO_COLOR)"
	$(RUNNER) $@
	@$(ECHO) ""
