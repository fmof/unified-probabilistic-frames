include Makefile.config

#default target
.PHONY: help
help:

################################################################

.PHONY: concrete ferrum clean env_set help test

.PHONY: minsky-cpp
minsky-cpp:
	$(MAKE) -C minsky all THRIFT="$(THRIFT_EXEC)" MINSKY_CXXFLAGS="$(CXXFLAGS)" MINSKY_SO="$(MINSKY_SO)"
$(MINSKY_SO):
	$(MAKE) -C minsky all THRIFT="$(THRIFT_EXEC)" MINSKY_CXXFLAGS="$(CXXFLAGS)" MINSKY_SO="$(MINSKY_SO)"

.PHONY: minsky-python
minsky-python:
	$(MAKE) -C minsky py THRIFT="$(THRIFT_EXEC)" MINSKY_CXXFLAGS="$(CXXFLAGS)" MINSKY_SO="$(MINSKY_SO)"

.PHONY: minsky
minsky: minsky-cpp $(MINSKY_SO) minsky-python

##################################################
##################################################
######## GENERATED/BASIC CONCRETE TARGETS ########
##################################################
##################################################

# each sub-dir will add to this
REPLACED_SRC :=
####################################
CONCRETE_CPP_CORE_MODULE_INCLUDE = $(CONCRETE_CORE_INCLUDE_DIR)/concrete
CONCRETE_CPP_CORE_MODULE_BUILD = $(SHARED_BUILD_WHERE)/concrete
CONCRETE_CPP_CORE_MODULE_THRIFT_WHERE = concrete-cpp/thrift/thrift
CONCRETE_CPP_CORE_MODULE_THRIFT_README = concrete-cpp/thrift/README.md
CONCRETE_CPP_CORE_THRIFT_EXEC = $(THRIFT_EXEC)
include concrete-cpp/core/Makefile.module
####################################
include concrete-cpp/util/Makefile.module
CONCRETE_CPP_UTIL_OBJ := $(patsubst %.cc,%.o, $(filter %.cc,$(CONCRETE_CPP_UTIL_REPLACED_SRC))) \
	$(patsubst %.cpp,%.o, $(filter %.cpp,$(CONCRETE_CPP_UTIL_REPLACED_SRC)))
CXXFLAGS += -I$(SHARED_INCLUDE_WHERE)
####################################
RELATIVE_PATH = ferrum
include ferrum/ferrum/Makefile.module
CXXFLAGS += -Iferrum
FERRUM_OBJ := $(patsubst %.cpp,%.o, $(filter %.cpp,$(FERRUM_REPLACED_SRC)))
####################################
include upf/Makefile
#CXXFLAGS += -Ihermes
UPF_OBJ := $(patsubst %.cpp,%.o, $(filter %.cpp,$(UPF_REPLACED_SRC)))
####################################
HPP_FILES_ = $(UPF_H_) $(FERRUM_H_)
TCC_FILES_ = $(UPF_TCC_) $(FERRUM_TCC_)
####################################
CXXFLAGS += -Imodels
include models/Makefile
####################################
RELATIVE_PATH = ferrum
include ferrum/test/Makefile.module
include test/Makefile
####################################
OBJ := 	$(CONCRETE_CORE_MODULE_GEN_OBJ_) \
	$(CONCRETE_CPP_UTIL_OBJ) \
	$(patsubst %.cpp,%.o, $(filter %.cpp,$(FERRUM_REPLACED_SRC) $(UPF_REPLACED_SRC)))

.PHONY: concrete_util_transfer
concrete_util_transfer: $(CONCRETE_CPP_UTIL_H_)
	@if [ ! -d $(SHARED_INCLUDE_WHERE)/concrete_util ]; then mkdir -p $(SHARED_INCLUDE_WHERE)/concrete_util; fi
	cp $^ $(SHARED_INCLUDE_WHERE)/concrete_util/

.PHONY: concrete
concrete: $(CONCRETE_CORE_MODULE_GEN_H) concrete_util_transfer $(CONCRETE_CORE_MODULE_GEN_OBJ_)

.PHONY: upf
upf: concrete minsky $(OBJ)
	@$(ECHO) "Done compiling"

.PHONY: minsky-path
minsky-path:
	export PYTHONPATH="$(shell pwd)/minsky/python:${PYTHONPATH}"

clean:
	rm -f $(DEPS)
	rm -f $(OBJ)
	rm -rf $(SHARED_BUILD_WHERE)
	rm -rf $(SHARED_EXC_WHERE)
	rm -rf $(SHARED_INCLUDE_WHERE)

.PHONY: full-clean test-clean
test-clean: clean
	rm -rf $(SHARED_TEST_WHERE)
full-clean: clean test-clean
	rm -rf $(SHARED_LIB_WHERE)


GTEST_DIR := gtest-1.7.0
GTEST_INCLUDE := -isystem $(GTEST_DIR)/include -I$(GTEST_DIR)/include -I$(GTEST_DIR)
GTEST_LIB := -L$(GTEST_DIR)/lib
GTEST_LIB_LINKING_FLAGS := -lpthread
GTEST_FLAGS += $(GTEST_INCLUDE) $(GTEST_LIB) -pthread $(GTEST_LIB_LINKING_FLAGS) -DGTEST_USE_OWN_TR1_TUPLE=0 -DUPF_UNIT_TEST -DFERRUM_UNIT_TEST
GTEST_INCLUDE += -I$(SHARED_INCLUDE_WHERE)


DEPS = $(OBJ:.o=.d)
-include $(DEPS)
$(SHARED_BUILD_WHERE)/test/%.d : test/%.cpp
	@if [ ! -d `dirname $@` ]; then mkdir -p `dirname $@`; fi
	@scripts/autodepend.sh \
            "$(CXX) -MG -MP -MM $(CXX_SHARED_FLAGS) -I$(GTEST_DIR)/include -I$(GTEST_DIR) $(CXXFLAGS)" \
            `dirname $*.cpp` \
            "$(SHARED_BUILD_WHERE)/test/" \
             $< > $@

$(SHARED_BUILD_WHERE)/%.d : %.cpp
	@if [ ! -d `dirname $@` ]; then mkdir -p `dirname $@`; fi
	@scripts/autodepend.sh \
           "$(CXX) -MG -MP -MM $(CXX_SHARED_FLAGS) $(CXXFLAGS)" \
           `dirname $*.cpp` \
           "$(SHARED_BUILD_WHERE)" \
            $< > $@

### HACK!!!
$(SHARED_BUILD_WHERE)/concrete_util/%.d : concrete-cpp/util/src/%.cc
	@if [ ! -d `dirname $@` ]; then mkdir -p `dirname $@`; fi
	@scripts/autodepend.sh \
           "$(CXX) -MG -MP -MM $(CXX_SHARED_FLAGS) $(CXXFLAGS) -Iconcrete-cpp/util/src" \
           concrete-cpp/util/src \
           "$(SHARED_BUILD_WHERE)" \
            $< > $@

$(SHARED_BUILD_WHERE)/test/%.o: test/%.cpp $(SHARED_BUILD_WHERE)/test/%.d 
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
	$(CXX) $(CXX_SHARED_FLAGS) $(CXXFLAGS) $(GTEST_FLAGS) -c $< -o $@


$(SHARED_BUILD_WHERE)/%.o: $(SHARED_BUILD_WHERE)/%.d %.cpp
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
	$(CXX) $(CXX_SHARED_FLAGS) $(CXXFLAGS) -c $*.cpp -o $@

### HACK!!!
$(SHARED_BUILD_WHERE)/concrete_util/%.o: $(SHARED_BUILD_WHERE)/concrete_util/%.d concrete-cpp/util/src/%.cc
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
	$(CXX) $(CXX_SHARED_FLAGS) $(CXXFLAGS) -Iconcrete-cpp/util/src -c concrete-cpp/util/src/$*.cc -o $@

#$(HPP_FILES_) $(TCC_FILES_) 
$(UNIT_MODELS_NO_FULL) : % : concrete $(OBJ) $(SHARED_BUILD_WHERE)/%.o
	@if [ ! -d $(SHARED_EXC_WHERE)/$(@D) ]; then mkdir -p $(SHARED_EXC_WHERE)/$(@D); fi
	$(CXX) $(CXXFLAGS_L) -Iferrum -g -o $(SHARED_EXC_WHERE)/$@ \
$(OBJ) $(SHARED_BUILD_WHERE)/$*.o $(MINSKY_LINKING_FLAGS) $(CXX_LINKING_FLAGS)

gtest_objects = $(SHARED_LIB_WHERE)/gtest_main-$(BUILD).o $(SHARED_LIB_WHERE)/gtest-all-$(BUILD).o

$(SHARED_TEST_WHERE)/% : concrete $(OBJ) $(SHARED_BUILD_WHERE)/test/%.o $(gtest_objects)
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
 # unfortunately, some libraries must be repeated (hence the multiple -l{thrift,gsl,gslcblas,m}
	$(CXX) $(CXXFLAGS_L) $(GTEST_LIB) -g -o $@ $(SHARED_BUILD_WHERE)/test/$*.o $(OBJ) $(gtest_objects) -lgsl -lgslcblas -lm $(GTEST_LIB_LINKING_FLAGS) $(MINSKY_LINKING_FLAGS) $(CXX_LINKING_FLAGS)

$(UNIT_TESTS_NO_FULL) : % : $(SHARED_TEST_WHERE)/%
	$(call verify_exist,$(SHARED_TEST_WHERE)/$*)
	@$(ECHO) "Running unit test $(STATUS_COLOR)$<$(NO_COLOR)"
	cd test && $(RUNNER) $(SHARED_TEST_WHERE)/$@
	@$(ECHO) ""

test: $(FERRUM_UNIT_TESTS_NO_FULL) $(UNIT_TESTS_NO_FULL)


#########################################################################################
#########################################################################################


GTEST_HEADERS_ = $(GTEST_DIR)/include/gtest/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS_)

AR=ar -rv
# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
$(SHARED_LIB_WHERE)/gtest-all-$(BUILD).o : $(GTEST_SRCS_)
	@$(ECHO) "$(STATUS_COLOR)Trying to compile and assemble gtest-all $(OK_COLOR) $(BUILD) $(STATUS_COLOR) version $(NO_COLOR)" 
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
	$(CXX) $(GTEST_FLAGS) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc -o $@ $(CXX_LINKING_FLAGS)

$(SHARED_LIB_WHERE)/gtest-$(BUILD).a : $(SHARED_LIB_WHERE)/gtest-all-$(BUILD).o
	@$(ECHO) "$(STATUS_COLOR)Creating static gtest-all $(OK_COLOR) $(BUILD) $(STATUS_COLOR) version $(NO_COLOR)" 
	$(AR) $@ $^

$(SHARED_LIB_WHERE)/gtest_main-$(BUILD).o : $(GTEST_SRCS_)
	@$(ECHO) "$(STATUS_COLOR)Trying to compile and assemble gtest-main $(OK_COLOR) $(BUILD) $(STATUS_COLOR) version $(NO_COLOR)" 
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
	$(CXX) $(GTEST_FLAGS) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc -o $@ $(CXX_LINKING_FLAGS)

$(SHARED_LIB_WHERE)/gtest_main-$(BUILD).a : $(SHARED_LIB_WHERE)/gtest-all-$(BUILD).o $(SHARED_LIB_WHERE)/gtest_main-$(BUILD).o
	@$(ECHO) "$(STATUS_COLOR)Creating static gtest-main $(OK_COLOR) $(BUILD) $(STATUS_COLOR) version $(NO_COLOR)" 
	$(AR) $@ $^

##################################################
##################################################
################# HELP TARGETS ###################
##################################################
##################################################

help:
	@echo 'upf makefile:'
	@echo '    Known targets:'
	@echo '     - upf:         Compile the basic upf (and dependent) objects.'
	@echo '     - test:	    Run all upf & ferrum tests'
	@echo '     - list-tests:  Find which tests can be run individually. For a test, e.g., test_foo, run "make test_foo"'
	@$(ECHO) $(foreach m,$(FERRUM_UNIT_TESTS_NO_FULL) $(UNIT_TESTS_NO_FULL), "\t\t* $(OK_COLOR)$(m)$(NO_COLOR)\n")
	@echo '     - list-models'
	@$(ECHO) $(foreach m,$(UNIT_MODELS_NO_FULL),"\t\t* $(OK_COLOR)$(m)$(NO_COLOR)\n")
	@echo '☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤☤'
