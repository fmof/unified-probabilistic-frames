# concrete-cpp/util/Makefile

#########################################

CONCRETE_CORE := $(realpath ../core)
CONCRETE_INCLUDE_DIR := $(CONCRETE_CORE)/include

########################################

# For generating the shared object
CONCRETE_LIB_DIR := $(CONCRETE_CORE)/lib
CONCRETE_SO := $(CONCRETE_LIB_DIR)/concrete$(SHLIB_SUFFIX)

########################################

CONCRETE_CXX_COMPILE := $(CXX) $(CXXFLAGS) -shared -fPIC

.PHONY: all clean base_env_set run_thrift so_env_set help test

##################################################
##################################################
########### GENERATED CONCRETE TARGETS ###########
##################################################
##################################################

BUILD_DIR := build
LIB_DIR := lib
INC_DIR := include

UTIL_NAME := concrete_util
NAMESPACE_NAME := $(UTIL_NAME)
LIB_NAME := lib$(UTIL_NAME)$(SHLIB_SUFFIX)
UTIL_SO_LIB := $(shell pwd)/$(LIB_DIR)
export CONCRETE_UTIL_SO := $(UTIL_SO_LIB)/$(LIB_NAME)
export CONCRETE_UTIL_STAT := $(UTIL_SO_LIB)/lib$(UTIL_NAME).a


SRC_H_FILES := $(wildcard src/*.h)
SRC_CPP_FILES := $(wildcard src/*.cc)

OBJ_FILES := $(addprefix $(BUILD_DIR)/,$(notdir $(SRC_CPP_FILES:.cc=.o)))
INC_H_FILES := $(subst src/,$(INC_DIR)/$(UTIL_NAME)/,$(SRC_H_FILES))

all: base_env_set $(INC_H_FILES) $(OBJ_FILES) $(CONCRETE_UTIL_SO) $(CONCRETE_UTIL_STAT)
	@echo "Done: please see $(GEN_CONCRETE)"

$(INC_DIR)/$(UTIL_NAME)/%.h: src/%.h
	cp $< $@


$(BUILD_DIR)/%.o: src/%.cc
	@if [ ! -d $(BUILD_DIR) ]; then mkdir $(BUILD_DIR); fi
	$(CXX) -shared -fPIC $(CXXFLAGS) -Isrc -c -o $@ $<

$(CONCRETE_UTIL_SO): $(OBJ_FILES)
	@if [ ! -d $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi 
	$(CXX) -shared -fPIC $(CXXFLAGS) -o $@ $^ -lthrift -lthriftz -lz -lboost_log -lboost_filesystem -lboost_system $(THRIFT_LIB_DIR)/libthrift$(SHLIB_SUFFIX)

$(CONCRETE_UTIL_STAT): $(OBJ_FILES)
	ar rcs $@ $(OBJ_FILES)

base_env_set:
	@if [ ! -d $(BUILD_DIR) ]; then mkdir $(BUILD_DIR); fi
	@if [ ! -d $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi
	@if [ ! -d $(INC_DIR) ]; then mkdir $(INC_DIR); fi
	@if [ ! -d $(INC_DIR)/$(UTIL_NAME) ]; then mkdir $(INC_DIR)/$(UTIL_NAME); fi

clean: clean-test
	rm -rf $(BUILD_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(INC_DIR)

clean-test:
	make -C test clean


##################################################
##################################################
################# TEST TARGETS ###################
##################################################
##################################################
test: $(INC_H_FILES) $(OBJ_FILES) $(CONCRETE_UTIL_SO)
	$(eval oldCXX = $(CXXFLAGS))
	$(eval CXXFLAGS += -isystem $(CONCRETE_INCLUDE_DIR) -I$(CONCRETE_INCLUDE_DIR) -L$(CONCRETE_LIB_DIR) -I$(realpath $(INC_DIR)) -L$(realpath $(LIB_DIR)) )
	make -C test all PASSED_CXX='$(CXXFLAGS)' UTIL_INCLUDE='$(realpath $(INC_DIR))' CORE_INCLUDE='$(CORE_INCLUDE)' CORE_LIB='$(CORE_LIB)' UTIL_LIB='$(realpath $(LIB_DIR))'
	$(eval CXXFLAGS = $(oldCXX))

.PHONY: install install-so install-headers
install-headers: $(INC_H_FILES)
	mkdir -p $(INSTALL_PREFIX_)/include/$(NAMESPACE_NAME)
	install -m 0644 $^ $(INSTALL_PREFIX_)/include/$(NAMESPACE_NAME)
install-so: $(CONCRETE_UTIL_SO) 
	mkdir -p $(INSTALL_PREFIX_)/lib
	install -m 0755 $< $(INSTALL_PREFIX_)/lib
install: install-headers install-so


help:
	@echo '--------------------------------------------------------------------------------'
	@echo 'concrete-cpp::util makefile'
	@echo '--------------------------------------------------------------------------------'


