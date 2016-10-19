include ../Makefile.config

#default target
help:

BUILD_WHERE=$(BUILD_DIR)/build
INCLUDE_WHERE=$(BUILD_DIR)/include
LIB_WHERE=$(BUILD_DIR)/lib
EXC_WHERE=$(BUILD_DIR)/exec

ALL_BUILD_DIRS := $(sort $(BUILD_DIR) $(BUILD_WHERE) \
	$(INCLUDE_WHERE) $(LIB_WHERE) $(EXC_WHERE) \
	$(foreach build,$(BUILD_WHERE) $(INCLUDE_WHERE), \
		$(subst $(TOP_LEVEL_SRC)/,$(build)/,$(SRC_SUB_DIRS))))

################################################################

.PHONY: all core clean env_set help test util

##################################################
##################################################
######## GENERATED/BASIC MINSKY TARGETS ########
##################################################
##################################################

GENERATED_TOP := $(shell pwd)
GEN_TYPE := cpp
GEN_DIR = $(GENERATED_TOP)/gen-$(GEN_TYPE)
NAMESPACE_NAME = minsky
MINSKY_SRC_DIR = thrift
GEN_MINSKY = $(GENERATED_TOP)/include/$(NAMESPACE_NAME)

BASE_GEN := $(subst thrift/,${GEN_DIR}/,$(shell find thrift/ -type f -name '*.thrift'))

# FILE_constants.{cpp,h}
# FILE_types.{cpp,h}
GEN_CONST_CPP := ${BASE_GEN:.thrift=_constants.cpp}
GEN_CONST_H := ${BASE_GEN:.thrift=_constants.h}
GEN_TYPES_CPP := ${BASE_GEN:.thrift=_types.cpp}
GEN_TYPES_H := ${BASE_GEN:.thrift=_types.h}
ALL_GEN_H = $(subst $(GEN_DIR)/,$(GEN_MINSKY)/,$(GEN_CONST_H) $(GEN_TYPES_H) )
ALL_GEN_CPP = $(subst $(GEN_DIR)/,$(GEN_MINSKY)/,$(GEN_CONST_CPP) $(GEN_TYPES_CPP))
ALL_GEN_CPP_H = $(ALL_GEN_H) $(ALL_GEN_CPP)

THRIFT ?= thrift

THRIFT_COMPILE = $(THRIFT) -gen cpp -I ${MINSKY_SRC_DIR} -o $(GENERATED_TOP)

########################################

# For generating the shared object
LIB_SHORTNAME := minskythrift
LIB_NAME := lib$(LIB_SHORTNAME)$(MINSKY_LINKING_TYPE)
export MINSKY_SO ?= $(LIB_WHERE)/$(LIB_NAME)

########################################

MINSKY_CXXFLAGS ?= $(CXXFLAGS)
MINSKY_CXX_COMPILE := $(CXX) $(MINSKY_CXXFLAGS) -shared -fPIC

##################################################
##################################################
########### GENERATED MINSKY TARGETS ###########
##################################################
##################################################
core: base_env_set $(MINSKY_SO)
	@$(ECHO) "$(OK_STRING) Done: please see $(OK_COLOR)$(GEN_MINSKY)$(NO_COLOR)"
	@$(ECHO) ""

.PHONY: py
py:
	bash build_python.sh

base_env_set:
	@if [ ! -d $(GENERATED_TOP) ]; then mkdir -p $(GENERATED_TOP); fi
	@if [ ! -d $(GEN_MINSKY) ]; then mkdir -p $(GEN_MINSKY); fi

so_env_set:
	@if [ ! -d $(LIB_WHERE) ]; then mkdir $(LIB_WHERE); fi

$(LIB_WHERE):
	@if [ ! -d $(LIB_WHERE) ]; then mkdir $(LIB_WHERE); fi

$(GEN_MINSKY)/%_constants.h: $(MINSKY_SRC_DIR)/%.thrift
	@echo "Compiling thrift $< to make $*_constants.h"
	$(THRIFT_COMPILE) $<
	mv $(GEN_DIR)/$*_constants.h $@
$(GEN_MINSKY)/%_constants.cpp: $(MINSKY_SRC_DIR)/%.thrift
	@echo "Compiling thrift $< to make $*_constants.cpp"
	$(THRIFT_COMPILE) $<
	mv $(GEN_DIR)/$*_constants.cpp $@
$(GEN_MINSKY)/%_types.h: $(MINSKY_SRC_DIR)/%.thrift
	@echo "Compiling thrift $< to make $*_types.h"
	$(THRIFT_COMPILE) $<
	mv $(GEN_DIR)/$*_types.h $@
$(GEN_MINSKY)/%_types.cpp: $(MINSKY_SRC_DIR)/%.thrift
	@echo "Compiling thrift $< to make $*_types.cpp"
	$(THRIFT_COMPILE) $<
	mv $(GEN_DIR)/$*_types.cpp $@

.PHONY: generated-thrift
generated-thrift: $(ALL_GEN_CPP_H)

$(ALL_BUILD_DIRS):
	$(foreach d,$@, if [ ! -d $(d) ]; then mkdir -p $(d); fi)

$(MINSKY_LIB_DIR)/libminskythrift.so: $(ALL_GEN_CPP_H) | $(LIB_WHERE) $(GEN_MINSKY) $(ALL_BUILD_DIRS)
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
	@$(ECHO) "$(STATUS_COLOR)Compiling generated thrift files to$(NO_COLOR) $(OK_COLOR)$(MINSKY_SO)$(NO_COLOR)"
	$(MINSKY_CXX_COMPILE) -o $(MINSKY_SO) $(subst $(GEN_DIR)/,$(GEN_MINSKY)/,$(GEN_CONST_CPP) $(GEN_TYPES_CPP))
	@$(ECHO) "$(OK_STRING) Made core library $(OK_COLOR)$@$(NO_COLOR)"
	@if [ ! -d $(shell pwd)/lib ]; then mkdir -p $(shell pwd)/lib; fi
	if [ ! $(MINSKY_SO) -eq $(shell pwd)/lib/$(LIB_NAME) ]; then cp $(MINSKY_SO) $(shell pwd)/lib/$(LIB_NAME); fi

$(MINSKY_LIB_DIR)/libminskythrift.a: $(ALL_GEN_CPP_H) | $(LIB_WHERE) $(GEN_MINSKY) $(ALL_BUILD_DIRS)
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
	@$(ECHO) "$(STATUS_COLOR)Compiling generated thrift files to$(NO_COLOR) $(OK_COLOR)static$(NO_COLOR)"
	$(MINSKY_CXX_COMPILE) -o $(MINSKY_SO) $(subst $(GEN_DIR)/,$(GEN_MINSKY)/,$(GEN_CONST_CPP) $(GEN_TYPES_CPP))
	@$(ECHO) "$(OK_STRING) Made core library $(OK_COLOR)$@$(NO_COLOR)"
	if [ ! -d $(shell pwd)/lib ]; then mkdir -p $(shell pwd)/lib; fi
	if [ ! $(MINSKY_SO) -eq $(shell pwd)/lib/$(LIB_NAME) ]; then cp $(MINSKY_SO) $(shell pwd)/lib/$(LIB_NAME); fi

clean: clean-test
	rm -rf $(GEN_DIR)
	rm -f $(subst ${GEN_DIR}/,${GEN_MINSKY}/,$(ALL_GEN_CPP_H))
	rm -rf $(GEN_MINSKY)
	rm -rf $(LIB_WHERE)
	rm -rf $(GEN_MINSKY)


all: core
	@if [ ! -d $(GEN_MINSKY) ]; then mkdir -p $(GEN_MINSKY); fi
	cp --preserve=all $(GEN_DIR)/* $(GEN_MINSKY)
	@echo "Done: please see $(GEN_MINSKY)"

##################################################
##################################################
################ CLEAN TARGETS ###################
##################################################
##################################################

help:
#	@echo '$(ALL_GEN_CPP_H)'
	@echo '--------------------------------------------------------------------------------'
#	@echo 'Command-line options:'
#	@echo '  - PREFIX:   path to installation prefix directory'
#	@echo '--------------------------------------------------------------------------------'

test-core:
	make -C core test
	@$(ECHO) "$(OK_STRING) All core unit tests passed"
test-util: 
	make -C util test CORE_INCLUDE='$(realpath core/include)' CORE_LIB='$(realpath core/lib)'
	@$(ECHO) "$(OK_STRING) All util unit tests passed"
test: test-core test-util
	@$(ECHO) "$(OK_STRING) All unit tests passed"

.PHONY: verify-prefix install install-core install-util
verify-prefix:
	@if [ -z "$(PREFIX)" ]; then echo "$(ERROR_STRING) Please set the $(ERROR_COLOR)PREFIX$(NO_COLOR) variable."; exit 1; fi

install-core: verify-prefix
	make -C core install
	@$(ECHO) "$(OK_STRING) Core library installed"
install-util: verify-prefix
	make -C util install CORE_INCLUDE='$(realpath core/include)' CORE_LIB='$(realpath core/lib)'
	@$(ECHO) "$(OK_STRING) Util library installed"
install: install-core install-util
	@$(ECHO) "$(OK_STRING) All libraries installed"
