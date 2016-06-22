
CONFIG_FILE := Makefile.config
include $(CONFIG_FILE)

#default target
help:

# If we want to debug the generated C++, then 
# explicitly set BUILD=debug
BUILD ?= release

# We're always going to compile to C++ 2011 standard
export RUNNER ?= $(SHELL)
CXX_BUILD_FLAGS_BASE = -std=c++11
CXX_BUILD_FLAGS_debug   = -g -O0
CXX_BUILD_FLAGS_release = -O3
export CXXFLAGS += $(CXX_BUILD_FLAGS_BASE) $(CXX_BUILD_FLAGS_$(BUILD))
ifeq ($(BUILD), "debug")
	RUNNER = gdb
else
        RUNNER += -c
endif

MAKE = make 

################################################################

.PHONY: concrete core clean env_set help test util

##################################################
##################################################
######## GENERATED/BASIC CONCRETE TARGETS ########
##################################################
##################################################
core:
	$(MAKE) -C core 

util:
	$(MAKE) -C util

concrete: core util
	@if [ ! -d $(GEN_CONCRETE) ]; then mkdir $(GEN_CONCRETE); fi
	@cp $(GEN_DIR)/* $(GEN_CONCRETE)
	@echo "Done: please see $(GEN_CONCRETE)"

##################################################
##################################################
################ CLEAN TARGETS ###################
##################################################
##################################################

clean:
	$(MAKE) -C core clean
	$(MAKE) -C util clean

clean-core:
	$(MAKE) -C core clean

clean-test:
	$(MAKE) -C test clean

clean-full: clean
	rm -rf $(GEN_CONCRETE)

##################################################
##################################################
################# TEST TARGETS ###################
##################################################
##################################################
help:
	@echo '--------------------------------------------------------------------------------'
	@echo 'concrete-cpp makefile:'
	@echo '  Targets:'
	@echo '  - core:         Compile the core concrete files'
	@echo '  - test-core:    Test the core concrete files'
	@echo '  - util:         Compile the utility library'
	@echo '  - test-util:    Test the utility library'
	@echo '  - all:          Same as running "make core util"'
	@echo '  - test:         Same as running "make test-core test-util"'
	@echo '  - install-core: Install the core files'
	@echo '  - install-util: Install the util files'
	@echo '  - install:      Install both the core and util files'
	@echo 'Command-line options:'
	@echo '  - PREFIX:   path to installation prefix directory'
	@echo '--------------------------------------------------------------------------------'

all: core util
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
