include ../Makefile.config

#########################################

#make sure CONC_SRC_DIR is set
GENERATED_TOP := $(shell pwd)
GEN_DIR = $(GENERATED_TOP)/gen-cpp
NAMESPACE_NAME = concrete
GEN_CONCRETE = $(GENERATED_TOP)/include/$(NAMESPACE_NAME)

BASE_GEN := $(subst ${CONC_SRC_DIR}/,${GEN_DIR}/,$(shell find ${CONC_SRC_DIR} -type f -name '*.thrift'))

# FILE_constants.{cpp,h}
# FILE_types.{cpp,h}
GEN_CONST_CPP := ${BASE_GEN:.thrift=_constants.cpp}
GEN_CONST_H := ${BASE_GEN:.thrift=_constants.h}
GEN_TYPES_CPP := ${BASE_GEN:.thrift=_types.cpp}
GEN_TYPES_H := ${BASE_GEN:.thrift=_types.h}
BASE_CONCRETE_CORE_VERSION_H_ = version.hpp
CONCRETE_CORE_VERSION_H_ = src/$(BASE_CONCRETE_CORE_VERSION_H_)
_CONCRETE_CORE_VERSION_CPP_ = $(GEN_DIR)/version.cpp
CONCRETE_CORE_VERSION_CPP_ = $(GEN_CONCRETE)/version.cpp
ALL_GEN_H = $(subst $(GEN_DIR)/,$(GEN_CONCRETE)/,$(GEN_CONST_H) $(GEN_TYPES_H) )
ALL_GEN_CPP = $(subst $(GEN_DIR)/,$(GEN_CONCRETE)/,$(GEN_CONST_CPP) $(GEN_TYPES_CPP))
ALL_GEN_OBJ_ = $(patsubst %.cpp,%.o,$(GEN_CONST_CPP) $(GEN_TYPES_CPP) $(_CONCRETE_CORE_VERSION_CPP_) )
ALL_GEN_CPP_H = $(ALL_GEN_H) $(ALL_GEN_CPP)

THRIFT_COMPILE = $(THRIFT) -gen cpp -I ${CONC_SRC_DIR} -o $(GENERATED_TOP)

########################################

# For generating the shared object
LIB_SHORTNAME := concrete
LIB_NAME := lib$(LIB_SHORTNAME)$(SHLIB_SUFFIX)
GEN_SO_LIB := $(GENERATED_TOP)/lib
export CONCRETE_SO := $(GEN_SO_LIB)/$(LIB_NAME)
export CONCRETE_STAT := $(GEN_SO_LIB)/lib$(LIB_SHORTNAME).a

########################################

CONCRETE_CXX_COMPILE := $(CXX) $(CXXFLAGS) -shared -fPIC

.PHONY: all clean base_env_set so_env_set help test
##################################################
##################################################
########### GENERATED CONCRETE TARGETS ###########
##################################################
##################################################
all: base_env_set $(CONCRETE_CORE_VERSION_CPP_) $(CONCRETE_SO) $(CONCRETE_STAT)
	@$(ECHO) "$(OK_STRING) Done: please see $(OK_COLOR)$(GEN_CONCRETE)$(NO_COLOR)"
	@$(ECHO) ""

base_env_set:
	@if [ ! -d $(GENERATED_TOP) ]; then mkdir $(GENERATED_TOP); fi
	@if [ ! -d $(GEN_CONCRETE) ]; then mkdir $(GEN_CONCRETE); fi
	@if [ ! -d $(GEN_CONCRETE)/concrete ]; then mkdir $(GEN_CONCRETE)/concrete; fi

so_env_set:
	@if [ ! -d $(GEN_SO_LIB) ]; then mkdir $(GEN_SO_LIB); fi


.CONC_VER = $(shell grep "Current version:" $(THRIFT_README) | tr -d '`' | cut -d: -f2- | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$$//')
.__version__file__.:
	@if [ -z "X$(.CONC_VER)" ]; then echo "Cannot determine version from Thrift README: Please verify that the readme exists (at $<)"; exit 1; else echo $(.CONC_VER) > $@; fi

$(GEN_SO_LIB):
	@if [ ! -d $(GEN_SO_LIB) ]; then mkdir $(GEN_SO_LIB); fi

$(GEN_CONCRETE)/%_constants.h: $(CONC_SRC_DIR)/%.thrift
	@echo "Compiling thrift $< to make $*_constants.h"
	$(THRIFT_COMPILE) $<
	cp $(GEN_DIR)/$*_constants.h $@
$(GEN_CONCRETE)/%_constants.cpp: $(CONC_SRC_DIR)/%.thrift
	@echo "Compiling thrift $< to make $*_constants.cpp"
	$(THRIFT_COMPILE) $<
	cp $(GEN_DIR)/$*_constants.cpp $@
$(GEN_CONCRETE)/%_types.h: $(CONC_SRC_DIR)/%.thrift
	@echo "Compiling thrift $< to make $*_types.h"
	$(THRIFT_COMPILE) $<
	cp $(GEN_DIR)/$*_types.h $@
$(GEN_CONCRETE)/%_types.cpp: $(CONC_SRC_DIR)/%.thrift
	@echo "Compiling thrift $< to make $*_types.cpp"
	$(THRIFT_COMPILE) $<
	cp $(GEN_DIR)/$*_types.cpp $@

.PHONY: generated-thrift
generated-thrift: $(ALL_GEN_CPP_H)

$(_CONCRETE_CORE_VERSION_CPP_): $(CONCRETE_CORE_VERSION_H_) .__version__file__.
	@rm -rf $@
	@touch $@
	@echo "#include \"$(BASE_CONCRETE_CORE_VERSION_H_)\"" >> $@
	@echo "" >> $@
	@echo "namespace concrete {" >> $@
	@echo "const std::string CONCRETE_VERSION=\"$(.CONC_VER)\";" >> $@
	@echo "}" >> $@
	@echo "" >> $@
	@rm .__version__file__.
	cp $< $(GEN_CONCRETE)/.
	cp $< $(GEN_DIR)/.

$(CONCRETE_CORE_VERSION_CPP_): $(_CONCRETE_CORE_VERSION_CPP_)
	cp $< $@

%.o: %.cpp
	$(CONCRETE_CXX_COMPILE) -o $@ -c $<

BOTH_LIB_DEPS = $(ALL_GEN_CPP_H) $(CONCRETE_CORE_VERSION_CPP_) $(ALL_GEN_OBJ_)

$(CONCRETE_SO): $(BOTH_LIB_DEPS) | $(GEN_SO_LIB)
	@$(ECHO) "$(STATUS_COLOR)Compiling generated thrift files to$(NO_COLOR) $(OK_COLOR)$(CONCRETE_SO)$(NO_COLOR)"
#	$(CONCRETE_CXX_COMPILE) -o $(CONCRETE_SO) $(subst $(GEN_DIR)/,$(GEN_CONCRETE)/,$(GEN_CONST_CPP) $(GEN_TYPES_CPP)) $(CONCRETE_CORE_VERSION_CPP_)
	$(CONCRETE_CXX_COMPILE) -o $(CONCRETE_SO) $(ALL_GEN_OBJ_)
	@$(ECHO) "$(OK_STRING) Made core SHARED library $(OK_COLOR)$@$(NO_COLOR)"

$(CONCRETE_STAT): $(BOTH_LIB_DEPS) | $(GEN_SO_LIB)
	ar rcs $@ $(ALL_GEN_OBJ_)
	@$(ECHO) "$(OK_STRING) Made core STATIC library $(OK_COLOR)$@$(NO_COLOR)"

clean: clean-test
	rm -rf $(GEN_DIR)
	rm -f $(subst ${GEN_DIR}/,${GEN_CONCRETE}/,$(ALL_GEN_CPP_H))
	rm -rf $(GEN_CONCRETE)
	rm -rf $(GEN_SO_LIB)
	rm -rf $(GEN_CONCRETE)
	rm -rf $(ALL_GEN_OBJ_)
	rm -rf $(CONCRETE_CORE_VERSION_CPP_)
	@rm -rf  .__version__file__.

clean-test:
	make -C test clean

clean-full: clean
	rm -rf $(GEN_CONCRETE)

##################################################
##################################################
################# TEST TARGETS ###################
##################################################
##################################################
GEN_CONC_INCLUDE = $(GEN_CONCRETE)/*.cpp $(GEN_CONCRETE)/*.h
test: 
	$(eval oldCXX = $(CXXFLAGS))
	$(eval CXXFLAGS += -isystem $(GENERATED_TOP)/include/ -I$(GENERATED_TOP)/include/ -L$(GEN_SO_LIB) )
	make -C test all PASSED_CXX='$(CXXFLAGS)'
	$(eval CXXFLAGS = $(oldCXX))

.PHONY: install install-so install-headers
install-headers: $(ALL_GEN_H) | $(GEN_CONCRETE)
	mkdir -p $(INSTALL_PREFIX_)/include/$(NAMESPACE_NAME)
	install -m 0644 $(ALL_GEN_H) $(INSTALL_PREFIX_)/include/$(NAMESPACE_NAME)
install-so: $(CONCRETE_SO) | $(GEN_SO_LIB)
	mkdir -p $(INSTALL_PREFIX_)/lib
	install -m 0755 $< $(INSTALL_PREFIX_)/lib
install: install-headers install-so


help:
	@echo $(ALL_GEN_OBJ_)
	@echo '--------------------------------------------------------------------------------'
	@echo 'concrete-cpp makefile'
	@echo ''
	@echo '--------------------------------------------------------------------------------'


