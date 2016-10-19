UPF_SRC_DIR = upf
UPF_SRC_ = $(shell find $(UPF_SRC_DIR) -type f -name '*.cpp') 
SRC += $(UPF_SRC_) $(SHARED_BUILD_WHERE)/upf/version.cpp
UPF_REPLACED_SRC += $(subst upf/,$(SHARED_BUILD_WHERE)/upf/,$(UPF_SRC_)) \
   $(SHARED_BUILD_WHERE)/upf/version.cpp
UPF_H_ = $(shell find $(UPF_SRC_DIR) -type f -name '*.hpp')
SRC += $(UPF_H_)
UPF_TCC_ = $(shell find $(UPF_SRC_DIR) -type f -name '*.tcc')
SRC += $(UPF_TCC_)

### SPECIAL CASE the version file
$(SHARED_BUILD_WHERE)/upf/version.cpp:
	@if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
	bash scripts/get_git_hash.sh $@ upf UPF "upf/version.hpp"
