MODEL_SRC_DIRS = models
MODEL_CPP_FILES_ = $(foreach subdir,$(MODEL_SRC_DIRS),$(wildcard $(subdir)/*.cpp))
SRC += $(MODEL_CPP_FILES_)
REPLACED_SRC += $(subst $(MODEL_SRC_DIRS)/,$(SHARED_BUILD_WHERE)/models/,$(MODEL_CPP_FILES_))
UNIT_MODELS_NO_FULL = $(filter-out $(MODEL_SRC_DIRS),$(foreach file,$(subst $(MODEL_SRC_DIRS)/,,$(MODEL_CPP_FILES_:.cpp=)),models/$(file)))

.PHONY: $(UNIT_MODELS_NO_FULL)
#models/annotated_nyt_stats-obj = minsky.o data_util.o logging.o crtlda_concrete.o crtlda_pruner_minsky.o crtlda_minsky.o crtlda_sage_variational.o redis_corpus.o
#models/annotated_nyt_stats-obj = minsky.o data_util.o logging.o crtlda_pruner_minsky.o crtlda_minsky.o crtlda_sage_variational.o 
#models/annotated_nyt_stats-libnar-obj = $(foreach f,$(models/annotated_nyt_stats-obj),$(SHARED_BUILD_WHERE)/libnar/$(f))
# main_libnar/%:
# 	$(MAKE) -C ../src %.o
# %-libnar-prereqs : %_libnar-obj
