#include "ferrum/minsky.hpp" // this needs to occur here
#include "ferrum/data_util.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/crtlda_defs.hpp"
#include "ferrum/crtlda_concrete.hpp"
#include "ferrum/crtlda_pruner_minsky.hpp"
#include "ferrum/crtlda_minsky.hpp"
#include "ferrum/crtlda_sampling.hpp"
#include "ferrum/redis_corpus.hpp"

#include "concrete_util/uuid_util.h"

#include "upf/upf_sampling.hpp"
#include "upf/upf_sampling.tcc"

#include <gsl/gsl_rng.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define CHECK_CORPUS_PTR(c) do {			\
    if( (c) == NULL ) {					\
      ERROR << "There was an error getting the corpus";	\
      exit(1);						\
    }							\
  } while(0)

// free if not null
#define FREE_IF_NNULL(c) do {			\
    if( (c) != NULL ) {				\
      delete (c);				\
    }						\
  } while(0)

void init_logging(char* program_name) {
#ifdef LOG_AS_BOOST
  boost::log::core::get()->set_filter
    (
     boost::log::trivial::severity >= boost::log::trivial::info
     );
#elif defined(LOG_AS_COUT)
#else
  //  google::InitGoogleLogging(program_name);
#endif
}

struct OptionClosureWrapper {
  int *num_templates;
  int *num_slots;
  int *num_epochs;
  ferrum::Toolnames* tools;
  int *num_docs_for_init;
};

struct OptionNames {
  const char* HELP = "help";
  const char* TRAIN_LOC = "train";
  const char* TRAIN_NAME = "train-name";
  const char* DEV_LOC = "dev";
  const char* DEV_NAME = "dev-name";
  const char* TRAIN_FROM_DB = "train-from-db";
  const char* DEV_FROM_DB = "dev-from-db";
  const char* SHUFFLE_TRAIN = "shuffle-train";
  const char* SHUFFLE_DEV = "shuffle-dev";
  const char* HOST = "host";
  const char* PORT = "port";
  const char* NO_RUN = "no-run";
  const char* UNIFY_VOCAB = "unify-vocab";
  const char* VOCAB_NAME = "vocab";
  const char* PRINT_AS_LOADING = "print-when-loading";
  ///////////////////////////////
  const char* NUM_TOPICS = "topics";
  const char* DESCENT = "descent";

  //////////////////////////////
  const char* PRUNER_TYPE = "pruner-type";
  const char* COMPUTE_BACKGROUND = "compute-background";
  const char* UNIFORM_BACKGROUND  = "uniform-background";
  const char* NUM_DOCS_INIT = "num-docs-init";

  const char* INFERENCER_SERIALIZATION = "inferencer-serialization";
  const char* SERIALIZED_INFERENCER = "serialized-inferencer";
};


int get_options(int n_args, char** args,
		po::variables_map* vm,
		ferrum::SamplingStrategy* strategy,
		ferrum::SymmetricHyperparams* hyperparams,
		OptionClosureWrapper* ocw,
		const OptionNames& names) {
  po::options_description desc("Allowed options");

  desc.add_options()
    ("help", "produce help message")
    ("train", po::value< std::string >(), 
     "input training location, either a tar.gz file on disk, or the name of a Redis list")
    ("train-name",
     po::value< std::string >()->default_value("train_corpus"),
     "name for the training corpus")
    ("train-number",
     po::value< int >()->default_value(-1),
     "only load N documents. For N < 0, load all. Default: -1")
    ("shuffle-train",
     po::bool_switch()->default_value(false),
     "shuffle the training corpus (in memory). default: false")
    ("dev", po::value< std::string >(), "input test path")
    ("dev-name",
     po::value< std::string >()->default_value("dev_corpus"),
     "name for the dev corpus")
    ("train-from-db", po::bool_switch()->default_value(false),
     "load the training set from an existing database")
    ("dev-from-db", po::bool_switch()->default_value(false),
     "load the dev set from an existing database")
    ("host", po::value< std::string >(), "redis host name")
    ("port", po::value< unsigned short >(), "redis port")
    ("no-run", po::bool_switch()->default_value(false),
     "do not actually run")
    ("num-dependency-hops", po::value<unsigned int>()->default_value(2),
     "number of dependency hops to make (default: 2)")
    ("dependency-name", po::value<std::string>(&(ocw->tools->dep_parse_tool)),
     "dependency parse tool name")
    ("unify-vocab", po::bool_switch()->default_value(false),
     "unify the vocabularies created in this program with a global one stored on the redis server")
    ("unified-doc-string", po::value<std::string>()->default_value(ferrum::REDIS_CORPUS_DOC_FIELD),
     "doc_str field name for unified vocabularies, if --unify-vocab")
    ("gov-vocab", po::value<std::string>()->default_value("gov_vocab"),
     "gov vocab name to load from the redis server")
    ("rel-vocab", po::value<std::string>()->default_value("rel_vocab"),
     "rel vocab name to load from the redis server")
    //////////////////////////
    ("templates", po::value<int>(ocw->num_templates)->default_value(10), 
     "number of templates to set")
    ("slots", po::value<int>(ocw->num_slots)->default_value(5), 
     "number of slots per template to use")

    ("hyper-template-usage", po::value<double>(&(hyperparams->h_theta)), 
     "template usage hyperparameter")
    ("hyper-slot-usage", po::value<double>(&(hyperparams->h_slot)), 
     "slot usage hyperparameter")
    ("hyper-verb", po::value<double>(&(hyperparams->h_gov)), 
     "verb hyperparameter")
    ("hyper-arc", po::value<double>(&(hyperparams->h_gov)), 
     "relation hyperparameter")

    // ("uniform-verb-background", po::value< double >(),
    //  "if provided, set all verb background values to this value")
    // ("uniform-arc-background", po::value< double >(),
    //  "if provided, set all arc background values to this value")

    ("pruner-type", 
     po::value< ferrum::MinskyPrunerEnum >()->default_value(ferrum::MinskyPrunerEnum::SituationGoverned),
     "which mention pruner to use; valid choices are a subset of ferrum::MinskyPrunerEnum. (choices: {VerbGoverned (== 1), *SituationGoverned (== 2),}")
    ("num-docs-init",
     po::value< int >(ocw->num_docs_for_init)->default_value(-1),
     "the number of documents to choose for subset initialization; default (-1) == use the Eisenstein et al. (2011) initialization technique, 10000 * #docs / #corpus tokens")

    // ("init-topics-by", 
    //  po::value< isage::wtm::TopicInitializerChoice >(&init_topics_how_)->default_value(isage::wtm::TopicInitializerChoice::SUBSET),
    //  "how to initialize the topic parameters; valid choices are given by isage::wtm::TopicInitializerChoice. (default: 1 == isage::wtm::TopicInitializerChoice::SUBSET)")
    // ("regularization", 
    //  po::value< isage::wtm::SageTopicRegularization >(&regularize_topics_how_)->default_value(isage::wtm::SageTopicRegularization::L2),
    //  "how to regularize the topics in MAP estimation; valid choices are given by isage::wtm::SageTopicRegularization. (default: 1 == isage::wtm::SageTopicRegularization::L2)")
    ("train-epochs", po::value<int>(ocw->num_epochs)->default_value(5),
     "Number of epochs to run")

    ("compute-elbo", po::bool_switch()->default_value(true),
     "compute ELBO during learning/inference (default: true)")
    ("no-compute-elbo", po::bool_switch()->default_value(false),
     "do not compute ELBO during learning/inference (default: false)")
    ("elbo-as-break", po::bool_switch()->default_value(true),
     "compute ELBO during learning/inference and use it for early stopping (default: true)")
    ("no-elbo-as-break", po::bool_switch()->default_value(false),
     "do not use ELBO during learning/inference as early stopping (default: false)")


    ("sparsity-threshold", po::value<double>()->default_value(0.0), "sparsity threshold")
    ("thrift-sparsity", po::bool_switch()->default_value(false), "store sparse SAGE distributions in serialized thrift")
    ("optimization-sparsity", po::bool_switch()->default_value(false), "enforce sparsity during optimization")
    ("num-samples", po::value<int>(&(strategy->num_iterations))->default_value(1000),
     "number of Gibbs sampling steps to take (default: 1000)")
    ("burnin", po::value<int>(&(strategy->burn_in))->default_value(100),
     "number of burn-in samples (default: 100)")
    // ("update-model-interval", po::value<int>(&(strategy->update_model_every))->default_value(5), "update the model every [some] number of EM steps (default: 5)")
    // ("never-update-model", po::bool_switch(&(strategy->never_update_model))->default_value(false),
    //  "never transfer parameters back to the model. This should be set for space-saving only")
    ("print-templates-every", po::value<int>()->default_value(25), "print templates every [some] number of sampling steps (default: 25)")
    ("print-usage-every", po::value<int>()->default_value(25), "print template usage every [some] number of sampling steps (default: 25)")
    // ("label-every", po::value<int>(&(strategy->label_every))->default_value(5), "label SitationMentions every [some] number of EM steps (default: 5)")
    ("label-host", po::value< std::string >(), "redis host name")
    ("label-port", po::value< unsigned short >(), "redis port")
    // ("top-k", po::value<int>(&(strategy->print_topics_k))->default_value(10), "number of words per topic to print (default: 10)")
    // ("em-verbosity", po::value<int>(&(strategy->em_verbosity))->default_value(1),
    //  "how verbose should EM output be (default: 1; higher == more verbose)")
    ////////////////////////////////
    (names.INFERENCER_SERIALIZATION, po::value<std::string>(), "filename to serialize inference state to")
    ("serialized-inferencer", po::value<std::string>(), "filename to READ serialized inference state from")
    ///////////////////////////////
    ("infer-on-test", po::value<bool>()->default_value(false),
     "run inference on the test set (default: true)")
    ("compute-dev-ll", po::bool_switch()->default_value(false),
     "compute log-likelihood on the heldout set")
    ("compute-dev-label", po::value<std::string>(),
     "label the heldout data")
    // ("compute-coherence", po::value<std::vector<std::string> >()->multitoken(),
    //  "compute topic coherence on all the space-separated distributions. Known values are gov_obs, rel_obs, gov_kind, rel_kind, template, slot")
    // // ("compute-coherence", po::bool_switch()->default_value(false),
    // //  "compute topic coherence. When true, both --coherence-file and --coherence-M must be set. (default: false)")
    // ("coherence-file", po::value<std::string>(),
    //  "file to write out coherence results")
    // ("coherence-M", po::value<int>(), "take top M for each distribution while computing coherences")
    // ("templates-to-file", po::value<std::string>(), "if given, print corpus (corpora) to file, otherwise do not print the corpora. these are read in from --serialized-inferencer")
    // //("lemma-mapper-delimiter", po::value<std::string>(), "Delimiter for lemma mapper file")
    ;

  po::store(po::parse_command_line(n_args, args, desc), *vm);
  if (vm->count("help")) {
    ERROR << desc << "\n";
    return 1;
  }
  po::notify(*vm);    
  return 0;
}

template <typename Vocab, typename Corpus, typename CDR>
bool read_dev
(
 const po::variables_map& vm,
 Vocab& gov_vocab,
 Vocab& rel_vocab,
 Corpus** corpus,
 const ferrum::Toolnames& tools,
 std::shared_ptr<ferrum::db::Redis> redis,
 int N,
 CDR corp_disk_reader
) {
  bool do_dev = false;
  if(vm.count("dev")) {
    if(vm["dev-from-db"].as<bool>()) {
      ERROR << "Cannot set both --dev and --dev-from-db";
      return 1;
    }
    //gov_kind_vocab.allow_new_words(false);
    //rel_kind_vocab.allow_new_words(false);
    ferrum::db::Address addr_h("localhost", 4532);
    FREE_IF_NNULL(*corpus);
    *corpus = new Corpus(vm["dev-name"].as<std::string>(), redis);
    corp_disk_reader
      (
       vm["dev"].as<std::string>(),
       N,
       vm["num-dependency-hops"].as<unsigned int>(),
       tools,
       gov_vocab,
       rel_vocab,
       *corpus
       );
    CHECK_CORPUS_PTR( *corpus );
    (*corpus)->multithreaded( redis->multithreaded() );
    do_dev = true;
  } else if(vm["dev-from-db"].as<bool>()) {
    INFO << "Going to read DEV from existing database";
    *corpus = new Corpus(vm["dev-name"].as<std::string>(), redis);
    CHECK_CORPUS_PTR( *corpus );
    (*corpus)->multithreaded( redis->multithreaded() );
    size_t num_loaded = (*corpus)->load_from_key(vm["dev-name"].as<std::string>());
    INFO << "Loaded " << num_loaded << " documents from dev corpus " << vm["dev-name"].as<std::string>();
    do_dev = true;
  }
  return do_dev;
}

template <typename Inferencer,
	  typename Vocab,
	  typename OutTar>
void save_state(const Inferencer& inf,
		const Vocab& gov_vocab, const Vocab& rel_vocab,
		OutTar& output_tar,
		const std::string& ser_name,
		int epoch) {
  minsky::residual::ResidualUniqueSlots rgs =
    inf.create_minsky_view(gov_vocab, rel_vocab);
  output_tar =
    std::shared_ptr<ferrum::CompressedTar>
    (new ferrum::CompressedTar
     (
      ser_name + 
      ".iteration" + std::to_string(epoch) + ".tar.gz",
      ferrum::CompressedTarOperation::WRITE
      )
     );
  std::string thrift_str
    (
     ferrum::thrift::thrift_struct_to_string<ferrum::thrift::TCompactProtocol>(rgs)
     );
  output_tar->write_data
    (
     "upf_syntax_only.epoch" + std::to_string(epoch) + ".thrift",
     &(thrift_str[0]),
     thrift_str.size()
     );
}

int main(int n_args, char** args) {
  init_logging(args[0]);

  int num_templates;
  int num_slots;
  int num_epochs_;
  int num_docs_for_init = -1;

  ferrum::Toolnames tools;
  OptionClosureWrapper ocw;
  ocw.num_templates = &num_templates;
  ocw.num_slots = &num_slots;
  ocw.num_epochs = &num_epochs_;
  ocw.tools = &tools;
  ocw.num_docs_for_init = &num_docs_for_init;

  ferrum::SampleEveryIter strategy(1000, 100);
  ferrum::SymmetricHyperparams hyperparams;

  po::variables_map vm;
  OptionNames names;
  int option_status = get_options(n_args, args, &vm, &strategy, &hyperparams, &ocw, names);
  if(option_status) {
    return option_status;
  }

  typedef std::string string;
  typedef ferrum::Vocabulary<string> SVocab;
  typedef minsky::EDoc Doc;
  typedef ferrum::CollapsedGibbsDMC Inferencer;

  // ferrum::StopWordList *gov_lemma_stops = NULL;
  // if(vm.count("gov-lemma-stoplist")) {
  //   gov_lemma_stops = new ferrum::StopWordList(vm["gov-lemma-stoplist"].as<std::string>());
  // }

  void (*corp_disk_reader)(const std::string&,  //path
			   int, // number to load
			   unsigned int,  // num dependency hops
			   const ferrum::Toolnames&,
			   ferrum::Vocabulary<std::string>&,
			   ferrum::Vocabulary<std::string>&,
			   ferrum::RedisCorpus<Doc>*) = NULL;
  switch(vm["pruner-type"].as< ferrum::MinskyPrunerEnum >()) {
  case ferrum::MinskyPrunerEnum::VerbGoverned:
    corp_disk_reader = &ferrum::set_archived_corpus<ferrum::RedisCorpus<Doc>,
						    ferrum::MinskyVerbGovernedPruner, // the pruner type
						    concrete::util::TCompactProtocol // how to read the Communications
						    >;
    break;
  case ferrum::MinskyPrunerEnum::SituationGoverned:
    corp_disk_reader = &ferrum::set_archived_corpus<ferrum::RedisCorpus<Doc>,
						    ferrum::MinskySituationGovernedPruner, // the pruner type
						    concrete::util::TCompactProtocol // how to read the Communications
						    >;
    break;
  case ferrum::MinskyPrunerEnum::UNSET:
  default:
    ERROR << "Unknown pruner type " << vm["pruner-type"].as< ferrum::MinskyPrunerEnum >();
    return 2;
  }
  assert(corp_disk_reader != NULL);

  if(! vm.count("host") || !vm.count("port") ) {
    ERROR << "Both --host and --port must be set";
    return 1;
  }
  std::string host(vm["host"].as<std::string>());
  unsigned short port(vm["port"].as<unsigned short>());

  ferrum::db::Address addr(host, port);
  std::shared_ptr<ferrum::db::Redis> main_redis_connection(new ferrum::db::Redis(addr));
  std::shared_ptr<ferrum::db::Redis> label_redis_connection;
  if(vm.count("label-port") == 0) {
    if(vm.count("label-host") == 0) {
      label_redis_connection = main_redis_connection;
    } else {
      ERROR << "Cannot have --label-host but not --label-port";
      return 1;
    }
  } else {
    if(vm.count("label-host") != 0) {
      label_redis_connection = std::shared_ptr<ferrum::db::Redis>
	(
	 new ferrum::db::Redis
	 (
	  ferrum::db::Address(vm["label-host"].as<std::string>(),
			       vm["label-port"].as<unsigned short>())
	  )
	 );
    } else {
      ERROR << "Cannot have --label-port but not --label-host";
      return 1;
    }
  }

  if(vm.count("train") || vm["train-from-db"].as<bool>() ) {
    // SVocab gov_kind_vocab("__GOV_KIND_OOV__");
    // SVocab rel_kind_vocab("__REL_KIND_OOV__");
    SVocab gov_vocab("__GOV_OOV__");
    SVocab rel_vocab("__REL_OOV__");

    ferrum::RedisCorpus<Doc>* corpus =
      new ferrum::RedisCorpus<Doc>(vm["train-name"].as<std::string>(),
				   main_redis_connection);
    const bool unify_train_vocab = vm["unify-vocab"].as<bool>();
    // if unifying, then change the docstr
    if(unify_train_vocab) {
      concrete::util::uuid_factory uf;
      corpus->doc_str( "temp::" + corpus->doc_str() + "::" +
		       vm["train"].as<std::string>() + "::" +
		       uf.get_uuid()
		       );
    } else {
      corpus->doc_str(vm["unified-doc-string"].as<std::string>());
    }
    if(vm.count("train")) {
      if(vm["train-from-db"].as<bool>()) {
	ERROR << "Cannot set both --train and --train-from-db";
	return 1;
      }
      INFO << "Going to read TRAINING from " << vm["train"].as<std::string>();
      corp_disk_reader
	(
	 vm["train"].as<std::string>(),
	 vm["train-number"].as<int>(),
	 vm["num-dependency-hops"].as<unsigned int>(),
	 tools,
	 gov_vocab,
	 rel_vocab,
	 corpus
	 );
      CHECK_CORPUS_PTR(corpus);
      corpus->multithreaded( main_redis_connection->multithreaded() );
      std::vector<std::string> vnames = {
	vm["gov-vocab"].as<std::string>(), 
	vm["rel-vocab"].as<std::string>()
      };
      if(unify_train_vocab) {
	std::vector<SVocab*> vocabs = {&gov_vocab, &rel_vocab};
	minsky::SyntacticEDocVocabUpdater sevu;
	corpus->unify_vocabs(vm["train-name"].as<std::string>() + "::unified",
			     vnames, vocabs, sevu,
			     vm["unified-doc-string"].as<std::string>());
      } else {
	corpus->save_vocab(vnames[0], gov_vocab);
	corpus->save_vocab(vnames[1], rel_vocab);
      }
    } else {
      INFO << "Connecting to existing database at " << host << ":" << port;
      CHECK_CORPUS_PTR(corpus);
      corpus->multithreaded( main_redis_connection->multithreaded() );
      INFO << "... done connecting to existing db";

      INFO << "Going to load gov vocab [" << vm["gov-vocab"].as<std::string>() << "]from existing database at " << host << ":" << port;
      corpus->_load_vocab(vm["gov-vocab"].as<std::string>(), gov_vocab);
      INFO << ".... done loading gov vocab";

      INFO << "Going to load rel vocab [" << vm["rel-vocab"].as<std::string>() << "]from existing database at " << host << ":" << port;
      corpus->_load_vocab(vm["rel-vocab"].as<std::string>(), rel_vocab);
      INFO << ".... done loading rel vocab";

      INFO << "Going to read TRAINING \"" << vm["train-name"].as<std::string>() << "\" from existing database at " << host << ":" << port << " under field " << corpus->doc_str();
      size_t num_loaded = corpus->load_from_key(vm["train-name"].as<std::string>(),
						true, // use the ferrum::REDIS_CORPUS_DOC_LIST as a prefix
						false, // don't print the doc ids of all the ones we load
						vm["train-number"].as<int>());
      INFO << "Loaded " << num_loaded << " documents from " << vm["train-name"].as<std::string>();
    }

    if(vm["shuffle-train"].as<bool>()) {
      INFO << "Shuffling the training corpus";
      corpus->shuffle();
    }


    if(vm["no-run"].as<bool>()) {
      // do nothing
    } else {
      ferrum::StringDiscreteKindPrinter inf_printer;
      //inf_printer.gk_vocab = &gov_kind_vocab;
      //inf_printer.rk_vocab = &rel_kind_vocab;
      inf_printer.g_vocab  = &gov_vocab;
      inf_printer.r_vocab  = &rel_vocab;
      double verb_bck_val = 0.0;
      if(vm.count("uniform-verb-background")) {
	verb_bck_val = vm["uniform-verb-background"].as<double>();
      }
      std::shared_ptr<std::vector<double> > gov_background(new std::vector<double>(gov_vocab.num_words(), verb_bck_val));
      double arc_bck_val = 0.0;
      if(vm.count("uniform-arc-background")) {
	arc_bck_val = vm["uniform-arc-background"].as<double>();
      }
      std::shared_ptr<std::vector<double> > rel_background(new std::vector<double>(rel_vocab.num_words(), arc_bck_val));
      ferrum::MinskyEntityCounter mec(minsky::AnnotationLevel::SYNTAX);
      ferrum::RedisBackgroundComputer<Doc, SVocab> rbc;
      ferrum::EDocBackgroundCounter<ferrum::RedisCorpus<Doc>, ferrum::RedisBackgroundComputer<Doc, SVocab> > bec;
      // if( ! vm.count("uniform-verb-background") ) {
      // 	INFO << "Computing predicate/governor background model...";
      // 	bec.in_memory_as(gov_background)
      // 	  .defined_by(&gov_vocab)
      // 	  .how(&mec)
      // 	  .with(&rbc)
      // 	  .over(corpus)
      // 	  .compute_predicate_background<std::string>();
      // 	INFO << "... done computing predicate/governor background model...";
      // } else {
      // 	INFO << "Keeping the predicate/governor background model at a constant " << verb_bck_val;
      // }

      // if( ! vm.count("uniform-arc-background") ) {
      // 	INFO << "Computing relation background model...";
      // 	bec.in_memory_as(rel_background)
      // 	  .defined_by(&rel_vocab)
      // 	  .compute_relation_background<std::string>();
      // 	INFO << "... done computing relation background model...";
      // } else {
      // 	INFO << "Keeping the relation background model at a constant " << arc_bck_val;
      // }

      Inferencer inf(num_templates, num_slots, hyperparams);
      inf.sampling_strategy(&strategy);

      INFO << "Initializing training model";
      // int num_mentions = ferrum::get_num_mentions(corpus);
      inf.init(corpus, gov_vocab, rel_vocab);
      INFO << "... done initializing training model";

      ferrum::RedisCorpus<Doc>* heldout_corpus = NULL;
      gov_vocab.allow_new_words(false);
      rel_vocab.allow_new_words(false);
      bool do_dev = read_dev(vm, gov_vocab, rel_vocab, &heldout_corpus, tools, main_redis_connection, -1, corp_disk_reader);
      // SAMPLING LOOP
      ferrum::DKVWriters sw_wrapper;
      inf_printer.g_vocab = &gov_vocab;
      inf_printer.print.tg = vm["print-templates-every"].as<int>();
      inf_printer.print.usage_t = vm["print-usage-every"].as<int>();
      //int h_epoch = 0;
      int epoch = 0;
      do {
	std::shared_ptr<ferrum::CompressedTar> output_tar;
	if(vm.count(names.INFERENCER_SERIALIZATION)) {
	  save_state(inf, gov_vocab, rel_vocab, output_tar,
		     vm[names.INFERENCER_SERIALIZATION].as<std::string>(),
		     epoch);
	}
	INFO << "Starting learning epoch " << epoch;
	if(output_tar == NULL) {
	  inf.learn< ferrum::RedisCorpus<Doc>,
		     ::ferrum::db::RedisThriftSmartWriter,
		     decltype(label_redis_connection) >
	    (
	     corpus,
	     epoch,
	     inf_printer,
	     &sw_wrapper,
	     false, // not heldout
	     label_redis_connection
	     );
	} else {
	  inf.learn< ferrum::RedisCorpus<Doc>,
		     ::ferrum::TarThriftSmartWriter,
		     std::shared_ptr<ferrum::CompressedTar> >
	    (
	     corpus,
	     epoch,
	     inf_printer,
	     &sw_wrapper,
	     false, // not heldout
	     output_tar
	     );
	}
	INFO << "Done with inference in epoch " << epoch;
	if(do_dev) {
	  ERROR << "Heldout not yet implemented :( ";
	  // ferrum::VStrategy heldout_strategy(strategy);
	  // heldout_strategy.heldout = true;
	  // INFO << "Starting heldout epoch " << h_epoch;
	  // ferrum::Timer h_timer("heldout_inference");
	  // var_inf.learn< ferrum::RedisCorpus<Doc>,
	  // 		 ::ferrum::db::RedisThriftSmartWriter,
	  // 		 decltype(main_redis_connection) >
	  //   (
	  //    heldout_corpus,
	  //    heldout_strategy,
	  //    h_epoch,
	  //    inf_printer,
	  //    &sw_wrapper,
	  //    label_redis_connection
	  //    );
	  // INFO << "Done with heldout inference";
	  // ++h_epoch;
	}
	++epoch;
      } while(epoch < num_epochs_); // end training epoch loop
      if(vm.count(names.INFERENCER_SERIALIZATION)) {
	std::shared_ptr<ferrum::CompressedTar> output_tar;
	save_state(inf, gov_vocab, rel_vocab, output_tar,
		   vm[names.INFERENCER_SERIALIZATION].as<std::string>(),
		   epoch);
      }
      FREE_IF_NNULL(heldout_corpus);
    }
    FREE_IF_NNULL(corpus);

  } // else if(vm.count("serialized-inferencer")) {
  //   // 0. create redis connection
  //   // if(! vm.count("host") || !vm.count("port") ) {
  //   //   ERROR << "Both --host and --port must be set";
  //   //   return 1;
  //   // }
  //   // std::string host(vm["host"].as<std::string>());
  //   // unsigned short port(vm["port"].as<unsigned short>());
  //   // ferrum::db::Address addr(host, port);
  //   // std::shared_ptr<ferrum::db::Redis> main_redis_connection(new ferrum::db::Redis(addr));

  //   // 1. read in model
  //   ferrum::CompressedTar ct(vm["serialized-inferencer"].as<std::string>());
  //   ct.read();
  //   minsky::residual::ResidualGlobalSlots rgs;
  //   for(const auto& entry : ct) {
  //     std::string name(ct.name(entry));
  //     size_t size = ct.size(entry);
  //     void* buffer = malloc(size);
  //     if(! buffer) {
  // 	ERROR << "Unable to malloc size " << size << " while reading " << name;
  // 	throw std::bad_alloc();
  //     }
  //     TRACE << "sizes: " << size;
  //     ct.read(buffer, size);
  //     if(size == 0) {
  // 	WARN << "Unable to read anything for inference state " << name;
  // 	break;
  //     }
  //     ferrum::thrift::thrift_struct_from_buffer<ferrum::thrift::TCompactProtocol>
  // 	(
  // 	 buffer,
  // 	 size,
  // 	 &rgs
  // 	 );
  //     free(buffer);
  //     break;
  //   }

  //   // 2. read-in vocabulary
  //   if(rgs.vocabularies.size() != 2) {
  //     ERROR << "Reading " << rgs.vocabularies.size() << " vocabularies, not 2";
  //     throw 5;
  //   }
  //   SVocab gov_vocab;
  //   (void)gov_vocab.reset_with(rgs.vocabularies[0]);
  //   SVocab rel_vocab;
  //   (void)rel_vocab.reset_with(rgs.vocabularies[1]);
    
  //   ferrum::StringDiscreteKindPrinter hinf_printer;
  //   //inf_printer.gk_vocab = &gov_kind_vocab;
  //   //inf_printer.rk_vocab = &rel_kind_vocab;
  //   hinf_printer.g_vocab  = &gov_vocab;
  //   hinf_printer.r_vocab  = &rel_vocab;

  //   // 3. read-in background counts
  //   std::shared_ptr<std::vector<double> > gov_background(new std::vector<double>(gov_vocab.num_words(), 0.0));
  //   std::copy(rgs.predicate_background.begin(), rgs.predicate_background.end(), gov_background->begin());
  //   std::shared_ptr<std::vector<double> > rel_background(new std::vector<double>(rel_vocab.num_words(), 0.0));
  //   std::copy(rgs.relation_background.begin(), rgs.relation_background.end(), rel_background->begin());

  //   // 4. fill shp with hyperparams from read-in model, and create model
  //   //Model dm(num_templates, num_slots);
  //   Inferencer heldout_inf(num_templates, num_slots);
  //   if(rgs.__isset.predicate_hyper) {
  //     heldout_inf.hyper_gov(rgs.predicate_hyper);
  //   } else {
  //     heldout_inf.hyper_gov(gov_vocab.num_words(), hyperparams.verb_hyper);
  //   }
  //   if(rgs.__isset.relation_hyper) {
  //     heldout_inf.hyper_rel(rgs.relation_hyper);
  //   } else {
  //     heldout_inf.hyper_rel(rel_vocab.num_words(), hyperparams.arc_hyper);
  //   }
  //   if(rgs.__isset.slot_hyper) {
  //     heldout_inf.hyper_slot(rgs.slot_hyper);
  //   } else {
  //     heldout_inf.slot_usage_hyper_val(hyperparams.slot_usage);
  //   }
  //   heldout_inf.template_usage_hyper_val(hyperparams.template_usage);
    

  //   // 5. create strategy
  //   ferrum::VStrategy heldout_strategy(strategy);
  //   heldout_strategy.heldout = true;

  //   // if(vm.count("templates-to-file")) {
  //   //   ferrum::DiscreteKindPrinter<ferrum::Vocabulary<string>, ferrum::Vocabulary<string>,
  //   // 				  ferrum::Vocabulary<string>, ferrum::Vocabulary<string> > dkp;
  //   //   dkp.gk_vocab = &gkv;
  //   //   dkp.rk_vocab = &rkv;
  //   //   dkp.g_vocab = &gv;
  //   //   dkp.r_vocab = &rv;
  //   //   dm.print_templates(vm["templates-to-file"].as<std::string>(), &dkp);
  //   // }
  //   ferrum::RedisCorpus<Doc>* heldout_corpus = NULL;
  //   bool do_dev = read_dev(vm, gov_vocab, rel_vocab, &heldout_corpus, tools, main_redis_connection, corp_disk_reader);
  //   if(do_dev) {
  //     // load corpus
      
  //     bool need_to_infer = vm.count("compute-dev-label");
  //     ferrum::DKVWriters sw_wrapper;
  //     //ferrum::StringDiscreteKindPrinter sw_wrapper;
  //     if(need_to_infer) {
  // 	// create sw_wrappers
  //     }
  //     need_to_infer |= vm["compute-dev-ll"].as<bool>();
  //     if(need_to_infer) {	
  // 	ferrum::Timer htimer("heldout_infer");
  // 	heldout_inf.learn< ferrum::RedisCorpus<Doc>,
  // 			   ::ferrum::db::RedisThriftSmartWriter,
  // 			   decltype(main_redis_connection) >
  // 	  (
  // 	   heldout_corpus,
  // 	   heldout_strategy,
  // 	   0, // epoch ID
  // 	   hinf_printer,
  // 	   &sw_wrapper,
  // 	   label_redis_connection
  // 	   );
  //   	INFO << "Done with sampling heldout";
  //     } // end if(need_to_sample)
  //     if(vm["compute-dev-ll"].as<bool>()) {
  // 	double heldout_ll = 0.0;
  //   	// heldout_ll =
  //   	//   dm.latent_and_kind_marginalized_loglikelihood(heldout_corpus, gv, rv,
  //   	// 						heldout_sampler.doc_template_params());
  //   	INFO << "Held-out marginalized ll: " << heldout_ll;
  //     }
  //     // // if(vm.count("compute-coherence")) {
  //     // // 	if(!vm.count("coherence-file")) {
  //     // // 	  ERROR << "You must supply a coherence file";
  //     // // 	  throw 5;
  //     // // 	}
  //     // // 	std::vector<std::string> which_coherences = vm["compute-coherence"].as<std::vector< std::string> >();
  //     // // 	std::ofstream myfile;
  //     // // 	myfile.open(vm["coherence-file"].as<std::string>());
  //     // // 	myfile << "which\tid\tcoherence\twhich_avg\n";
  //     // // 	for(const std::string& which : which_coherences) {
  //     // // 	  ERROR << which;
  //     // // 	  std::vector<double> coherences;
  //     // // 	  if(which == "gov_kind") {
  //     // // 	    coherences = 
  //     // // 	      dm.compute_coherences< Doc, string>(vm["coherence-M"].as<int>(), which, heldout_corpus, gkv);
  //     // // 	    int id = 0;
  //     // // 	    double avg = ferrum::average(coherences);
  //     // // 	    for(const double coher : coherences) {
  //     // // 	      myfile << which << "\t" << (id++) << "\t" << coher << "\t" << avg << "\n";
  //     // // 	    }
  //     // // 	  }
  //     // // 	  if(which == "gov_obs") {
  //     // // 	    coherences = 
  //     // // 	      dm.compute_coherences< Doc, string>(vm["coherence-M"].as<int>(), which, heldout_corpus, gv);
  //     // // 	    int id = 0;
  //     // // 	    double avg = ferrum::average(coherences);
  //     // // 	    for(const double coher : coherences) {
  //     // // 	      myfile << which << "\t" << (id++) << "\t" << coher << "\t" << avg << "\n";
  //     // // 	    }
  //     // // 	  }
  //     // // 	  if(which == "gov_obs_kind_marginalized") {
  //     // // 	    coherences = 
  //     // // 	      dm.compute_coherences< Doc, string>(vm["coherence-M"].as<int>(), which, heldout_corpus, gv);
  //     // // 	    int id = 0;
  //     // // 	    double avg = ferrum::average(coherences);
  //     // // 	    for(const double coher : coherences) {
  //     // // 	      myfile << which << "\t" << (id++) << "\t" << coher << "\t" << avg << "\n";
  //     // // 	    }
  //     // // 	  }
  //     // // 	}
  //     // // 	myfile.close();
  //     // // } // end if(compute-coherences)
  //   } // end if(dev)
  // } // end if(serialized-inferencer)
  return 0;
}
