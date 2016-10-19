# UPF --- A Unified Bayesian Model of Scripts, Frames and Language

Version 1.0.0

This is a C++ 11 implementation of [A Unified Bayesian Model of Scripts, Frames and Language](http://cs.jhu.edu/~ferraro/papers/ferraro-upf-2016.pdf) (Ferraro and Van Durme, 2016).
It has been developed and tested on Linux x86_64, under G++ >= 4.8.4.
It should also compile on a Mac, though a couple Makefile changes may need to be made.

## Dependencies

### Hard Dependencies
The following are required to get both the UPF and baseline model to compile.
* boost (recent, works with >= 1.56)
* [Thrift](https://thrift.apache.org/download) >= 0.9.3
* gsl == 1.16
* cblas
* [libarchive](https://github.com/libarchive/libarchive/releases) == [3.1.2](https://github.com/libarchive/libarchive/archive/v3.1.2.tar.gz)
* [redis](http://redis.io/) >= [3.0.0](http://download.redis.io/releases/)
* [hiredis](https://github.com/redis/hiredis) == [0.13.3](https://github.com/redis/hiredis/archive/v0.13.3.tar.gz)
* GoogleLOG (GLOG) == [0.3.3](https://github.com/google/glog/archive/v0.3.3.tar.gz)

### Optional Dependencies
As this project shares code with other projects, there are some optional dependencies.
* [libLBFGS](http://www.chokkan.org/software/liblbfgs/) == [1.10](https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz)
  - [Github](https://github.com/chokkan/liblbfgs)
* [hihiredis](https://gitlab.hltcoe.jhu.edu/fferraro/hihiredis)
* atlas

There are a number of dependencies that are included in this repo directly.
They are compiled on-demand.
* Google test (included in repo)
* Eigen (included in repo)
* [Concrete](http://hltcoe.github.io/) >= 4.8 < 5
  - Concrete is a data schema, described in [Ferraro et al., 2014](http://cs.jhu.edu/~ferraro/papers/ferraro-concrete-2014.pdf).

## Building

### Local Configurations
Depending on where you installed the dependencies, you may need to update `Makefile.config`.
* If the headers are installed in `/usr/local/include`, and the shared objects are installed in `/usr/local/lib`, then you do _not_ have to change anything.
* If the headers (and shared objects) are installed in the same place (but not `/usr/local/{include,lib}`), then change lines 5 and 6 of `Makefile.config`.
* If the headers (shared objects) are installed in separate directories, then change lines 11-61 as appropriate.

### Building the targets

`make help` will display all known targets.
To build the UPF model, run `make models/upf_cgibbs_driver`.
To build the baseline model, run `make models/crtlda_cgibbs_basic_driver`.

There are a number of Makefile ENV variables that can be set to change compilation.
Some major ones are:
* `DEBUG={0,1}`: turn on debug compilation. `DEBUG=1` turns ON debugging, i.e., `-g -O0`). To run quickly, `DEBUG=0`. Default: 0
* `LINK_HOW={dynamic,static}`: use dynamic or static linking for certain libraries. Default: dynamic
* `LOG_AS_COUT`: This is meant to change what logging is used. It is by default undefined. Set it to anything to use stdout logging instead of Google Log.

## Data

Due to licensing issues, I cannot release the fully annotated input files.
They are Concrete Communications with:
* a Stanford dependency parse (collapsed-cc)
* Stanford part-of-speech and lemmatization tags
* Stanford entity coreference
* Semafor frame semantic parses.

The list of training ids is in `data/nyt10k.id_list.txt`.


## Licensing

This code is released under GPL v3.0.
Please contact me (ferraro [at] cs [dot] jhu [dot] edu) with any questions.
