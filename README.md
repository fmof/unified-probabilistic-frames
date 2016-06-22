## Dependencies

* boost (recent, works with >= 1.56)
* Thrift >= 0.9.3
* concrete >= 4.8 < 5
* gsl == 1.16
* [libLBFGS](http://www.chokkan.org/software/liblbfgs/) == [1.10](https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz)
  - [Github](https://github.com/chokkan/liblbfgs)
* GoogleLOG (GLOG) == [0.3.3](https://github.com/google/glog/archive/v0.3.3.tar.gz)
* [libarchive](https://github.com/libarchive/libarchive/releases) == [3.1.2](https://github.com/libarchive/libarchive/archive/v3.1.2.tar.gz)
* [hiredis](https://github.com/redis/hiredis) == [0.13.3](https://github.com/redis/hiredis/archive/v0.13.3.tar.gz)
* [hihiredis](https://gitlab.hltcoe.jhu.edu/fferraro/hihiredis)


* google test (included in repo)
* eigen (included in repo)

Other:
* [redis](http://redis.io/) >= [3.0.0](http://download.redis.io/releases/redis-3.0.5.tar.gz)

## Subtree stuff

### Prereqs
```
git remote add bb-ferrum git@bitbucket.org:fmof/ferrum.git
git remote add bb-minsky git@bitbucket.org:fmof/minsky.git
```

### Contributing upstream
```
git subtree push --prefix=ferrum bb-ferrum master
git subtree push --prefix=minsky bb-minsky master
```

### Pulling from upstream
```
git subtree pull --prefix=ferrum bb-ferrum master
git subtree pull --prefix=minsky bb-minsky master
```

### If stuff gets messed up and I need to resplit
```
git subtree split --prefix ferrum -b sync-ferrum
git remote add bb-ferrum git@bitbucket.org:fmof/ferrum.git
git push bb-ferrum sync-ferrum:master
git rm -r ferrum/
rm -rf ferrum/
git add -A
git commit -am "removed ferrum/ (temp)"
git subtree add --prefix=ferrum bb-ferrum master
```

```
git subtree split --prefix minsky -b sync-minsky
git remote add bb-minsky git@bitbucket.org:fmof/minsky.git
git push bb-minsky sync-minsky:master
git rm -r minsky/
rm -rf minsky/
git add -A
git commit -am "removed minsky/ (temp)"
git subtree add --prefix=minsky bb-minsky master
```
