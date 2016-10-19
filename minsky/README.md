# minsky
`minsky` is a schema and multi-language library designed around event models.
It is based on [Thrift](https://thrift.apache.org/).

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
