#!/bin/bash

set -e
set -o xtrace

git subtree push -P ferrum bb-ferrum master
git subtree push -P minsky bb-minsky master
