#!/bin/bash

set -o nounset

echo $1 > /dev/null

NAMESPACE=${2:-}
if [ -z "${NAMESPACE}" ]; then
    NAMESPACE=ferrum
fi
LIBNAME=${3:-}
if [ -z "${LIBNAME}" ]; then
    LIBNAME=FERRUM
fi
INCLUDE_NAME=${4:-}
if [ -z "${INCLUDE_NAME}" ]; then
    INCLUDE_NAME="version.hpp"
fi

git diff --exit-code > /dev/null
is_unstaged=$?

git diff --cached --exit-code > /dev/null
is_ncommittted=$?

: ${REQUIRE_STAGED=1}
: ${REQUIRE_COMMITTED=0}

# if [ ! $is_ncommittted == 0 ] && [ $REQUIRE_STAGED == 1 ]; then
#     echo "Cannot build an executable because there are unstaged changes."
#     exit 1
# fi
if [ ! $is_ncommittted == 0 ] && [ $REQUIRE_COMMITTED == 1 ]; then
    echo "Cannot build an executable because there are non-committed changes."
    exit 1
fi

curr_hash=$(git rev-parse HEAD)
curr_date=$(date +%Y%m%d-%H%M%S)

cat <<EOF > $1
#include "${INCLUDE_NAME}"

namespace ${NAMESPACE} {
  const std::string ${LIBNAME}_BUILD_DATE="$curr_date";
  const std::string ${LIBNAME}_GIT_SHA="$curr_hash";
}

EOF
