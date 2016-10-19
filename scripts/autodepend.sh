#!/bin/sh
COMPILER="$1"
shift 1
#echo "${COMPILER}"
export DIR="$1"
#echo "DIR is ${DIR}"
shift 1
export PREFIX="$1"
shift 1
#echo "FILES are $@"
case "$DIR" in
"" | ".")
${COMPILER} "$@" #|
#sed -e ’s@ˆ\(.*\)\.o:@\1.d \1.o:@’
;;
*)
${COMPILER} "$@" |
perl -pe 's/^(.*)\.o\s*:/$ENV{"PREFIX"}\/$1\.d $ENV{"PREFIX"}\/$1\.o:/;'
#sed -e "s@ˆ\(.*\)\.o:@$DIR/\1.d $DIR/\1.o:@"
;;
esac
