#!/usr/bin/env python

import argparse
#import minsky
from autoload_minsky import minsky
import codecs
import concrete
import concrete.util
import redis
import os
from tube.util import parse_redis_addr
import random
import sys

from thrift.transport.TTransport import TMemoryBuffer
from concrete.util.thrift_factory import factory
import collections

EDOC = 'edoc'
SIMPLEDOC = 'simpledoc'

parser = argparse.ArgumentParser()
parser.add_argument('addr', type = parse_redis_addr, help = "Redis host:port")
parser.add_argument('which', choices = (EDOC, SIMPLEDOC))
parser.add_argument('name', nargs = '+', type = str, help = "Redis key")
parser.add_argument('--N', type = int, default = 10)
parser.add_argument('--from-list', action = 'store_true')
parser.add_argument('--name-prefix', type = str, default = 'document:')
parser.add_argument('--random', action = 'store_true')
parser.add_argument('--labels', type = str, default = None)
parser.add_argument('--labels-sep', type = str, default = '\t')
parser.add_argument('--field', type = str, default = None)
parser.add_argument('--predicate-vocab', type = str, default = None)
parser.add_argument('--relation-vocab', type = str, default = None)
parser.add_argument('--sem-predicate-vocab', type = str, default = None)
parser.add_argument('--sem-relation-vocab', type = str, default = None)
parser.add_argument('--word-vocab', type = str, default = None)
parser.add_argument('--feature-value-directory', type = str, default = None,
                    help = 'directory to store liblinear-style output; '
                           'currently only applicable when `which` == %s' % (SIMPLEDOC))
parser.add_argument('--feature-value-data', type = str, default = 'train.data',
                    help = 'basename of file to store liblinear-style output; '
                           'currently only applicable when `which` == %s' % (SIMPLEDOC))
# parser.add_argument('--feature-value-weight', choices = ('count','binary','tfidf'),
#                     default = 'count',
#                     help = 'feature values to produce for liblinear-style output')
args = parser.parse_args()

labels = collections.defaultdict(lambda : '1')
cats = {}
label_arr = [] 
if not args.labels is None:
    with codecs.open(args.labels, 'r', 'utf-8') as fp:
        for line in fp.xreadlines():
            xline = line.strip().split(args.labels_sep)
            if not len(xline) == 2:
                sys.stderr.write("Found line with more than one separating \"%s\"\n" % args.labels_sep)
                sys.stderr.write("%s" % line)
                exit(1)
            if not xline[1] in cats:
                label_arr.append(xline[1])
                cats[xline[1]] = str(len(label_arr)) 
            labels[xline[0]] = cats[xline[1]]
#            print "%s --> %s" % (xline[0], xline[1])
    if args.feature_value_directory:
        lname = os.path.join(args.feature_value_directory, 'labels')
        print "Writing labels to %s" % (lname)
        with codecs.open(lname, 'w', 'utf-8') as lp:
            for idx, l in enumerate(label_arr):
                lp.write(l)
                lp.write("\n")
                

def show_edoc(mdoc, ow, gvocab = None, rvocab = None, gkv = None, rkv = None, fv = False):
    if mdoc.entities is None:
        ow.write("EDoc id %s has 0 entities\n" % ( mdoc.id ))
        return
    ow.write("EDoc id %s has %d entities\n" % ( mdoc.id, len(mdoc.entities) ))
    
    for entity in mdoc.entities:
        if entity.mentions is None:
            continue
        ow.write("{ %d\n" % (len(entity.mentions)))
        for mention in entity.mentions:
            tups = []
            for pa in mention.structures:
                wgv = None
                wrv = None
                if pa.annot_level == minsky.AnnotationLevel.SYNTAX:
                    wgv = gvocab
                    wrv = rvocab
                elif pa.annot_level == minsky.AnnotationLevel.SEMANTIC:
                    wgv = gkv
                    wrv = rkv
                t = (pa.predicate.word if wgv is None else wgv[pa.predicate.word],
                     pa.relation if wrv is None else wrv[pa.relation],
                     minsky.AnnotationLevel._VALUES_TO_NAMES[pa.annot_level])
                tups.append(t)
            ow.write("\t%s %s\n" % (str(tups), mention.location_id if mention.location_id else ""))
        ow.write("}\n")


def show_sdoc(mdoc, ow, vocab = None, fv = False):
    if mdoc.sentences is None:
        if fv:
            pass
        else:
            ow.write("SimpleDoc id %s has 0 \"sentences\"\n" % ( mdoc.id ))
        return
    if not fv:
        ow.write("SimpleDoc id %s has %d sentences\n" % ( mdoc.id, len(mdoc.sentences) ))
    if fv:
        ow.write("%s" % (labels[mdoc.id]))
    for area in mdoc.sentences:
        if not area.words is None:
            if not fv:
                ow.write("[ %d\n" % (len(area.words)))
            if fv:
                from collections import Counter
                c = Counter(area.words)
                for w in sorted(c.keys()):
                    ow.write(" %d:%d" % (w, c[w]))
                ow.write(" # doc %s\n" % mdoc.id)
            else:
                tups = [word if vocab is None else vocab[word] for word in area.words]
                ow.write("\t%s\n" % (str(tups)))
                
        if not area.counts is None and not area.counts.icounts is None:
            if fv:
                for word in sorted(area.counts.icounts.keys()):
                    ow.write(" %d:%d" % (word, area.counts.icounts[word]))
                ow.write(" # doc %s\n" % mdoc.id)
            else:
                tups = [word if vocab is None else vocab[word] for (word, count) in area.counts.icounts.iteritems() for i in xrange(count)]
                ow.write( "[ %d\n" % (len(tups)))
                ow.write("\t%s\n" % (str(tups)))
        if not fv:
            ow.write("]\n")

def read_from_buffer(buf, which, name = None):
    transport_in = TMemoryBuffer(buf)
    protocol_in = factory.createProtocol(transport_in)
    obj = None
    if which == EDOC:
        obj = minsky.EDoc()
    elif which == SIMPLEDOC:
        obj = minsky.SimpleDoc()
    try:
        obj.read(protocol_in)
    except Exception as e:
        import sys
        sys.stderr.write("%s Error reading buffer %s\n" % (("With key " + name) if not name is None else "", str(buf)))
        raise e
    return obj

def get_vocab(redis_db, key):
    if rdb.exists(key):
        voc = rdb.hgetall(key)
        return (voc['word_list'].strip().split(), False)
    else:
        with codecs.open(key, 'r', 'utf-8') as vfp:
            return ([x.strip() for x in vfp.xreadlines()], True)
    
rdb = redis.Redis(*args.addr)
vocabs = []
existing_vocab = False
if args.predicate_vocab:
    (v, existing_vocab) = get_vocab(rdb, args.predicate_vocab)
    vocabs.append(v)
if args.relation_vocab:
    (v, existing_vocab) = get_vocab(rdb, args.relation_vocab)
    vocabs.append(v)
if args.sem_predicate_vocab:
    (v, existing_vocab) = get_vocab(rdb, args.sem_predicate_vocab)
    vocabs.append(v)
if args.sem_relation_vocab:
    (v, existing_vocab) = get_vocab(rdb, args.sem_relation_vocab)
    vocabs.append(v)
if args.word_vocab:
    (v, existing_vocab) = get_vocab(rdb, args.word_vocab)
    vocabs.append(v)
vocab_writer = None
out_writer = sys.stdout
fv = False
if not args.feature_value_directory is None and not existing_vocab:
    vname = os.path.join(args.feature_value_directory,
                         'vocab')
    vocab_writer = codecs.open(vname, 'w', 'utf-8')
    print "Writing vocab list to %s" % vname
    for (vidx, word) in enumerate(vocabs[0]):
        if vidx == 0:
            continue
        vocab_writer.write(word.decode('utf-8'))
        vocab_writer.write("\n")
    out_writer = codecs.open(os.path.join(args.feature_value_directory,
                                          args.feature_value_data), 'w', 'utf-8')
    fv = True

def process_mobj(obj, args, out_writer, vocabs, fv):
    if args.which == EDOC:
        show_edoc(obj, out_writer, *vocabs, fv = fv)
    elif args.which == SIMPLEDOC:
        show_sdoc(obj, out_writer, *vocabs, fv = fv)

def process_input(func):
    num_seen = 0
    if args.from_list:
        lst = []
        if args.field is None:
            pass
        else:
            lst = [x for name in args.name for x in rdb.lrange(name, 0, min(rdb.llen(name), args.N)) ]
        for d in lst:
            dd = args.name_prefix + d
            obuff = rdb.get(dd) if args.field is None else rdb.hget(dd, args.field)
            if obuff is None:
                continue
            obj = read_from_buffer( obuff, args.which, args.name_prefix + d )
            #print "Processing %s" % obj.id
            func(obj, args, out_writer, vocabs, fv)
            num_seen += 1
    else:
        for x in args.name:
            if num_seen >= args.N:
                break
            key_list = rdb.keys(x)
            if args.random:
                random.shuffle(key_list)
            for y in key_list:
                if num_seen >= args.N:
                    break
                ty = rdb.type(y)
                #print y, ty
                if ty == 'string' or \
                   (ty == 'hash' and not args.field is None):
                    if num_seen >= args.N:
                        break
                    obuff = rdb.get(y) if ty == 'string' else rdb.hget(y, args.field)
                    if obuff is None:
                        continue
                    obj = read_from_buffer( obuff, args.which, y )
                    #print "Processing %s" % obj.id
                    func(obj, args, out_writer, vocabs, fv)
                    num_seen += 1
                else:
                    rr = concrete.util.RedisReader(rdb, y, key_type = ty)
                    for obj_buf in rr:
                        if num_seen >= args.N:
                            break
                        obj = read_from_buffer(obj_buf, args.which, y)
                        func(obj, args, out_writer, vocabs, fv)
                        num_seen += 1

process_input(process_mobj)

