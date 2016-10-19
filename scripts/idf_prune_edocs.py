#!/usr/bin/env python

import argparse
import minsky
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
import math

parser = argparse.ArgumentParser()
parser.add_argument('addr', type = parse_redis_addr, help = "Redis host:port")
parser.add_argument('name', nargs = '+', type = str, help = "Redis key")
parser.add_argument('--threshold', type = float, default = -1.0, help = 'idf threshold. idf(predicate) >= threshold for the mention to be kept')
parser.add_argument('--N', type = int, default = -1)
parser.add_argument('--from-list', action = 'store_true')
parser.add_argument('--name-prefix', type = str, default = 'document:')
parser.add_argument('--field', type = str, default = None)
parser.add_argument('--predicate-vocab', type = str, default = None)
parser.add_argument('--relation-vocab', type = str, default = None)
args = parser.parse_args()

def prune_edoc(mdoc, ow, al, weight, threshold, gvocab = None, rvocab = None):
  if mdoc.entities is None:
    ow.write("EDoc id %s has 0 entities\n" % ( mdoc.id ))
    return
  nents = []
  for entity in mdoc.entities:
    if entity.mentions is None:
      continue
    nments = []
    for mention in entity.mentions:
      found_level = False
      for struct in mention.structures:
        if struct.annot_level == al:
          found_level = True
          #print "Mention %s has weight %f vs threshold %f" % (str(mention), weight(struct.predicate.word), threshold)
          if weight(struct.predicate.word) >= threshold:
            #print "Adding mention %s" % (mention.id)
            nments.append(mention)
            break
    if len(nments) > 0:
      entity.mentions = nments
      nents.append(entity)
      #print "adding entity"
    if not found_level:
      print "No level %s found" % (str(al))
  mdoc.entities = nents

def read_from_buffer(buf, name = None):
    transport_in = TMemoryBuffer(buf)
    protocol_in = factory.createProtocol(transport_in)
    obj = minsky.EDoc()
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
        print "cannot find vocab key ", key
    
rdb = redis.Redis(*args.addr)
vocabs = []
existing_vocab = False
if args.predicate_vocab:
    (v, existing_vocab) = get_vocab(rdb, args.predicate_vocab)
    vocabs.append(v)
if args.relation_vocab:
    (v, existing_vocab) = get_vocab(rdb, args.relation_vocab)
    vocabs.append(v)
vocab_writer = None
out_writer = sys.stdout

def process_mobj(obj, out_writer, al, weight, threshold, vocabs):
    prune_edoc(obj, out_writer, al, weight, threshold, *vocabs)

annotation_level = minsky.AnnotationLevel.SYNTAX
    
def process_input(func):
    num_seen = 0
    if args.from_list:
        lst = []
        if args.field is None:
          print "Please provide a --field argument when joining from a list"
        else:
          lst = [(name, x) for name in args.name for x in rdb.lrange(name, 0, min(rdb.llen(name), args.N)) ]

        def computed_idf():
          N = 0
          idf = collections.defaultdict(float)
          print "computing IDF weights...."
          for (name, d) in lst:
            dd = args.name_prefix + d
            obuff = rdb.get(dd) if args.field is None else rdb.hget(dd, args.field)
            if obuff is None:
              continue
            obj = read_from_buffer(obuff, args.name_prefix + d)
            N += 1
            if obj.entities is None:
              continue
            words = set([s.predicate.word for entity in obj.entities for mention in entity.mentions for s in mention.structures if s.annot_level == annotation_level])
            for w in words:
              idf[w] += 1
          print "done computing idf weights"
          def curried(w):
            return math.log(N / (idf[w] if w in idf else 1E-8))
          #print [(i, curried(i)) for i in idf.keys()]
          return curried
        idf_func = computed_idf()
        for (di, (name, d)) in enumerate(lst):
          dd = args.name_prefix + d
          obuff = rdb.get(dd) if args.field is None else rdb.hget(dd, args.field)
          if obuff is None:
            continue
          obj = read_from_buffer(obuff, args.name_prefix + d)
          #print "Processing %s" % obj.id
          pnuments = 0 if obj.entities is None else len(obj.entities)
          func(obj, out_writer, annotation_level, \
               idf_func, args.threshold, \
               vocabs)
          if obj.entities is None or len(obj.entities) == 0:
            print "EDoc %s no longer has any entities; removing from corpus" % (obj.id)
            continue
          if (di + 1) % 100 == 0:
            print "Saving the %dth document" % (di + 1)
          print "Saving EDoc %s with %d entities (%d removed)" % (obj.id, len(obj.entities), pnuments - len(obj.entities))
          transport_out = TMemoryBuffer()
          protocol_out = factory.createProtocol(transport_out)
          obj.write(protocol_out)
          suffix = ":idf" + str(args.threshold)
          if args.field is None:
            rdb.set(args.name_prefix + d + suffix,
                    transport_out.getvalue())
          else:
            rdb.hset(dd, args.field + suffix, transport_out.getvalue())
            rdb.rpush(name + suffix, d)
          num_seen += 1
    else:
        print "ERROR: Not available"
        exit(1)
        for x in args.name:
            if args.N >= 0 and num_seen >= args.N:
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
                    func(obj, args, out_writer, vocabs)
                    num_seen += 1
                else:
                    rr = concrete.util.RedisReader(rdb, y, key_type = ty)
                    for obj_buf in rr:
                        if num_seen >= args.N:
                            break
                        obj = read_from_buffer(obj_buf, args.which, y)
                        func(obj, args, out_writer, vocabs)
                        num_seen += 1

process_input(process_mobj)

