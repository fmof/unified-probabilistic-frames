#!/usr/bin/env python

import argparse
import codecs
from autoload_minsky import minsky
import concrete
import concrete.util
import tarfile
import sys
import numpy

from thrift import TSerialization
from concrete.util.thrift_factory import factory as thrift_factory
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('name', type = str, help = "path to Thrift RGS file")
parser.add_argument('entry', type = str, help = 'name of entry in tgz to load')
parser.add_argument('out_m', type = str, help = 'name of output matrix')
parser.add_argument('out_l', type = str, help = 'name of output labels')
parser.add_argument('--N', type = int, default = 10)
parser.add_argument('--M', type = int, default = 10)
parser.add_argument('--sep', type = str, default = ' ')
parser.add_argument('--with-values', dest = 'print_vals', action = 'store_true', default = True)
parser.add_argument('--without-values', dest = 'print_vals', action = 'store_false')
parser.add_argument('--highest-weight', dest = 'top_weights', action = 'store_true', default = True)
parser.add_argument('--lowest-weight', dest = 'top_weights', action = 'store_false')
parser.add_argument('--with-slot-dists', action = 'store_true')
parser.add_argument('--select-slots', type = int, default = 0,
                    help = 'print the N most probable slot emissions for the top K slots (per template) ')
parser.add_argument('--data-frame', action = 'store_true',
                    help = "if true, don't print out intermediate loading messages (default: false)")
parser.add_argument('--no-sort', dest = 'sort', action = 'store_false')
args = parser.parse_args()

rus = minsky.ResidualUniqueSlots()
fname = args.name
N = args.N
if fname.endswith('.tar.gz') or fname.endswith('.tgz'):
    tar = tarfile.open(fname, 'r|gz')
    while True:
        tarinfo = tar.next()
        if tarinfo is None:
            sys.stderr.write("No entries in archive %s\n" % fname)
            exit(1)
        if not tarinfo.isfile():
            continue
        if not tarinfo.name == args.entry:
            continue
        if not args.data_frame:
            print "Loading %s from %s..." % (args.entry, fname)
        rus = TSerialization.deserialize(rus,
                                         tar.extractfile(tarinfo).read(),
                                         protocol_factory=thrift_factory.protocolFactory)
        break
else:
    if not args.data_frame:
        print "Loading %s..." % (fname)
    rus = concrete.util.read_thrift_from_file(rus, fname)
if not args.data_frame:
    print "Done"

def which(t, d):
    if d == 'predicates':
        return t.predicate_frame
    else:
        return t.role



        
#back
#writer = codecs.getwriter('utf-8')(sys.stdout)
writer = sys.stdout
num_in_m = defaultdict(int)
sorter = (lambda x : -x[1]) if args.top_weights else (lambda x : x[1])

if rus.semantics is None:
    mat = []
    voc_idx = None
    for template in rus.thematics:
        phi = template.predicate_frame.distr.weights
        if phi.normalized is not None:
            phi = phi.normalized
        elif phi.unnormalized is not None:
            sphi = sum(phi.unnormalized)
            phi = [float(x)/sphi for x in phi]
        else:
            print "unsure how to handle neither weights.{normalized,unnormalized} being set"
            exit(1)
        if voc_idx is None:
            voc_idx = template.predicate_frame.distr.vocab_idx
        else:
            if not template.predicate_frame.distr.vocab_idx == voc_idx:
                print "vocab indices don't match up"
                exit(1)
        mat.append(phi)
    mat = numpy.array(mat).transpose()
    numpy.savetxt(args.out_m, mat)
    vocab = rus.vocabularies[voc_idx].words
    with codecs.open(args.out_l, 'w', 'utf-8') as ol:
        for (wi, w) in enumerate(vocab):
            ol.write("%d\n" % (wi))
else:
    tf = []
    voc_idx = None
    for template in rus.thematics:
        phi = template.predicate_frame.distr.weights
        if phi.normalized is not None:
            phi = phi.normalized
        elif phi.unnormalized is not None:
            sphi = sum(phi.unnormalized)
            phi = [float(x)/sphi for x in phi]
        else:
            print "unsure how to handle neither weights.{normalized,unnormalized} being set"
            exit(1)
        tf.append(phi)
    tf = numpy.array(tf)
    fv = []
    for semantic in rus.semantics.preds:
        phi = semantic.predicate_frame.distr.weights
        if phi.normalized is not None:
            phi = phi.normalized
        elif phi.unnormalized is not None:
            sphi = sum(phi.unnormalized)
            phi = [float(x)/sphi for x in phi]
        else:
            print "unsure how to handle neither weights.{normalized,unnormalized} being set"
            exit(1)
        fv.append(phi)
        if voc_idx is None:
            voc_idx = semantic.predicate_frame.distr.vocab_idx
        else:
            if not semantic.predicate_frame.distr.vocab_idx == voc_idx:
                print "vocab indices don't match up"
                exit(1)
    mat = tf.dot(numpy.array(fv)).transpose()
    numpy.savetxt(args.out_m, mat)
    vocab = rus.vocabularies[voc_idx].words
    with codecs.open(args.out_l, 'w', 'utf-8') as ol:
        for (wi, w) in enumerate(vocab):
            ol.write("%d\n" % (wi))
    for template in rus.thematics:
        pass
