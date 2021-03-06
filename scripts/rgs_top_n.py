#!/usr/bin/env python

import argparse
import codecs
from autoload_minsky import minsky
import concrete
import concrete.util
import tarfile
import sys

from thrift import TSerialization
from concrete.util.thrift_factory import factory as thrift_factory
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('name', type = str, help = "path to Thrift RGS file")
parser.add_argument('entry', type = str, help = 'name of entry in tgz to load')
parser.add_argument('--N', type = int, default = 10)
parser.add_argument('--M', type = int, default = 10)
parser.add_argument('--sep', type = str, default = ' ')
parser.add_argument('--with-values', dest = 'print_vals', action = 'store_true', default = True)
parser.add_argument('--without-values', dest = 'print_vals', action = 'store_false')
parser.add_argument('--highest-weight', dest = 'top_weights', action = 'store_true', default = True)
parser.add_argument('--lowest-weight', dest = 'top_weights', action = 'store_false')
parser.add_argument('--dist', choices = ['predicates', 'arcs'], default = 'predicates')
parser.add_argument('--with-slot-dists', action = 'store_true')
parser.add_argument('--select-slots', type = int, default = 0)
parser.add_argument('--data-frame', action = 'store_true',
                    help = "if true, don't print out intermediate loading messages (default: false)")
parser.add_argument('--no-sort', dest = 'sort', action = 'store_false')
args = parser.parse_args()

rgs = minsky.ResidualGlobalSlots()
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
        rgs = TSerialization.deserialize(rgs,
                                         tar.extractfile(tarinfo).read(),
                                         protocol_factory=thrift_factory.protocolFactory)
        break
else:
    if not args.data_frame:
        print "Loading %s..." % (fname)
    rgs = concrete.util.read_thrift_from_file(rgs, fname)
if not args.data_frame:
    print "Done"

def which(t, d):
    if d == 'predicates':
        return t.predicate_frame
    else:
        return t.role

vocab = (rgs.vocabularies[0 if args.dist == 'predicates' else 1]).words
#back
item = rgs.thematics if args.dist == 'predicates' else rgs.slots
#writer = codecs.getwriter('utf-8')(sys.stdout)
writer = sys.stdout
if not args.data_frame:
    print "vocab has %d elements" % (len(vocab))
num_in_m = defaultdict(int)
sorter = (lambda x : -x[1]) if args.top_weights else (lambda x : x[1])

for (topic_idx, topic) in enumerate(item):
    choice = which(topic, args.dist)
    if choice is None or choice.distr is None:
        continue
    tpred_sorted = list(enumerate(choice.distr.weights.residual))
    if args.sort:
        tpred_sorted.sort(key = sorter)
    max_to_pick = min(args.M, len(tpred_sorted)) if args.M >= 0 else len(tpred_sorted)
    for x in tpred_sorted[0:max_to_pick]:
        num_in_m[vocab[x[0]]] += 1
for (template_idx, template) in enumerate(item):
    choice = which(template, args.dist)
    if choice is None:
        if not args.data_frame:
            print "listing %d doesn't have a set frame" % (template_idx)
        continue
    if choice.distr is None:
        if not args.data_frame:
            print "listing %d doesn't have a set frame distribution" % (template_idx)
        continue
    lst = choice.distr.weights.residual
    if not len(lst) == len(vocab):
        if choice.distr.items is None or len(choice.dist.items) == 0:
            if not args.data_frame:
                print "listing %d has a residual not of the size of vocab (%d vs %d), but not Frame.Distribution.items set either" % (template_idx, len(lst), len(vocab))
            continue
        lst = zip(choice.distr.items, lst)
    else:
        lst = enumerate(lst)
    tpred_sorted = list(lst)
    if args.sort:
        tpred_sorted.sort(key = sorter)
    if not args.data_frame:
        print "%s %d" % ("Template" if args.dist == 'predicates' else "Slot", template_idx)
    if args.with_slot_dists and not args.data_frame:
        slot_d = template.slot.distr.weights.unnormalized
        norm = sum(slot_d)        
        writer.write("Slot Distribution:")
        for x in slot_d:
            writer.write(" %f" % (float(x)/float(norm)))
        writer.write("\n")
    if not args.data_frame:
        writer.write("{} for template {}".format(args.dist, template_idx))
    max_to_pick = min(N, len(tpred_sorted)) if args.N >= 0 else len(tpred_sorted)
    for (xrank,x) in enumerate(tpred_sorted[0:max_to_pick],1):
        if args.data_frame:
            writer.write("{}\t{}\t{}\t{}\t{}\n".format(vocab[x[0]].encode('utf-8'), x[1], template_idx, xrank, num_in_m[vocab[x[0]]] ))
        else:
            if args.print_vals:
                writer.write("({}, {}, {}){}".format(vocab[x[0]].encode('utf-8'), x[1], num_in_m[vocab[x[0]]], args.sep))
            else:
                writer.write("({}, , {}){}".format(vocab[x[0]].encode('utf-8'), num_in_m[vocab[x[0]]], args.sep))
    if not args.data_frame:
        writer.write("\n")
    if args.select_slots > 0 and not args.data_frame:
        slot_d = template.slot.distr.weights.unnormalized
        slot_norm = sum(slot_d)
        ssd = sorted(list(enumerate(slot_d)), key = lambda x : -x[1])
        max_to_pick = min(args.select_slots, len(slot_d))
        slot_vocab = rgs.vocabularies[1].words
        for (slot_id, slot_weight) in ssd[0:max_to_pick]:
            slot_item = rgs.slots[slot_id]
            choice = which(slot_item, 'arcs')
            if choice is None:
                print "listing %d doesn't have a set frame" % (slot_id)
                continue
            if choice.distr is None:
                print "listing %d doesn't have a set frame distribution" % (slot_id)
                continue
            lst = choice.distr.weights.residual
            if not len(lst) == len(slot_vocab):
                if choice.distr.items is None or len(choice.dist.items) == 0:
                    print "listing %d has a residual not of the size of vocab (%d vs %d), but not Frame.Distribution.items set either" % (slot_id, len(lst), len(slot_vocab))
                continue
                lst = zip(choice.distr.items, lst)
            else:
                lst = enumerate(lst)
            spred_sorted = list(lst)
            spred_sorted.sort(key = sorter)
            smax_to_pick = min(N, len(spred_sorted))
            writer.write("Slot %d (weight %f for template %d):%s" % (slot_id, slot_weight/slot_norm, template_idx, args.sep))
            for x in spred_sorted[0:smax_to_pick]:
                # if args.print_vals:
                #     writer.write("({}, {}, {}){}".format(slot_vocab[x[0]], x[1], num_in_m[slot_vocab[x[0]]], args.sep))
                # else:
                writer.write("({},){}".format(slot_vocab[x[0]],  args.sep))
            writer.write("\n")
        writer.write("\n")
