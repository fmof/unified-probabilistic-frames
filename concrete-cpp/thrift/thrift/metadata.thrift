/*
 * Copyright 2012-2014 Johns Hopkins University HLTCOE. All rights reserved.
 * This software is released under the 2-clause BSD license.
 * See LICENSE in the project root directory.
 */

namespace java edu.jhu.hlt.concrete
namespace py concrete.metadata
namespace cpp concrete
#@namespace scala edu.jhu.hlt.miser

include "uuid.thrift"
include "twitter.thrift"
include "email.thrift"
include "nitf.thrift"

/**
 * A struct that holds UUIDs for all theories that a particular
 * annotation was based upon (and presumably requires).
 *
 * Producers of TheoryDependencies should list all stages that they
 * used in constructing their particular annotation. They do not, 
 * however, need to explicitly label *each* stage; they can label
 * only the immediate stage before them.
 * 
 * Examples:
 * 
 * If you are producing a Tokenization, and only used the
 * SentenceSegmentation in order to produce that Tokenization, list
 * only the single SentenceSegmentation UUID in sentenceTheoryList.
 *
 * In this example, even though the SentenceSegmentation will have
 * a dependency on some SectionSegmentation, it is not necessary
 * for the Tokenization to list the SectionSegmentation UUID as a
 * dependency. 
 *
 * If you are a producer of EntityMentions, and you use two
 * POSTokenTagging and one NERTokenTagging objects, add the UUIDs for
 * the POSTokenTagging objects to posTagTheoryList, and the UUID of
 * the NER TokenTagging to the nerTagTheoryList.
 *
 * In this example, because multiple annotations influenced the 
 * new annotation, they should all be listed as dependencies.
 */
struct TheoryDependencies {
  1: optional list<uuid.UUID> sectionTheoryList
  2: optional list<uuid.UUID> sentenceTheoryList
  3: optional list<uuid.UUID> tokenizationTheoryList
  4: optional list<uuid.UUID> posTagTheoryList
  5: optional list<uuid.UUID> nerTagTheoryList
  6: optional list<uuid.UUID> lemmaTheoryList
  7: optional list<uuid.UUID> langIdTheoryList
  8: optional list<uuid.UUID> parseTheoryList
  9: optional list<uuid.UUID> dependencyParseTheoryList
  10: optional list<uuid.UUID> tokenAnnotationTheoryList
  11: optional list<uuid.UUID> entityMentionSetTheoryList
  12: optional list<uuid.UUID> entitySetTheoryList
  13: optional list<uuid.UUID> situationMentionSetTheoryList
  14: optional list<uuid.UUID> situationSetTheoryList
  15: optional list<uuid.UUID> communicationsList
}

//===========================================================================
// Metadata
//===========================================================================

/** 
 * Analytic-specific information about an attribute or edge. Digests
 * are used to combine information from multiple sources to generate a
 * unified value. The digests generated by an analytic will only ever
 * be used by that same analytic, so analytics can feel free to encode
 * information in whatever way is convenient. 
 */
struct Digest {
  /** 
   * The following fields define various ways you can store the
   * digest data (for convenience). If none of these meets your
   * needs, then serialize the digest to a byte sequence and store it
   * in bytesValue. 
   */
  1: optional binary bytesValue
  2: optional i64 int64Value
  3: optional double doubleValue
  4: optional string stringValue
  5: optional list<i64> int64List
  6: optional list<double> doubleList
  7: optional list<string> stringList
}

/** 
 * Metadata associated with an annotation or a set of annotations,
 * that identifies where those annotations came from. 
 */
struct AnnotationMetadata {
  /** 
   * The name of the tool that generated this annotation. 
   */
  1: required string tool

  /** 
   * The time at which this annotation was generated (in unix time
   * UTC -- i.e., seconds since January 1, 1970). 
   */
  2: required i64 timestamp

  /** 
   * A Digest, carrying over any information the annotation metadata
   * wishes to carry over.
   */
  4: optional Digest digest

  /**
   * The theories that supported this annotation. 
   * 
   * An empty field indicates that the theory has no 
   * dependencies (e.g., an ingester).
   */
  5: optional TheoryDependencies dependencies
  
  /**
   * An integer that represents a ranking for systems
   * that output k-best lists. 
   * 
   * For systems that do not output k-best lists, 
   * the default value (1) should suffice.
   */
  6: required i32 kBest = 1
}

/**
 * Metadata specific to a particular Communication object.
 * This might include corpus-specific metadata (from the Twitter API),
 * attributes associated with the Communication (the author),
 * or other information about the Communication.
 */
struct CommunicationMetadata {
  /** 
   * Extra information for communications where kind==TWEET:
   * Information about this tweet that is provided by the Twitter
   * API.  For information about the Twitter API, see:
   * https://dev.twitter.com/docs/platform-objects
   */
  1: optional twitter.TweetInfo tweetInfo

  /**
   * Extra information for communications where kind==EMAIL
   */
  2: optional email.EmailCommunicationInfo emailInfo

  /**
   * Extra information that may come from the NITF
   * (News Industry Text Format) schema. See 'nitf.thrift'.
   */

  3: optional nitf.NITFInfo nitfInfo
}