package com.matthicks.genomesequencing

import java.io.File

import com.outr.scribe.Logging

import scala.annotation.tailrec
import scala.io.Source

/**
  * Genome Assembly as Shortest Superstring
  * http://rosalind.info/problems/long/
  *
  * Uses asynchronous operations to improve performance.
  *
  * @author Matt Hicks
  */
object GenomeAssemblyAsShortestSuperstring extends Logging {
  def main(args: Array[String]): Unit = {
    val strands = load(new File("rosalind_long.fasta"))
    val shortest = superstringFor(strands)
    logger.info(s"Shortest Superstring: $shortest (${shortest.length})")
  }

  /**
    * Loads a FASTA format file returning a collection of DNA strands.
    *
    * @param file the .fas or .fasta file to load
    * @return set of DNA strands
    */
  def load(file: File): Set[String] = {
    var strands = Set.empty[String]
    val b = new StringBuilder
    val source = Source.fromFile(file)
    try {
      source.getLines().foreach { line =>
        if (line.startsWith(">")) {
          if (b.nonEmpty) {
            strands += b.toString()
            b.clear()
          }
        } else {
          b.append(line.trim)
        }
      }
      val s = b.toString().trim
      if (s.nonEmpty) {
        strands += s
      }
      strands
    } finally {
      source.close()
    }
  }

  /**
    * Analyzes the strands and builds the shortest superstring from it.
    *
    * @param strands represents a collection of DNA strands
    * @return the shortest superstring combining all strands
    */
  def superstringFor(strands: Set[String]): String = {
    val dnaMatches = strands.map { strand =>
      val set = strands - strand
      strand -> set.flatMap(findOverlap(strand, _))
    }.toMap

    val singleConnection = dnaMatches.find {
      case (strand, matches) => matches.tail.isEmpty
    }
    singleConnection match {
      case Some((strand, matches)) => {
        val dnaMatch = matches.head
        val updated = strands - strand - dnaMatch.strand
        val next = dnaMatches(dnaMatch.strand)
        analyze(dnaMatches, dnaMatch.combined, next, updated)
      }
      case None => throw new RuntimeException("No known endpoint!")   // TODO: evaluate if support for this is necessary
    }
  }

  /**
    * Analyzes the paths to build a superstring and returns the shortest result.
    *
    * WARNING: This uses non-tail recursion and adds levels to the stack equal to the number of strands being processed.
    * As a result, this can potentially result in a stack overflow if too many strands are being analyzed.
    *
    * @param dnaMatches the pre-built dna matches to determine paths
    * @param dna the currently combined DNA
    * @param next the possible matches of strands that can be connected
    * @param remaining the unused strands that must be connected before the path can be resolved
    * @return shortest superstring as a String
    */
  private def analyze(dnaMatches: Map[String, Set[DNAMatch]],
                      dna: String,
                      next: Set[DNAMatch],
                      remaining: Set[String]): String = {
    if (remaining.isEmpty) {
      dna
    } else {
      val potential: Set[DNAMatch] = next.collect {
        case o if remaining.contains(o.strand) => findOverlap(dna, o.strand)
      }.flatten
      val results: Set[String] = potential.map { dnaMatch =>
        analyze(dnaMatches, dnaMatch.combined, dnaMatches(dnaMatch.strand), remaining - dnaMatch.strand)
      }
      results.foldLeft("")((shortest, dna) => if (shortest.isEmpty || dna.length < shortest.length) dna else shortest)
    }
  }

  /**
    * Finds the overlap of a combined DNA with a strand. Only returns a match if the overlap is 50% or greater.
    *
    * @param dna the combined DNA to match the strand to
    * @param strand the DNA strand to assemble
    * @param offset the offset, since this is a recursive method this will increment until it reaches 50% or
    *               finds a match
    * @return DNAMatch if one is available
    */
  @tailrec
  final def findOverlap(dna: String, strand: String, offset: Int = 1): Option[DNAMatch] = {
    if (offset > (math.max(dna.length, strand.length) / 2) - 1) {
      None
    } else if (offset == 1 && dna.indexOf(strand) != -1) {     // Matched entire strand
      Some(new DNAMatch(dna, strand))
    } else if (dna.startsWith(strand.substring(offset))) {
      Some(new DNAMatch(s"${strand.substring(0, offset)}$dna", strand))
    } else if (dna.endsWith(strand.substring(0, strand.length - (offset + 1)))) {
      Some(new DNAMatch(s"$dna${strand.substring(strand.length - (offset + 1))}", strand))
    } else  {
      findOverlap(dna, strand, offset + 1)
    }
  }

  /**
    * Placeholder internally for an overlapping match.
    *
    * @param combined the combination of the strand and the dna
    * @param strand the strand used to derive combined
    */
  class DNAMatch(val combined: String, val strand: String)
}