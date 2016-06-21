package spec

import java.io.File

import com.matthicks.genomesequencing.GenomeAssemblyAsShortestSuperstring
import org.scalatest.{Matchers, WordSpec}

class SuperstringSpec extends WordSpec with Matchers {
  "Superstring" when {
    "using the short file" should {
      var strands = Set.empty[String]
      "load the file" in {
        strands = GenomeAssemblyAsShortestSuperstring.load(new File("rosalind_short.fasta"))
      }
      "verify the correct number of strands" in {
        strands.size should equal(4)
      }
      "verify the shortest superstring is correct" in {
        val shortest = GenomeAssemblyAsShortestSuperstring.superstringFor(strands)
        shortest should equal("ATTAGACCTGCCGGAATAC")
      }
    }
    "using the long file" should {
      var strands = Set.empty[String]
      "load the file" in {
        strands = GenomeAssemblyAsShortestSuperstring.load(new File("rosalind_long.fasta"))
      }
      "verify the correct number of strands" in {
        strands.size should equal(50)
      }
      "verify the shortest superstring is the correct length" in {
        val shortest = GenomeAssemblyAsShortestSuperstring.superstringFor(strands)
        shortest.length should equal(19914)
      }
    }
  }
}
