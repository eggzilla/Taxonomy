-- | This module contains data structures for
--   taxonomy data

module Bio.TaxonomyData where
import Data.Tree

-- | 
data SimpleTaxonomyEntry = SimpleTaxonomyEntry 
 {
  -- node id in GenBank
   entryId :: Int,
   -- rank of this node (superkingdom, kingdom, ...) 
   entryRank :: String,
   -- locus-name prefix; not unique
   entryName :: String,
   -- 
   entryGenomeAccession :: Maybe String
 }
 deriving (Show, Read, Eq)


-- | NCBI Taxonomy database dump hierachichal data structure
-- as defined in ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
data NCBITaxDump = NCBITaxDump
  {
    taxDumpCitations :: [TaxDumpCitation],
    taxDumpDelNodes :: [TaxDumpDelNode],
    taxDumpDivisions :: [TaxDumpDivision],
    taxDumpGenCodes :: [TaxDumpGenCode],
    taxDumpMergedNodes :: [TaxDumpMergedNode],
    taxDumpNames :: [TaxDumpName],
    taxDumpNodes :: [TaxDumpNode]
  }
  deriving (Show, Read, Eq)

data TaxDumpCitation = TaxDumpCitation
  {
   -- the unique id of citation
   citId :: Int,
   -- citation key
   citKey :: Maybe String,
   -- unique id in PubMed database (0 if not in PubMed)
   pubmedId :: Maybe Int,
   -- unique id in MedLine database (0 if not in MedLine)
   medlineId :: Maybe Int,
   -- URL associated with citation
   url :: Maybe String,
   -- any text (usually article name and authors)
   -- The following characters are escaped in this text by a backslash:
						-- newline (appear as "\n"),
						-- tab character ("\t"),
						-- double quotes ('\"'),
						-- backslash character ("\\").
   text :: Maybe String,
   -- list of node ids separated by a single space
   taxIdList :: Maybe [Int]
  }
  deriving (Show, Read, Eq)

data TaxDumpDelNode = TaxDumpDelNode
  {
   -- deleted node id
   delTaxId :: Int
  }
  deriving (Show, Read, Eq)

data TaxDumpDivision = TaxDumpDivision
  {
   -- taxonomy database division id
   divisionId :: Int,
   -- GenBank division code (three characters)
   divisionCDE :: String,
   -- e.g. BCT, PLN, VRT, MAM, PRI...
   divisonName :: String,
   divisionComments :: Maybe String
  }
  deriving (Show, Read, Eq)

data TaxDumpGenCode = TaxDumpGenCode
  {
   -- GenBank genetic code id
   geneticCodeId :: Int,
   -- genetic code name abbreviation
   abbreviation :: Maybe String,
   -- genetic code name
   geneCodeName :: String,
   -- translation table for this genetic code
   cde :: String,
   -- start codons for this genetic code
   starts :: String
  }
  deriving (Show, Read, Eq)

data TaxDumpMergedNode = TaxDumpMergedNode
  {
   -- id of nodes which has been merged
   oldTaxId :: Int,
   -- id of nodes which is result of merging
   newTaxId :: Int
  }
  deriving (Show, Read, Eq)

data TaxDumpName = TaxDumpName
  {
   -- the id of node associated with this name
   nameTaxId :: Int,
   -- name itself
   nameTxt :: String,
   -- the unique variant of this name if name not unique
   uniqueName :: Maybe String,
   -- (synonym, common name, ...)
   nameClass :: String
  }
  deriving (Show, Read, Eq)

data TaxDumpNode = TaxDumpNode
  {
   -- node id in GenBank
   taxId :: Int,
   -- parent node id in GenBank taxonomy database
   parentTaxId :: Int,
   -- rank of this node (superkingdom, kingdom, ...) 
   rank :: String,
   -- locus-name prefix; not unique
   emblCode :: Maybe String,
   -- see division.dmp file
   nodeDivisionId :: String,
   -- 1 if node inherits division from parent
   inheritedDivFlag :: Bool,
   -- see gencode.dmp file
   nodeGeneticCodeId :: String,
   -- 1 if node inherits genetic code from parent
   inheritedGCFlag :: Bool,
   -- see gencode.dmp file 
   mitochondrialGeneticCodeId :: String,
   -- 1 if node inherits mitochondrial gencode from parent
   inheritedMGCFlag :: Bool,
   -- 1 if name is suppressed in GenBank entry lineage
   genBankHiddenFlag :: Bool,
   -- 1 if this subtree has no sequence data yet
   hiddenSubtreeRootFlag :: Bool,
   -- free-text comments and citations
   nodeComments :: Maybe String
  }
  deriving (Show, Read, Eq)
