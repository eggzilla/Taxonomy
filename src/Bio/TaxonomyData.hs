-- | This module contains data structures for
--   taxonomy data

module Bio.TaxonomyData where
    
-- | 
data Taxonomy = Taxonomy
  {
    taxData :: String
  }
  deriving (Show, Eq)

-- | NCBI Taxonomy database dump hierachichal data structure
-- as defined in ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
data NCBITaxDump = Taxonomy
  {
    TaxDumpCitations :: [TaxDumpCitation],
    TaxDumpDelNodes :: [TaxDumpDelNode],
    TaxDumpDivisions :: [TaxDumpDivison],
    TaxDumpGencodes :: [TaxDumpGencode],
    TaxDumpMergedNodes :: [TaxDumpMergedNode],
    TaxDumpNames :: [TaxDumpName],
    TaxDumpNodes :: [TaxDumpNode]
  }
  deriving (Show, Eq)

data TaxDumpCitation = Taxonomy
  {
   -- the unique id of citation
   citId :: String,
   -- citation key
   citKey :: String,
   -- unique id in PubMed database (0 if not in PubMed)
   pubmedId :: Maybe String,
   -- unique id in MedLine database (0 if not in MedLine)
   medlineId :: Maybe String,
   -- URL associated with citation
   url :: String,
   -- any text (usually article name and authors)
   -- The following characters are escaped in this text by a backslash:
						-- newline (appear as "\n"),
						-- tab character ("\t"),
						-- double quotes ('\"'),
						-- backslash character ("\\").
   text :: String,
   -- list of node ids separated by a single space
   taxidList :: [String]
  }
  deriving (Show, Eq)

data TaxDumpDelNode = Taxonomy
  {
   -- deleted node id
   taxID :: Int
  }
  deriving (Show, Eq)

data TaxDumpDivision = Taxonomy
  {
   -- taxonomy database division id
   divisionId :: String,
   -- GenBank division code (three characters)
   divisionCDE :: String,
   -- e.g. BCT, PLN, VRT, MAM, PRI...
   divisonName :: String,
   comments :: String
  }
  deriving (Show, Eq)

data TaxDumpGencode = Taxonomy
  {
   -- GenBank genetic code id
   geneticCodeId :: String,
   -- genetic code name abbreviation
   abbreviation :: String,
   -- genetic code name
   name :: String,
   -- translation table for this genetic code
   cde :: String,
   -- start codons for this genetic code
   starts :: String
  }
  deriving (Show, Eq)

data TaxDumpMergedNode = Taxonomy
  {
   -- id of nodes which has been merged
   oldTaxId :: Int,
   -- id of nodes which is result of merging
   newTaxId :: Int
  }
  deriving (Show, Eq)

data TaxDumpName = Taxonomy
  {
   -- the id of node associated with this name
   taxID :: Int,
   -- name itself
   nameTxt :: String,
   -- the unique variant of this name if name not unique
   uniqueName :: String,
   -- (synonym, common name, ...)
   nameClass :: String
  }
  deriving (Show, Eq)

data TaxDumpNode = Taxonomy
  {
   -- node id in GenBank
   taxId :: Int,
   -- parent node id in GenBank taxonomy database
   parentTaxId :: Int,
   -- rank of this node (superkingdom, kingdom, ...) 
   rank :: String,
   -- locus-name prefix; not unique
   emblCode :: String,
   -- see division.dmp file
   divisionId :: String,
   -- 1 if node inherits division from parent
   inheritiedDivFlag :: Bool,
   -- see gencode.dmp file
   geneticCodeId :: String,
   -- 1 if node inherits genetic code from parent
   inheritedGCFlag :: Bool,
   -- see gencode.dmp file 
   mitochondrialGeneticCodeId :: String,
   -- 1 if node inherits mitochondrial gencode from parent
   inheritedMGCFlag :: Bool,
   -- 1 if name is suppressed in GenBank entry lineage
   genBankHiddenFlag :: Boolean,
   -- 1 if this subtree has no sequence data yet
   hiddenSubtreeRootFlag :: Bool,
   -- free-text comments and citations
   comments :: String
  }
  deriving (Show, Eq)
