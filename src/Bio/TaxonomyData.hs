-- | This module contains data structures for
--   taxonomy data

module Bio.TaxonomyData where

data SimpleTaxon = SimpleTaxon
  {
   -- node id in GenBank
   simpleTaxId :: Int,
   simpleScientificName :: String,
   -- parent node id in GenBank taxonomy database               
   simpleParentTaxId :: Int,
   -- rank of this node (superkingdom, kingdom, ...) 
   simpleRank :: Rank
  }
  deriving (Show, Read, Eq)

-- | NCBI Taxonomy database dump hierachichal data structure
-- as defined in ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
data NCBITaxDump = NCBITaxDump
  {
    taxCitations :: [TaxCitation],
    taxDelNodes :: [TaxDelNode],
    taxDivisions :: [TaxDivision],
    taxGenCodes :: [TaxGenCode],
    taxMergedNodes :: [TaxMergedNode],
    taxNames :: [TaxName],
    taxNodes :: [TaxNode]
  }
  deriving (Show, Read, Eq)

data TaxCitation = TaxCitation
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

data TaxDelNode = TaxDelNode
  {
   -- deleted node id
   delTaxId :: Int
  }
  deriving (Show, Read, Eq)

data TaxDivision = TaxDivision
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

data TaxGenCode = TaxGenCode
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

data TaxMergedNode = TaxMergedNode
  {
   -- id of nodes which has been merged
   oldTaxId :: Int,
   -- id of nodes which is result of merging
   newTaxId :: Int
  }
  deriving (Show, Read, Eq)

data TaxName = TaxName
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

-- | Taxonomic ranks: NCBI uses the uncommon Speciessubgroup 
data Rank = Norank | Form | Variety | Infraspecies | Subspecies | Speciessubgroup | Species | Speciesgroup | Superspecies | Series | Section | Subgenus | Genus | Subtribe | Tribe | Supertribe | Subfamily | Family | Superfamily | Parvorder | Infraorder | Suborder | Order | Superorder | Magnorder | Cohort | Legion | Parvclass | Infraclass | Subclass | Class | Superclass | Microphylum | Infraphylum | Subphylum | Phylum | Superphylum | Infrakingdom | Subkingdom | Kingdom | Superkingdom | Domain deriving (Eq, Ord, Show, Bounded, Enum)

readsRank :: [Char] -> [(Rank, [Char])]
instance Read Rank where
  readsPrec _ input = readsRank input

readsRank input -- = [(Domain x)| x <- reads input ]
   | input == "domain" = [(Domain,"")]
   | input == "superkingdom" = [(Superkingdom,"")]
   | input == "kingdom" = [(Kingdom,"")]
   | input == "subkingdom"  = [(Subkingdom,"")]
   | input == "infrakingdom" = [(Infrakingdom,"")]
   | input == "superphylum" = [(Superphylum,"")]
   | input == "phylum" = [(Phylum,"")]
   | input == "subphylum" = [(Subphylum,"")]
   | input == "infraphylum" = [(Infraphylum,"")]
   | input == "microphylum" = [(Microphylum,"")]
   | input == "superclass" = [(Superclass,"")]
   | input == "class" = [(Class,"")]
   | input == "subclass" = [(Subclass,"")]
   | input == "infraclass" = [(Infraclass,"")]
   | input == "parvclass " = [(Parvclass ,"")] 
   | input == "legion" = [(Legion,"")] 
   | input == "cohort" = [(Cohort,"")] 
   | input == "magnorder " = [(Magnorder ,"")] 
   | input == "superorder" = [(Superorder,"")] 
   | input == "order" = [(Order,"")]
   | input == "suborder" = [(Suborder,"")]
   | input == "infraorder" = [(Infraorder,"")] 
   | input == "parvorder" = [(Parvorder,"")] 
   | input == "superfamily" = [(Superfamily,"")]
   | input == "family" = [(Family,"")]
   | input == "subfamily" = [(Subfamily,"")]
   | input == "supertribe" = [(Supertribe,"")]
   | input == "tribe" = [(Tribe,"")] 
   | input == "subtribe" = [(Subtribe,"")] 
   | input == "genus" = [(Genus,"")]
   | input == "subgenus" = [(Subgenus,"")] 
   | input == "section" = [(Section,"")] 
   | input == "series" = [(Series,"")] 
   | input == "superspecies" = [(Superspecies,"")] 
   | input == "species group" = [(Speciesgroup,"")]
   | input == "species" = [(Species,"")]
   | input == "species subgroup" = [(Speciessubgroup,"")]
   | input == "subspecies" = [(Subspecies,"")] 
   | input == "infraspecies" = [(Infraspecies,"")]
   | input == "varietas" = [(Variety,"")]
   | input == "forma" = [(Form,"")]
   | input == "no rank" = [(Norank,"")]
   | otherwise = [(Norank,"")]  

data TaxNode = TaxNode
  {
   -- node id in GenBank
   taxId :: Int,
   -- parent node id in GenBank taxonomy database
   parentTaxId :: Int,
   -- rank of this node (superkingdom, kingdom, ...) 
   rank :: Rank,
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
