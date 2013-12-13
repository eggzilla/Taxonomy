-- | This module contains data structures for
--   taxonomy data

module Bio.TaxonomyData where
    
-- | RNAplex output consists of a set of interactions 
data Taxonomy = Taxonomy
  {
    taxData :: String
  }
  deriving (Show, Eq)

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
   citId :: String,
   citKey :: String,
   pubmedId :: String,
   medlineId :: String,
   url :: String,
   text :: String,
   taxidList :: [String]
  }
  deriving (Show, Eq)

data TaxDumpDelNode = Taxonomy
  {
   taxID :: String
  }
  deriving (Show, Eq)

data TaxDumpDivision = Taxonomy
  {
   divisionId :: String,
   divisionCDE :: String,
   divisonName :: String,
   comments :: String
  }
  deriving (Show, Eq)

data TaxDumpGencode = Taxonomy
  {
   geneticCodeId :: String,
   abbreviation :: String,
   name :: String,
   cde :: String,
   starts :: String
  }
  deriving (Show, Eq)

data TaxDumpMergedNode = Taxonomy
  {
   oldTaxId :: String,
   newTaxId :: String
  }
  deriving (Show, Eq)

data TaxDumpName = Taxonomy
  {
   taxID :: String,
   nameTxt :: String,
   uniqueName :: String,
   nameClass :: String
  }
  deriving (Show, Eq)

data TaxDumpNode = Taxonomy
  {
   taxId :: String,
   parentTaxId :: String,
   rank :: String,
   emblCode :: String,
   divisionId :: String,
   inheritiedDivFlag :: String,
   geneticCodeId :: String,
   inheritedGCFlag :: String,
   mitochondrialGeneticCodeId :: String,
   inheritedMGCFlag :: String,
   genBankHiddenFlag :: String,
   hiddenSubtreeRootFlag :: String,
   comments :: String
  }
  deriving (Show, Eq)
