-- | This module contains data structures for
--   taxonomy data

{-# LANGUAGE FlexibleInstances #-}

module Bio.TaxonomyData where
import Prelude
--import qualified Data.ByteString as B
import qualified Data.Aeson as A
import qualified Data.Vector as V
import Data.Graph.Inductive
import qualified Data.Text as T
import qualified Data.Text.Lazy as TL   
--import qualified Data.Text.Encoding

-- | SimpleTaxon only contains the most relevant fields of a taxonomy entry.
--   For all annotaded fields use the Taxon datatype and its associated functions
data SimpleTaxon = SimpleTaxon
  {
   -- node id in GenBank
   simpleTaxId :: Int,
   simpleScientificName :: TL.Text,
   -- parent node id in GenBank taxonomy database               
   simpleParentTaxId :: Int,
   -- rank of this node (superkingdom, kingdom, ...) 
   simpleRank :: Rank
  }
  deriving (Show, Read, Eq)

-- | Datastructure for tree comparisons
data CompareTaxon = CompareTaxon
  {
   compareScientificName :: TL.Text,
   compareRank :: Rank,
   -- number indicating in which trees, 
   inTree :: [Int]
  }
  deriving (Show, Read, Eq)

-- | Data structure for Entrez taxonomy fetch result
data Taxon = Taxon
  {  taxonTaxId :: Int
  ,  taxonScientificName :: String
  ,  taxonParentTaxId :: Int
  ,  taxonRank :: Rank
  ,  division :: String
  ,  geneticCode :: TaxGenCode
  ,  mitoGeneticCode :: TaxGenCode
  ,  lineage :: String
  ,  lineageEx :: [LineageTaxon]
  ,  createDate :: String
  ,  updateDate :: String
  ,  pubDate :: String
  } deriving (Show, Eq)


data TaxonName = TaxonName
  {  classCDE :: String
  ,  dispName :: String
  } deriving (Show, Eq)

-- | Lineage Taxons denote all parent Taxonomy nodes of a node retrieved by Entrez fetch
data LineageTaxon = LineageTaxon
  {  lineageTaxId :: Int
  ,  lineageScienticName :: String
  ,  lineageRank :: Rank}
  deriving (Show, Eq)
           
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

-- | Datastructure for entries of Taxonomy database dump citations file
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

-- | Datastructure for entries of Taxonomy database dump deleted nodes file
data TaxDelNode = TaxDelNode
  {
   -- deleted node id
   delTaxId :: Int
  }
  deriving (Show, Read, Eq)

-- | Datastructure for entries of Taxonomy database dump division file
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

-- | Datastructure for entries of Taxonomy database dump gencode file
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

-- | Datastructure for entries of Taxonomy database dump mergednodes file
data TaxMergedNode = TaxMergedNode
  {
   -- id of nodes which has been merged
   oldTaxId :: Int,
   -- id of nodes which is result of merging
   newTaxId :: Int
  }
  deriving (Show, Read, Eq)

-- | Datastructure for entries of Taxonomy database dump names file
data TaxName = TaxName
  {
   -- the id of node associated with this name
   nameTaxId :: Int,
   -- name itself
   nameTxt :: TL.Text,
   -- the unique variant of this name if name not unique
   uniqueName :: TL.Text,
   -- (synonym, common name, ...)
   nameClass :: TL.Text
  }
  deriving (Show, Read, Eq)

-- | Taxonomic ranks: NCBI uses the uncommon Speciessubgroup 
data Rank = Norank | Form | Variety | Infraspecies | Subspecies | Speciessubgroup | Species | Speciesgroup | Superspecies | Series | Section | Subgenus | Genus | Subtribe | Tribe | Supertribe | Subfamily | Family | Superfamily | Parvorder | Infraorder | Suborder | Order | Superorder | Magnorder | Cohort | Legion | Parvclass | Infraclass | Subclass | Class | Superclass | Microphylum | Infraphylum | Subphylum | Phylum | Superphylum | Infrakingdom | Subkingdom | Kingdom | Superkingdom | Domain deriving (Eq, Ord, Show, Bounded, Enum)

readsRank :: String -> [(Rank, String)]
instance Read Rank where
  readsPrec _ = readsRank 

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

-- | Datastructure for entries of Taxonomy database dump nodes file
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

-- | Simple Gene2Accession table 
data SimpleGene2Accession = SimpleGene2Accession
  { simpleTaxIdEntry :: Int,
    simpleGenomicNucleotideAccessionVersion :: String
  } deriving (Show, Eq, Read) 

-- | Datastructure for Gene2Accession table
data Gene2Accession = Gene2Accession
  { taxIdEntry :: Int,
    geneID :: Int,
    status :: String,
    rnaNucleotideAccessionVersion :: String,
    rnaNucleotideGi :: String,
    proteinAccessionVersion :: String,
    proteinGi :: String,
    genomicNucleotideAccessionVersion :: String,
    genomicNucleotideGi :: String,
    startPositionOnTheGenomicAccession :: String,
    endPositionOnTheGenomicAccession ::  String,
    orientation :: String,
    assembly :: String,
    maturePeptideAccessionVersion :: String,
    maturePeptideGi :: String
  } deriving (Show, Eq, Read)  

instance A.ToJSON (Gr SimpleTaxon Double) where
  toJSON inputGraph = simpleTaxonJSONValue inputGraph 1

simpleTaxonJSONValue :: Gr SimpleTaxon Double -> Node -> A.Value
simpleTaxonJSONValue inputGraph node = jsonValue
  where jsonValue = A.object [currentScientificName,T.pack "children" A..= children]
        childNodes = suc inputGraph node
        currentLabel = lab inputGraph node
        currentScientificName = T.pack "name" A..= maybe (T.pack "notFound") (TL.toStrict  . simpleScientificName) currentLabel
        children = A.Array (V.fromList (map (simpleTaxonJSONValue inputGraph) childNodes))
        --jsonValue = A.object [currentScientificName,currentId,currentRank,(T.pack "children") A..= children]
        --currentId = (T.pack "id") A..= (maybe (T.pack "notFound") (\a -> T.pack (show (simpleTaxId a))) currentLabel)
        --currentRank = (T.pack "rank") A..= (maybe (T.pack "notFound") (\a -> T.pack (show (simpleRank a))) currentLabel)
