-- | Parse and process taxonomy data

module Bio.Taxonomy (                      
                       module Bio.TaxonomyData
                      ) where

import Bio.TaxonomyData
import Text.ParserCombinators.Parsec
import Text.ParserCombinators.Parsec.Token
import Text.ParserCombinators.Parsec.Language (emptyDef)    
import Control.Monad

--Auxiliary functions

readDouble :: String -> Double
readDouble = read              

readInt :: String -> Int
readInt = read

readBool :: String -> Bool
readBool 0 = False
readBool 1 = True

--readNCBITaxonomyDatabaseDump
--taxdump consists of several files

-- | Parse the input as NCBITaxDump datatype
--parseNCBITaxonomyDatabaseDump :: GenParser Char st NCBITaxDump
--parseNCBITaxonomyDatabaseDump = do
--  citations <- parseNCBITaxDumpCitations
--  delNodes <- parse 
--  divisons <-
--  genCodes <-
--  mergedNodes <-
--  names <-
--  nodes <- 
--  eof  
--  return $ NCBITaxDump citations delNodes 

parseNCBITaxonomyDatabaseDumpCitations :: GenParser Char st TaxDumpCitations
parseNCBITaxonomyDatabaseDumpCitations = do
  citations <- many1 parseNCBITaxDumpCitations
  return $ TaxDumpCitations citations

parseNCBITaxonomyDatabaseDumpDelNodes :: GenParser Char st TaxDumpDelNodes
parseNCBITaxonomyDatabaseDumpDelNodes = do
  delNodes <- many1 parseNCBITaxDumpCitations
  return $ TaxDelNodes delNodes
  
parseNCBITaxonomyDatabaseDumpDivisons :: GenParser Char st TaxDumpDivisions
parseNCBITaxonomyDatabaseDumpDivisons = do
  divisions <- parseNCBITaxDumpCitations
  return $ TaxDumpDivisions divisions

parseNCBITaxonomyDatabaseDumpGenCodes :: GenParser Char st TaxDumpGenCodes
parseNCBITaxonomyDatabaseDumpGenCodes = do
  genCodes <- parseNCBITaxDumpCitations
  return $ TaxDumpGenCodes genCodes

parseNCBITaxonomyDatabaseDumpMergedNodes :: GenParser Char st TaxDumpMergedNodes
parseNCBITaxonomyDatabaseDumpMergedNodes = do
  mergedNodes <- parseNCBITaxDumpCitations
  return $ TaxDumpMergedNodes mergedNodes

parseNCBITaxonomyDatabaseDumpNames :: GenParser Char st TaxDumpNames
parseNCBITaxonomyDatabaseDumpNames = do
  names <- parseNCBITaxDumpCitations
  return $ TaxDumpNames names

parseNCBITaxonomyDatabaseDumpNodes :: GenParser Char st TaxDumpNodes
parseNCBITaxonomyDatabaseDumpNodes = do
  nodes <- parseNCBITaxDumpCitations
  return $ TaxDumpNodes nodes

----------------------------

parseNCBITaxonomyDatabaseDumpCitation :: GenParser Char st TaxDumpCitation
parseNCBITaxonomyDatabaseDumpCitation = do
  citId <- many1 digit
  tab
  char ('|')
  tab 
  citKey <- many1 (noneOf '\t')
  tab
  char ('|')
  tab 
  pubmedId <- optionMaybe (many1 digit)
  tab
  char ('|')
  tab  
  medlineId <- optionMaybe (many1 digit)
  tab
  char ('|')
  tab 
  url <- optionMaybe (many1 (noneOf ("\t")))
  tab
  char ('|')
  tab
  text <- optionMaybe (many1 (noneOf '\t'))
  tab
  char ('|')
  tab
  taxIdList optionMaybe (many1 parseTaxIdList)
  tab
  char ('|')
  eol
  return $ TaxDumpCitation (readInt citId) citKey (liftM readInt pubmedId) (liftM readInt medlineId) url text taxIdList

parseTaxIdList :: GenParser Char st Int
  optional space
  taxId many1 digits
  noneof ('\n')
  return (readInt taxId)

parseNCBITaxonomyDatabaseDumpDelNode :: GenParser Char st TaxDumpDelNode
parseNCBITaxonomyDatabaseDumpDelNode = do
  delNode <- many1 digit
  space
  char ('|')
  eol
  return $ (readInt delNode)
  
parseNCBITaxonomyDatabaseDumpDivison :: GenParser Char st TaxDumpDivision
parseNCBITaxonomyDatabaseDumpDivison = do
  divisionId <- many1 digit
  tab
  char ('|')
  tab
  divisionCDE <- many1 upper
  tab
  char ('|')
  tab
  divisionName <- many1 char
  tab
  char ('|')
  tab
  comments <- optionMaybe (many1 char)
  tab
  char ('|')
  eol
  return $ TaxDumpDivision (readInt divisionId) divisionCDE divisionName comments 

parseNCBITaxonomyDatabaseDumpGenCode :: GenParser Char st TaxDumpGenCode
parseNCBITaxonomyDatabaseDumpGenCode = do
  geneticCodeId <- many1 digit 
  tab
  char ('|')
  abbreviation <- optionMaybe (many1 char)
  tab
  char ('|')
  tab
  name <- many1 char
  tab
  char ('|')
  tab
  cde <- many1 char
  tab
  char ('|')
  tab
  starts <- many1 char
  eol
  return $ TaxDumpGenCode (readInt geneticCodeId) abbreviation name cde starts

parseNCBITaxonomyDatabaseDumpMergedNode :: GenParser Char st TaxDumpMergedNode
parseNCBITaxonomyDatabaseDumpMergedNode = do
  oldTaxId <- many1 digit
  tab
  char ('|')
  tab  
  newTaxId <- many1 digit
  tab
  char ('|')
  eol
  return $ TaxDumpMergedNode (readInt oldTaxId) (readInt newTaxId)

parseNCBITaxonomyDatabaseDumpName :: GenParser Char st TaxDumpName
parseNCBITaxonomyDatabaseDumpName = do
  taxId <- many1 digit
  tab
  char ('|')
  tab
  nameTxt <- many1 char
  tab
  char ('|')
  tab
  uniqueName <- optionMaybe (many1 noneOf ('\t'))
  tab
  char ('|')
  tab
  nameClass <- many1 noneOf ('\t')
  return $ TaxDumpName (readInt taxId) nameTxt uniqueName nameClass

parseNCBITaxonomyDatabaseDumpNode :: GenParser Char st TaxDumpNode
parseNCBITaxonomyDatabaseDumpNode = do
  taxId <- many1 digit
  tab
  char ('|') 
  tab
  parentTaxId <- many1 digit
  tab
  char ('|')
  tab
  rank <- many1 (noneOf ('\t'))
  tab
  char ('|')
  tab 
  emblCode <- optionMaybe (noneOf ('\t'))
  tab
  char ('|')
  tab
  divisionId <- many1 digit
  tab
  char ('|')
  tab
  inheritedDivFlag <- many1 digit
  tab
  char ('|')
  tab 
  geneticCodeId <- many1 digit
  tab
  char ('|')
  tab
  inheritedGCFlag <- many1 digit
  tab
  char ('|')
  tab
  mitochondrialGeneticCodeId <- many1 digit
  tab
  char ('|')
  tab
  inheritedMGCFlag <- many1 digit
  tab
  char ('|')
  tab
  genBankHiddenFlag <- many1 digit
  tab
  char ('|')
  tab
  hiddenSubtreeRootFlag <- many1 digit 
  tab
  char ('|')
  tab
  comments <- many1 (noneOf ('\t'))
  tab
  char ('|')
  eol
  return $ TaxDumpNode (readInt taxId) (readInt parentTaxId) rank emblCode divisionId (readBool inheritedDivFlag) geneticCodeId (readBool inheritedGCFlag) mitochondrialGeneticCodeId (readBool inheritedMGCFlag) (readBool genBankHiddenFlag) (readBool hiddenSubtreeRootFlag) comments

readBool :: String -> Bool
readBool 0 = False
readBool 1 = True

--parent
--returns the parent or Nothing

--childen
--returns direct children of this node

--allchildren
--returns all children of this node

--sequenced children
--returns all nodes with avialable genomes

--parent_rank
--returns the parent with the specified rank or Nothing
