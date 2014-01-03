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
readBool "0" = False
readBool "1" = True

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

parseNCBITaxonomyDatabaseDumpCitations :: GenParser Char st [TaxDumpCitation]
parseNCBITaxonomyDatabaseDumpCitations = do
  citations <- many1 parseNCBITaxonomyDatabaseDumpCitation
  return $ citations

parseNCBITaxonomyDatabaseDumpDelNodes :: GenParser Char st [TaxDumpDelNode]
parseNCBITaxonomyDatabaseDumpDelNodes = do
  delNodes <- many1 parseNCBITaxonomyDatabaseDumpDelNode
  return $ delNodes
  
parseNCBITaxonomyDatabaseDumpDivisons :: GenParser Char st [TaxDumpDivision]
parseNCBITaxonomyDatabaseDumpDivisons = do
  divisions <- many1 parseNCBITaxonomyDatabaseDumpDivision
  return $ divisions

parseNCBITaxonomyDatabaseDumpGenCodes :: GenParser Char st [TaxDumpGenCode]
parseNCBITaxonomyDatabaseDumpGenCodes = do
  genCodes <- many1 parseNCBITaxonomyDatabaseDumpGenCode
  return $ genCodes

parseNCBITaxonomyDatabaseDumpMergedNodes :: GenParser Char st [TaxDumpMergedNode]
parseNCBITaxonomyDatabaseDumpMergedNodes = do
  mergedNodes <- many1 parseNCBITaxonomyDatabaseDumpMergedNode
  return $ mergedNodes

parseNCBITaxonomyDatabaseDumpNames :: GenParser Char st [TaxDumpName]
parseNCBITaxonomyDatabaseDumpNames = do
  names <- many1 parseNCBITaxonomyDatabaseDumpName
  return $ names

parseNCBITaxonomyDatabaseDumpNodes :: GenParser Char st [TaxDumpNode]
parseNCBITaxonomyDatabaseDumpNodes = do
  nodes <- many1 parseNCBITaxonomyDatabaseDumpNode
  return $ nodes

----------------------------

parseNCBITaxonomyDatabaseDumpCitation :: GenParser Char st TaxDumpCitation
parseNCBITaxonomyDatabaseDumpCitation = do
  citId <- many1 digit
  tab
  char ('|')
  tab 
  citKey <- many1 (noneOf "\t")
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
  url <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab
  text <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab
  taxIdList <- optionMaybe (many1 parseTaxIdList)
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpCitation (readInt citId) citKey (liftM readInt pubmedId) (liftM readInt medlineId) url text taxIdList

parseTaxIdList :: GenParser Char st Int
parseTaxIdList = do
  optional space
  taxId <- many1 digit
  return $ (readInt taxId)

parseNCBITaxonomyDatabaseDumpDelNode :: GenParser Char st TaxDumpDelNode
parseNCBITaxonomyDatabaseDumpDelNode = do
  delNode <- many1 digit
  space
  char ('|')
  char ('\n')
  return $ TaxDumpDelNode (readInt delNode)
  
parseNCBITaxonomyDatabaseDumpDivision :: GenParser Char st TaxDumpDivision
parseNCBITaxonomyDatabaseDumpDivision = do
  divisionId <- many1 digit
  tab
  char ('|')
  tab
  divisionCDE <- many1 upper
  tab
  char ('|')
  tab
  divisionName <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  comments <- optionMaybe (many1 (noneOf ("\t")))
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpDivision (readInt divisionId) divisionCDE divisionName comments 

parseNCBITaxonomyDatabaseDumpGenCode :: GenParser Char st TaxDumpGenCode
parseNCBITaxonomyDatabaseDumpGenCode = do
  geneticCodeId <- many1 digit 
  tab
  char ('|')
  abbreviation <- optionMaybe (many1 (noneOf ("\t")))
  tab
  char ('|')
  tab
  name <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  cde <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  starts <- many1 (noneOf ("\t"))
  char ('\n')
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
  char ('\n')
  return $ TaxDumpMergedNode (readInt oldTaxId) (readInt newTaxId)

parseNCBITaxonomyDatabaseDumpName :: GenParser Char st TaxDumpName
parseNCBITaxonomyDatabaseDumpName = do
  taxId <- many1 digit
  tab
  char ('|')
  tab
  nameTxt <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  uniqueName <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab
  nameClass <- many1 (noneOf ("\t"))
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
  rank <- many1 (noneOf "\t")
  tab
  char ('|')
  tab 
  emblCode <- optionMaybe (many1 (noneOf "\t"))
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
  comments <- many1 (noneOf "\t")
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpNode (readInt taxId) (readInt parentTaxId) rank emblCode divisionId (readBool inheritedDivFlag) geneticCodeId (readBool inheritedGCFlag) mitochondrialGeneticCodeId (readBool inheritedMGCFlag) (readBool genBankHiddenFlag) (readBool hiddenSubtreeRootFlag) comments



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
