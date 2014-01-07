-- | Parse and process taxonomy data

module Bio.Taxonomy (                      
                       module Bio.TaxonomyData,
                       parseNCBITaxDumpCitations,
                       readNCBITaxDumpCitations,
                       parseNCBITaxDumpDelNodes,
                       readNCBITaxDumpDelNodes,
                       parseNCBITaxDumpDivisons,
                       readNCBITaxDumpDivisons,
                       parseNCBITaxDumpGenCodes,
                       readNCBITaxDumpGenCodes,
                       parseNCBITaxDumpMergedNodes,
                       readNCBITaxDumpMergedNodes,
                       parseNCBITaxDumpNames,
                       readNCBITaxDumpNames,
                       parseNCBITaxDumpNodes,
                       readNCBITaxDumpNodes
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

genParserTaxIdList :: GenParser Char st Int
genParserTaxIdList = do
  optional (char ' ')
  taxId <- many1 digit
  optional (char ' ')
  return $ (readInt taxId)

genParserTaxURL :: GenParser Char st String
genParserTaxURL = do
  --url1 <- (noneOf "\t") `sepBy` (optional (string "\tN"))
  url1 <- many1 (noneOf "\t")
  
  -- some URL fields contain \t characters
  --string "\t"
  --notFollowedBy char '|'
  --url2 <- (many1 (noneOf ("\t")))
  return $ url1 -- ++ url2)


-- | parse NCBITaxDumpCitations from input string
parseNCBITaxDumpCitations input = parse genParserNCBITaxDumpCitations "parseTaxDumpCitations" input

-- | parse NCBITaxDumpCitations from input filePath                      
readNCBITaxDumpCitations :: String -> IO (Either ParseError [TaxDumpCitation])  
readNCBITaxDumpCitations filePath = parseFromFile genParserNCBITaxDumpCitations filePath

-- | parse NCBITaxDumpDelNodes from input string
parseNCBITaxDumpDelNodes input = parse genParserNCBITaxDumpDelNodes "parseTaxDumpDelNodes" input

-- | parse NCBITaxDumpDelNodes from input filePath                      
readNCBITaxDumpDelNodes :: String -> IO (Either ParseError [TaxDumpDelNode])  
readNCBITaxDumpDelNodes filePath = parseFromFile genParserNCBITaxDumpDelNodes filePath

-- | parse NCBITaxDumpDivisons from input string
parseNCBITaxDumpDivisons input = parse genParserNCBITaxDumpDivisons "parseTaxDumpDivisons" input

-- | parse NCBITaxDumpDivisons from input filePath                      
readNCBITaxDumpDivisons :: String -> IO (Either ParseError [TaxDumpDivision])  
readNCBITaxDumpDivisons filePath = parseFromFile genParserNCBITaxDumpDivisons filePath

-- | parse NCBITaxDumpGenCodes from input string
parseNCBITaxDumpGenCodes input = parse genParserNCBITaxDumpGenCodes "parseTaxDumpGenCodes" input

-- | parse NCBITaxDumpGenCodes from input filePath                      
readNCBITaxDumpGenCodes :: String -> IO (Either ParseError [TaxDumpGenCode])  
readNCBITaxDumpGenCodes filePath = parseFromFile genParserNCBITaxDumpGenCodes filePath

-- | parse NCBITaxDumpMergedNodes from input string
parseNCBITaxDumpMergedNodes input = parse genParserNCBITaxDumpMergedNodes "parseTaxDumpMergedNodes" input

-- | parse NCBITaxDumpMergedNodes from input filePath                      
readNCBITaxDumpMergedNodes :: String -> IO (Either ParseError [TaxDumpMergedNode])  
readNCBITaxDumpMergedNodes filePath = parseFromFile genParserNCBITaxDumpMergedNodes filePath

-- | parse NCBITaxDumpNames from input string
parseNCBITaxDumpNames input = parse genParserNCBITaxDumpNames "parseTaxDumpNames" input

-- | parse NCBITaxDumpNames from input filePath                      
readNCBITaxDumpNames :: String -> IO (Either ParseError [TaxDumpName])  
readNCBITaxDumpNames filePath = parseFromFile genParserNCBITaxDumpNames filePath

-- | parse NCBITaxDumpNames from input string
parseNCBITaxDumpNodes input = parse genParserNCBITaxDumpNode "parseTaxDumpNode" input

-- | parse NCBITaxDumpCitations from input filePath                      
readNCBITaxDumpNodes :: String -> IO (Either ParseError [TaxDumpNode])  
readNCBITaxDumpNodes filePath = parseFromFile genParserNCBITaxDumpNodes filePath


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

genParserNCBITaxDumpCitations :: GenParser Char st [TaxDumpCitation]
genParserNCBITaxDumpCitations = do
  citations <- many1 genParserNCBITaxDumpCitation
  return $ citations

genParserNCBITaxDumpDelNodes :: GenParser Char st [TaxDumpDelNode]
genParserNCBITaxDumpDelNodes = do
  delNodes <- many1 genParserNCBITaxDumpDelNode
  return $ delNodes
  
genParserNCBITaxDumpDivisons :: GenParser Char st [TaxDumpDivision]
genParserNCBITaxDumpDivisons = do
  divisions <- many1 genParserNCBITaxDumpDivision
  return $ divisions

genParserNCBITaxDumpGenCodes :: GenParser Char st [TaxDumpGenCode]
genParserNCBITaxDumpGenCodes = do
  genCodes <- many1 genParserNCBITaxDumpGenCode
  return $ genCodes

genParserNCBITaxDumpMergedNodes :: GenParser Char st [TaxDumpMergedNode]
genParserNCBITaxDumpMergedNodes = do
  mergedNodes <- many1 genParserNCBITaxDumpMergedNode
  return $ mergedNodes

genParserNCBITaxDumpNames :: GenParser Char st [TaxDumpName]
genParserNCBITaxDumpNames = do
  names <- many1 genParserNCBITaxDumpName
  return $ names

genParserNCBITaxDumpNodes :: GenParser Char st [TaxDumpNode]
genParserNCBITaxDumpNodes = do
  nodes <- many1 genParserNCBITaxDumpNode
  return $ nodes

----------------------------

genParserNCBITaxDumpCitation :: GenParser Char st TaxDumpCitation
genParserNCBITaxDumpCitation = do
  citId <- many1 digit
  tab
  char ('|')
  tab 
  citKey <- optionMaybe (many1 (noneOf "\t"))
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
  url <- optionMaybe genParserTaxURL
  --optional (string "ë|")
  tab
  optional (skipMany (noneOf ("|")))
  char ('|')
  tab
  text <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab
  taxIdList <- optionMaybe (many1 genParserTaxIdList)
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpCitation (readInt citId) citKey (liftM readInt pubmedId) (liftM readInt medlineId) url text taxIdList

genParserNCBITaxDumpDelNode :: GenParser Char st TaxDumpDelNode
genParserNCBITaxDumpDelNode = do
  delNode <- many1 digit
  space
  char ('|')
  char ('\n')
  return $ TaxDumpDelNode (readInt delNode)
  
genParserNCBITaxDumpDivision :: GenParser Char st TaxDumpDivision
genParserNCBITaxDumpDivision = do
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

genParserNCBITaxDumpGenCode :: GenParser Char st TaxDumpGenCode
genParserNCBITaxDumpGenCode = do
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

genParserNCBITaxDumpMergedNode :: GenParser Char st TaxDumpMergedNode
genParserNCBITaxDumpMergedNode = do
  oldTaxId <- many1 digit
  tab
  char ('|')
  tab  
  newTaxId <- many1 digit
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpMergedNode (readInt oldTaxId) (readInt newTaxId)

genParserNCBITaxDumpName :: GenParser Char st TaxDumpName
genParserNCBITaxDumpName = do
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

genParserNCBITaxDumpNode :: GenParser Char st TaxDumpNode
genParserNCBITaxDumpNode = do
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
