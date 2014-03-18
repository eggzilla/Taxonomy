-- | Parse and process taxonomy data

module Bio.Taxonomy (                      
                       module Bio.TaxonomyData,
                       parseNCBITaxDumpCitations,
                       readNCBITaxDumpCitations,
                       parseNCBITaxDumpDelNodes,
                       readNCBITaxDumpDelNodes,
                       parseNCBITaxDumpDivisions,
                       readNCBITaxDumpDivisions,
                       parseNCBITaxDumpGenCodes,
                       readNCBITaxDumpGenCodes,
                       parseNCBITaxDumpMergedNodes,
                       readNCBITaxDumpMergedNodes,
                       parseNCBITaxDumpNames,
                       readNCBITaxDumpNames,
                       parseNCBITaxDumpNodes,
                       readNCBITaxDumpNodes,
                       readNCBITaxonomyDatabaseDump,
                       constructTaxTree
                      ) where
import Prelude 
import System.IO 
import Text.Parsec.Prim
import Bio.TaxonomyData
import Data.Maybe
import Text.ParserCombinators.Parsec
import Text.ParserCombinators.Parsec.Token
import Text.ParserCombinators.Parsec.Language (emptyDef)    
import Control.Monad
import Data.Tree
import Data.List
import Data.Either
import Data.Either.Unwrap 

--------------------------------------------------------
--data TaxTree = TaxLeaf TaxDumpNode | TaxNode TaxDumpNode [TaxTree] deriving (Eq,Read,Show)

constructTaxTree :: [TaxDumpNode] -> Tree TaxDumpNode
constructTaxTree (node:nodes) = Node node (concat (addChildElements (taxId node) nodes))

addChildElements :: Int -> [TaxDumpNode] -> [[Tree TaxDumpNode]]
addChildElements currentTaxId nodes = do
  let (childElements, remainingElements) = partition (\x -> parentTaxId x == currentTaxId) nodes
  let subtreeLists = map (\x -> (x:remainingElements)) childElements
  let subtrees = constructSubTrees subtreeLists
  return subtrees
         
constructSubTrees :: [[TaxDumpNode]] -> [Tree TaxDumpNode]
constructSubTrees subtreeLists =  map constructTaxTree subtreeLists

-- | parse NCBITaxDumpCitations from input string
parseNCBITaxDumpCitations input = parse genParserNCBITaxDumpCitations "parseTaxDumpCitations" input

-- | parse NCBITaxDumpCitations from input filePath                      
readNCBITaxDumpCitations :: String -> IO (Either ParseError [TaxDumpCitation])  
readNCBITaxDumpCitations filePath = parseFromFileEncISO88591 genParserNCBITaxDumpCitations filePath

-- | parse NCBITaxDumpDelNodes from input string
parseNCBITaxDumpDelNodes input = parse genParserNCBITaxDumpDelNodes "parseTaxDumpDelNodes" input

-- | parse NCBITaxDumpDelNodes from input filePath                      
readNCBITaxDumpDelNodes :: String -> IO (Either ParseError [TaxDumpDelNode])  
readNCBITaxDumpDelNodes filePath = parseFromFile genParserNCBITaxDumpDelNodes filePath

-- | parse NCBITaxDumpDivisons from input string
parseNCBITaxDumpDivisions input = parse genParserNCBITaxDumpDivisons "parseTaxDumpDivisons" input

-- | parse NCBITaxDumpDivisons from input filePath                      
readNCBITaxDumpDivisions :: String -> IO (Either ParseError [TaxDumpDivision])  
readNCBITaxDumpDivisions filePath = parseFromFile genParserNCBITaxDumpDivisons filePath

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

-- | parse NCBISimpleTaxDumpNames from input string
parseNCBISimpleTaxDumpNodes input = parse genParserNCBISimpleTaxDumpNode "parseSimpleTaxDumpNode" input

-- | parse NCBITaxDumpCitations from input filePath                      
readNCBISimpleTaxDumpNodes :: String -> IO (Either ParseError [SimpleTaxDumpNode])  
readNCBISimpleTaxDumpNodes filePath = parseFromFile genParserNCBISimpleTaxDumpNodes filePath

-- | Parse the input as NCBITaxDump datatype
readNCBITaxonomyDatabaseDump :: String -> IO (Either [[Char]] NCBITaxDump)
readNCBITaxonomyDatabaseDump folder = do
  citations <- readNCBITaxDumpCitations (folder ++ "citations.dmp")
  let citationsError = extractParseError citations
  delNodes <- readNCBITaxDumpDelNodes (folder ++ "delnodes.dmp")
  let delNodesError = extractParseError delNodes
  divisons <- readNCBITaxDumpDivisions (folder ++ "division.dmp")
  let divisonsError = extractParseError divisons
  genCodes <- readNCBITaxDumpGenCodes (folder ++ "gencode.dmp")
  let genCodesError = extractParseError genCodes
  mergedNodes <- readNCBITaxDumpMergedNodes (folder ++ "merged.dmp")
  let mergedNodesError = extractParseError mergedNodes
  names <- readNCBITaxDumpNames (folder ++ "names.dmp")
  let namesError = extractParseError names
  nodes <- readNCBITaxDumpNodes (folder ++ "nodes.dmp") 
  let nodesError = extractParseError nodes
  let parseErrors =  [citationsError, delNodesError, divisonsError, genCodesError, mergedNodesError, namesError, nodesError]
  return $ (checkParsing parseErrors citations delNodes divisons genCodes mergedNodes names nodes)

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

genParserNCBISimpleTaxDumpNodes :: GenParser Char st [SimpleTaxDumpNode]
genParserNCBISimpleTaxDumpNodes = do
  nodes <- many1 genParserNCBISimpleTaxDumpNode
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
  url <- genParserTaxURL
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
  tab
  abbreviation <- optionMaybe (many1 (noneOf ("\t")))
  tab
  char ('|')
  tab
  genCodeName <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  cde <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  starts <- many1 (noneOf ("\t"))
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpGenCode (readInt geneticCodeId) abbreviation genCodeName cde starts

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

genParserNCBISimpleTaxDumpNode :: GenParser Char st SimpleTaxDumpNode
genParserNCBISimpleTaxDumpNode = do
  simpleTaxId <- many1 digit
  tab
  char ('|') 
  tab
  simpleParentTaxId <- many1 digit
  tab
  char ('|')
  tab
  simpleRank <- many1 (noneOf "\t")
  tab
  char ('|')
  tab 
  optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab
  many1 digit
  tab
  char ('|')
  tab
  many1 digit
  tab
  char ('|')
  tab 
  many1 digit
  tab
  char ('|')
  tab
  many1 digit
  tab
  char ('|')
  tab
  many1 digit
  tab
  char ('|')
  tab
  many1 digit
  tab
  char ('|')
  tab
  many1 digit
  tab
  char ('|')
  tab
  many1 digit 
  tab
  char ('|')
  tab
  optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  char ('\n')
  return $ SimpleTaxDumpNode (readInt simpleTaxId) (readInt simpleParentTaxId) (readRank simpleRank) 

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
  comments <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpNode (readInt taxId) (readInt parentTaxId) (readRank rank) emblCode divisionId (readBool inheritedDivFlag) geneticCodeId (readBool inheritedGCFlag) mitochondrialGeneticCodeId (readBool inheritedMGCFlag) (readBool genBankHiddenFlag) (readBool hiddenSubtreeRootFlag) comments

--Auxiliary functions
readDouble :: String -> Double
readDouble = read              

readInt :: String -> Int
readInt = read

readBool :: String -> Bool
readBool "0" = False
readBool "1" = True

readRank :: String -> Rank
readRank a = read  a :: Rank

genParserTaxIdList :: GenParser Char st Int
genParserTaxIdList = do
  optional (char ' ')
  taxId <- many1 digit
  optional (char ' ')
  return $ (readInt taxId)

genParserTaxURL :: GenParser Char st (Maybe String)
genParserTaxURL = do
  tab 
  url1 <- optionMaybe (many1 (noneOf "\t"))
  tab
  url2 <- optionMaybe (many1 (noneOf ("|")))
  return $ (concatenateURLParts url1 url2)

concatenateURLParts :: Maybe String -> Maybe String -> Maybe String
concatenateURLParts url1 url2 
  | (isJust url1) && (isJust url2) = maybeStringConcat url1 url2
  | (isJust url1) && (isNothing url2) = url1
  | otherwise = Nothing 

maybeStringConcat :: Maybe String -> Maybe String -> Maybe String
maybeStringConcat = liftM2 (++)

readEncodedFile encoding name = do 
  handle <- openFile name ReadMode
  hSetEncoding handle encoding
  hGetContents handle

parseFromFileEncISO88591 :: Parser a -> String -> IO (Either ParseError a)
parseFromFileEncISO88591 parser fname = do 
         input <- readEncodedFile latin1 fname
         return (runP parser () fname input)

-- | check a list of parsing results for presence of Left aka Parse error
checkParsing :: [[Char]] -> Either ParseError [TaxDumpCitation] -> Either ParseError [TaxDumpDelNode] -> Either ParseError [TaxDumpDivision] -> Either ParseError [TaxDumpGenCode] -> Either ParseError [TaxDumpMergedNode] -> Either ParseError [TaxDumpName] -> Either ParseError [TaxDumpNode]-> Either [[Char]] NCBITaxDump
checkParsing parseErrors citations delNodes divisons genCodes mergedNodes names nodes
  | join (parseErrors) == "" = Right (NCBITaxDump (fromRight citations) (fromRight delNodes) (fromRight divisons) (fromRight genCodes) (fromRight mergedNodes) (fromRight names) (fromRight nodes))
  | otherwise = Left (parseErrors)

extractParseError :: Either ParseError a -> String
extractParseError parse
  | isLeft parse = show (fromLeft parse)
  | otherwise = ""
