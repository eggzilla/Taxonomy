-- | Parse and process taxonomy data

module Bio.Taxonomy (                      
                       module Bio.TaxonomyData,
                       genParserTaxonomyGraph,
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
                       parseNCBISimpleTaxons,
                       readNCBISimpleTaxons,
                       readNCBITaxonomyDatabaseDump,
                       constructTaxTree,
                       constructSimpleTaxTree
                      ) where
import Prelude 
import System.IO 
import Bio.TaxonomyData
import Text.Parsec.Prim (runP)
import Text.ParserCombinators.Parsec
import Control.Monad
import Data.Tree
import Data.List
import Data.Maybe    
import Data.Either
import qualified Data.Either.Unwrap as E
import Data.Graph.Inductive
--import qualified Data.GraphViz as GV
--import qualified Data.GraphViz.Printing as GVP

--------------------------------------------------------
--fgl graph representation
genParserTaxonomyGraph :: GenParser Char Int (Gr SimpleTaxon Double)
genParserTaxonomyGraph = do
  nodesEdges <- many1 (try (genParserGraphNodeEdge))
  optional eof
  let (nodesList,edgesList) =  unzip nodesEdges
  let taxedges = concat edgesList
  let taxnodes = concat nodesList
  return (mkGraph taxnodes taxedges)

genParserGraphNodeEdge :: GenParser Char st ([(Int,SimpleTaxon)],[(Int,Int,Double)])
genParserGraphNodeEdge = do
  _simpleTaxId <- many1 digit
  tab
  char ('|') 
  tab
  _simpleParentTaxId <- many1 digit
  tab
  char ('|')
  tab
  _simpleRank <- many1 (noneOf "\t")
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
  return $ ([((readInt _simpleTaxId),SimpleTaxon (readInt _simpleTaxId) [] (readInt _simpleParentTaxId) (readRank _simpleRank))],[((readInt _simpleTaxId),(readInt _simpleParentTaxId),(1 :: Double))])

----------------------------
-- Data.Tree representation
constructSimpleTaxTree :: [SimpleTaxon] -> Tree SimpleTaxon
constructSimpleTaxTree (taxnode:taxnodes) = Node taxnode (concat (addSimpleChildElements (simpleTaxId taxnode) taxnodes))
              

addSimpleChildElements :: Int -> [SimpleTaxon] -> [[Tree SimpleTaxon]]
addSimpleChildElements currentTaxId taxnodes = do
  let (childElements, remainingElements) = partition (\x -> simpleParentTaxId x == currentTaxId) taxnodes
  let subtreeLists = map (\x -> (x:remainingElements)) childElements
  let subtrees = constructSimpleSubTrees subtreeLists
  return subtrees
         
constructSimpleSubTrees :: [[SimpleTaxon]] -> [Tree SimpleTaxon]
constructSimpleSubTrees subtreeLists =  map constructSimpleTaxTree subtreeLists

constructTaxTree :: [TaxDumpNode] -> Tree TaxDumpNode
constructTaxTree (taxnode:taxnodes) = Node taxnode (concat (addChildElements (taxId taxnode) taxnodes))

addChildElements :: Int -> [TaxDumpNode] -> [[Tree TaxDumpNode]]
addChildElements currentTaxId taxnodes = do
  let (childElements, remainingElements) = partition (\x -> parentTaxId x == currentTaxId) taxnodes
  let subtreeLists = map (\x -> (x:remainingElements)) childElements
  let subtrees = constructSubTrees subtreeLists
  return subtrees
         
constructSubTrees :: [[TaxDumpNode]] -> [Tree TaxDumpNode]
constructSubTrees subtreeLists =  map constructTaxTree subtreeLists

-- | parse NCBITaxDumpCitations from input string
parseNCBITaxDumpCitations :: [Char] -> Either ParseError [TaxDumpCitation]
parseNCBITaxDumpCitations input = parse genParserNCBITaxDumpCitations "parseTaxDumpCitations" input

-- | parse NCBITaxDumpCitations from input filePath                      
readNCBITaxDumpCitations :: String -> IO (Either ParseError [TaxDumpCitation])  
readNCBITaxDumpCitations filePath = parseFromFileEncISO88591 genParserNCBITaxDumpCitations filePath

-- | parse NCBITaxDumpDelNodes from input string
parseNCBITaxDumpDelNodes :: [Char] -> Either ParseError [TaxDumpDelNode]
parseNCBITaxDumpDelNodes input = parse genParserNCBITaxDumpDelNodes "parseTaxDumpDelNodes" input

-- | parse NCBITaxDumpDelNodes from input filePath                      
readNCBITaxDumpDelNodes :: String -> IO (Either ParseError [TaxDumpDelNode])  
readNCBITaxDumpDelNodes filePath = parseFromFile genParserNCBITaxDumpDelNodes filePath

-- | parse NCBITaxDumpDivisons from input string
parseNCBITaxDumpDivisions :: [Char] -> Either ParseError [TaxDumpDivision]
parseNCBITaxDumpDivisions input = parse genParserNCBITaxDumpDivisons "parseTaxDumpDivisons" input

-- | parse NCBITaxDumpDivisons from input filePath                      
readNCBITaxDumpDivisions :: String -> IO (Either ParseError [TaxDumpDivision])  
readNCBITaxDumpDivisions filePath = parseFromFile genParserNCBITaxDumpDivisons filePath

-- | parse NCBITaxDumpGenCodes from input string
parseNCBITaxDumpGenCodes :: [Char] -> Either ParseError [TaxDumpGenCode]
parseNCBITaxDumpGenCodes input = parse genParserNCBITaxDumpGenCodes "parseTaxDumpGenCodes" input

-- | parse NCBITaxDumpGenCodes from input filePath                      
readNCBITaxDumpGenCodes :: String -> IO (Either ParseError [TaxDumpGenCode])  
readNCBITaxDumpGenCodes filePath = parseFromFile genParserNCBITaxDumpGenCodes filePath

-- | parse NCBITaxDumpMergedNodes from input string
parseNCBITaxDumpMergedNodes :: [Char] -> Either ParseError [TaxDumpMergedNode]
parseNCBITaxDumpMergedNodes input = parse genParserNCBITaxDumpMergedNodes "parseTaxDumpMergedNodes" input

-- | parse NCBITaxDumpMergedNodes from input filePath                      
readNCBITaxDumpMergedNodes :: String -> IO (Either ParseError [TaxDumpMergedNode])  
readNCBITaxDumpMergedNodes filePath = parseFromFile genParserNCBITaxDumpMergedNodes filePath

-- | parse NCBITaxDumpNames from input string
parseNCBITaxDumpNames :: [Char] -> Either ParseError [TaxDumpName]
parseNCBITaxDumpNames input = parse genParserNCBITaxDumpNames "parseTaxDumpNames" input

-- | parse NCBITaxDumpNames from input filePath                      
readNCBITaxDumpNames :: String -> IO (Either ParseError [TaxDumpName])  
readNCBITaxDumpNames filePath = parseFromFile genParserNCBITaxDumpNames filePath

-- | parse NCBITaxDumpNames from input string
parseNCBITaxDumpNodes :: [Char] -> Either ParseError TaxDumpNode
parseNCBITaxDumpNodes input = parse genParserNCBITaxDumpNode "parseTaxDumpNode" input

-- | parse NCBITaxDumpCitations from input filePath                      
readNCBITaxDumpNodes :: String -> IO (Either ParseError [TaxDumpNode])  
readNCBITaxDumpNodes filePath = parseFromFile genParserNCBITaxDumpNodes filePath

-- | parse NCBISimpleTaxDumpNames from input string
parseNCBISimpleTaxons :: [Char] -> Either ParseError SimpleTaxon
parseNCBISimpleTaxons input = parse genParserNCBISimpleTaxon "parseSimpleTaxon" input

-- | parse NCBITaxDumpCitations from input filePath                      
readNCBISimpleTaxons :: String -> IO (Either ParseError [SimpleTaxon])  
readNCBISimpleTaxons filePath = parseFromFile genParserNCBISimpleTaxons filePath

-- | Parse the input as NCBITaxDump datatype
readNCBITaxonomyDatabaseDump :: String -> IO (Either [[Char]] NCBITaxDump)
readNCBITaxonomyDatabaseDump folder = do
  citations <- readNCBITaxDumpCitations (folder ++ "citations.dmp")
  let citationsError = extractParseError citations
  taxdelNodes <- readNCBITaxDumpDelNodes (folder ++ "delnodes.dmp")
  let delNodesError = extractParseError taxdelNodes
  divisons <- readNCBITaxDumpDivisions (folder ++ "division.dmp")
  let divisonsError = extractParseError divisons
  genCodes <- readNCBITaxDumpGenCodes (folder ++ "gencode.dmp")
  let genCodesError = extractParseError genCodes
  mergedNodes <- readNCBITaxDumpMergedNodes (folder ++ "merged.dmp")
  let mergedNodesError = extractParseError mergedNodes
  names <- readNCBITaxDumpNames (folder ++ "names.dmp")
  let namesError = extractParseError names
  taxnodes <- readNCBITaxDumpNodes (folder ++ "nodes.dmp") 
  let nodesError = extractParseError taxnodes
  let parseErrors =  [citationsError, delNodesError, divisonsError, genCodesError, mergedNodesError, namesError, nodesError]
  return $ (checkParsing parseErrors citations taxdelNodes divisons genCodes mergedNodes names taxnodes)

genParserNCBITaxDumpCitations :: GenParser Char st [TaxDumpCitation]
genParserNCBITaxDumpCitations = do
  citations <- many1 genParserNCBITaxDumpCitation
  return $ citations

genParserNCBITaxDumpDelNodes :: GenParser Char st [TaxDumpDelNode]
genParserNCBITaxDumpDelNodes = do
  taxdelNodes <- many1 genParserNCBITaxDumpDelNode
  return $ taxdelNodes
  
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
  taxnodes <- many1 genParserNCBITaxDumpNode
  return $ taxnodes

genParserNCBISimpleTaxons :: GenParser Char st [SimpleTaxon]
genParserNCBISimpleTaxons = do
  taxnodes <- many1 genParserNCBISimpleTaxon
  return $ taxnodes
----------------------------

genParserNCBITaxDumpCitation :: GenParser Char st TaxDumpCitation
genParserNCBITaxDumpCitation = do
  _citId <- many1 digit
  tab
  char ('|')
  tab 
  _citKey <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab 
  _pubmedId <- optionMaybe (many1 digit)
  tab
  char ('|')
  tab  
  _medlineId <- optionMaybe (many1 digit)
  tab
  char ('|') 
  _url <- genParserTaxURL
  char ('|')
  tab
  _text <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab
  _taxIdList <- optionMaybe (many1 genParserTaxIdList)
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpCitation (readInt _citId) _citKey (liftM readInt _pubmedId) (liftM readInt _medlineId) _url _text _taxIdList

genParserNCBITaxDumpDelNode :: GenParser Char st TaxDumpDelNode
genParserNCBITaxDumpDelNode = do
  taxdelNode <- many1 digit
  space
  char ('|')
  char ('\n')
  return $ TaxDumpDelNode (readInt taxdelNode)
  
genParserNCBITaxDumpDivision :: GenParser Char st TaxDumpDivision
genParserNCBITaxDumpDivision = do
  _divisionId <- many1 digit
  tab
  char ('|')
  tab
  _divisionCDE <- many1 upper
  tab
  char ('|')
  tab
  _divisionName <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  _comments <- optionMaybe (many1 (noneOf ("\t")))
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpDivision (readInt _divisionId) _divisionCDE _divisionName _comments 

genParserNCBITaxDumpGenCode :: GenParser Char st TaxDumpGenCode
genParserNCBITaxDumpGenCode = do
  _geneticCodeId <- many1 digit 
  tab
  char ('|')
  tab
  _abbreviation <- optionMaybe (many1 (noneOf ("\t")))
  tab
  char ('|')
  tab
  _genCodeName <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  _cde <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  _starts <- many1 (noneOf ("\t"))
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpGenCode (readInt _geneticCodeId) _abbreviation _genCodeName _cde _starts

genParserNCBITaxDumpMergedNode :: GenParser Char st TaxDumpMergedNode
genParserNCBITaxDumpMergedNode = do
  _oldTaxId <- many1 digit
  tab
  char ('|')
  tab  
  _newTaxId <- many1 digit
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpMergedNode (readInt _oldTaxId) (readInt _newTaxId)

genParserNCBITaxDumpName :: GenParser Char st TaxDumpName
genParserNCBITaxDumpName = do
  _taxId <- many1 digit
  tab
  char ('|')
  tab
  _nameTxt <- many1 (noneOf ("\t"))
  tab
  char ('|')
  tab
  _uniqueName <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab
  _nameClass <- many1 (noneOf ("\t"))
  return $ TaxDumpName (readInt _taxId) _nameTxt _uniqueName _nameClass

genParserNCBISimpleTaxon :: GenParser Char st SimpleTaxon
genParserNCBISimpleTaxon = do
  _simpleTaxId <- many1 digit
  tab
  char ('|') 
  tab
  _simpleParentTaxId <- many1 digit
  tab
  char ('|')
  tab
  _simpleRank <- many1 (noneOf "\t")
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
  return $ SimpleTaxon (readInt _simpleTaxId) [] (readInt _simpleParentTaxId) (readRank _simpleRank) 

genParserNCBITaxDumpNode :: GenParser Char st TaxDumpNode
genParserNCBITaxDumpNode = do
  _taxId <- many1 digit
  tab
  char ('|') 
  tab
  _parentTaxId <- many1 digit
  tab
  char ('|')
  tab
  _rank <- many1 (noneOf "\t")
  tab
  char ('|')
  tab 
  _emblCode <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  tab
  _divisionId <- many1 digit
  tab
  char ('|')
  tab
  _inheritedDivFlag <- many1 digit
  tab
  char ('|')
  tab 
  _geneticCodeId <- many1 digit
  tab
  char ('|')
  tab
  _inheritedGCFlag <- many1 digit
  tab
  char ('|')
  tab
  _mitochondrialGeneticCodeId <- many1 digit
  tab
  char ('|')
  tab
  _inheritedMGCFlag <- many1 digit
  tab
  char ('|')
  tab
  _genBankHiddenFlag <- many1 digit
  tab
  char ('|')
  tab
  _hiddenSubtreeRootFlag <- many1 digit 
  tab
  char ('|')
  tab
  _comments <- optionMaybe (many1 (noneOf "\t"))
  tab
  char ('|')
  char ('\n')
  return $ TaxDumpNode (readInt _taxId) (readInt _parentTaxId) (readRank _rank) _emblCode _divisionId (readBool _inheritedDivFlag) _geneticCodeId (readBool _inheritedGCFlag) _mitochondrialGeneticCodeId (readBool _inheritedMGCFlag) (readBool _genBankHiddenFlag) (readBool _hiddenSubtreeRootFlag) _comments

--Auxiliary functions
readInt :: String -> Int
readInt = read

readBool :: String -> Bool
readBool "0" = False
readBool "1" = True
readBool _ = False               

readRank :: String -> Rank
readRank a = read  a :: Rank

genParserTaxIdList :: GenParser Char st Int
genParserTaxIdList = do
  optional (char ' ')
  _taxId <- many1 digit
  optional (char ' ')
  return $ (readInt _taxId)

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

readEncodedFile :: TextEncoding -> FilePath -> IO String                    
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
checkParsing parseErrors citations taxdelNodes divisons genCodes mergedNodes names taxnodes
  | join (parseErrors) == "" = Right (NCBITaxDump (E.fromRight citations) (E.fromRight taxdelNodes) (E.fromRight divisons) (E.fromRight genCodes) (E.fromRight mergedNodes) (E.fromRight names) (E.fromRight taxnodes))
  | otherwise = Left (parseErrors)

extractParseError :: Either ParseError a -> String
extractParseError _parse
  | isLeft _parse = show (E.fromLeft _parse)
  | otherwise = ""
