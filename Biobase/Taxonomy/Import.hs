-- | Functions for parsing, processing and visualization of taxonomy data.
--
-- === Usage example:
-- * Read in taxonomy data
--
--     > eitherTaxtree <- readNamedTaxonomy "/path/to/NCBI_taxonomydump_directory"
--
-- * Process data
--
--     > let subtree = extractTaxonomySubTreebyLevel [562] (fromRight eitherTaxTree) (Just 4)
--
-- * Visualize result
--
--     tput "/path/to/dotdirectory" subtree
module Biobase.Taxonomy.Import (  -- * Datatypes
                       -- Datatypes used to represent taxonomy data
                       module Biobase.Taxonomy.Types,
                       -- * Parsing
                       -- Functions prefixed with "read" read from filepaths, functions with parse from Haskell Strings.
                       readTaxonomy,
                       readNamedTaxonomy,
                       parseTaxonomy,
                       parseNCBITaxCitations,
                       readNCBITaxCitations,
                       parseNCBITaxDelNodes,
                       readNCBITaxDelNodes,
                       parseNCBITaxDivisions,
                       readNCBITaxDivisions,
                       parseNCBITaxGenCodes,
                       readNCBITaxGenCodes,
                       parseNCBITaxMergedNodes,
                       readNCBITaxMergedNodes,
                       parseNCBITaxNames,
                       readNCBITaxNames,
                       parseNCBITaxNodes,
                       readNCBITaxNodes,
                       parseNCBISimpleTaxons,
                       readNCBISimpleTaxons,
                       readNCBITaxonomyDatabase
                      ) where
import Prelude
import System.IO
import Biobase.Taxonomy.Types
import Text.Parsec.Prim (runP)
import Text.ParserCombinators.Parsec
import Control.Monad
import Data.List
import Data.Maybe
import qualified Data.Either.Unwrap as E
import Data.Graph.Inductive.Graph
import Data.Graph.Inductive.Tree
import qualified Data.ByteString.Char8 as B
import qualified Data.Text.Lazy as T
--------------------------------------------------------

---------------------------------------
-- Parsing functions

-- | NCBI taxonomy dump nodes and names in the input directory path are parsed and a SimpleTaxon tree with scientific names for each node is generated.
readNamedTaxonomy :: String -> IO (Either ParseError (Gr SimpleTaxon Double))
readNamedTaxonomy directoryPath = do
  nodeNames <- readNCBITaxNames (directoryPath ++ "names.dmp")
  if E.isLeft nodeNames
     then return (Left (E.fromLeft nodeNames))
     else do
       let rightNodeNames = E.fromRight nodeNames
       let filteredNodeNames = filter isScientificName rightNodeNames
       let namedTaxonomyGraph = genParserNamedTaxonomyGraph filteredNodeNames
       parseFromFileEncISO88591 namedTaxonomyGraph (directoryPath ++ "nodes.dmp")

isScientificName :: TaxName -> Bool
isScientificName name = nameClass name == scientificNameT
  where scientificNameT = B.pack "scientific name"

-- | NCBI taxonomy dump nodes and names in the input directory path are parsed and a SimpleTaxon tree is generated.
readTaxonomy :: String -> IO (Either ParseError (Gr SimpleTaxon Double))
readTaxonomy = parseFromFileEncISO88591 genParserTaxonomyGraph

-- | NCBI taxonomy dump nodes and names in the input directory path are parsed and a SimpleTaxon tree is generated.
parseTaxonomy :: String -> Either ParseError (Gr SimpleTaxon Double)
parseTaxonomy = parse genParserTaxonomyGraph "parseTaxonomy"

genParserTaxonomyGraph :: GenParser Char st (Gr SimpleTaxon Double)
genParserTaxonomyGraph = do
  nodesEdges <- many1 (try genParserGraphNodeEdge)
  optional eof
  let (nodesList,edgesList) =  unzip nodesEdges
  --let taxedges = filter (\(a,b,_) -> a /= b) edgesList
  let taxedges = filter notLoopEdge  edgesList
  --let taxnodes = concat nodesList
  --return (mkGraph taxnodes taxedges)
  let currentGraph = mkGraph nodesList taxedges
  return currentGraph


notLoopEdge :: (Int,Int,a) -> Bool
notLoopEdge (a,b,_) = a /= b

--genParserNodeEdges :: [TaxName] -> GenParser Char st [(Int,SimpleTaxon),(Int,Int,Double)]
--genParserNodeEdges = do
--  nodesEdges <- (many1 (try genParserGraphNodeEdge))
--  optional eof
--  return (nodesList,edgesList)


  --let taxedges = filter notLoopEdge edgesList
  --let taxnamednodes = map (setNodeScientificName filteredNodeNames) nodesList
  --let currentGraph = mkGraph taxnamednodes taxedges
  --return currentGraph

genParserNamedTaxonomyGraph :: [TaxName] -> GenParser Char st (Gr SimpleTaxon Double)
genParserNamedTaxonomyGraph filteredNodeNames = do
  nodesEdges <- (many1 (try genParserGraphNodeEdge))
  optional eof
  let (nodesList,edgesList) = unzip nodesEdges
  let taxedges = filter notLoopEdge edgesList
  let taxnamednodes = map (setNodeScientificName filteredNodeNames) nodesList
  let currentGraph = mkGraph taxnamednodes taxedges
  return currentGraph

setNodeScientificName :: [TaxName] -> (t, SimpleTaxon) -> (t, SimpleTaxon)
setNodeScientificName inputTaxNames (inputNode,inputTaxon) = outputNode
  where maybeRetrievedName = find (isTaxNameIdSimpleTaxid inputTaxon) inputTaxNames
        retrievedName = maybe (T.pack "no name") nameTxt maybeRetrievedName
        outputNode = (inputNode,inputTaxon{simpleScientificName = retrievedName})

isTaxNameIdSimpleTaxid :: SimpleTaxon -> TaxName -> Bool
isTaxNameIdSimpleTaxid inputTaxon inputTaxName = nameTaxId inputTaxName == simpleTaxId inputTaxon


genParserGraphNodeEdge :: GenParser Char st ((Int,SimpleTaxon),(Int,Int,Double))
genParserGraphNodeEdge = do
  _simpleTaxId <- many1 digit
  string "\t|\t"
  _simpleParentTaxId <- many1 digit
  string "\t|\t"
  _simpleRank <- many1 (noneOf "\t")
  many1 (noneOf "\n")
  char '\n'
  let _simpleTaxIdInt = readInt _simpleTaxId
  let _simpleParentTaxIdInt = readInt _simpleParentTaxId
  return ((_simpleTaxIdInt,SimpleTaxon _simpleTaxIdInt T.empty _simpleParentTaxIdInt (readRank _simpleRank)),(_simpleTaxIdInt,_simpleParentTaxIdInt,1 :: Double))

-- | parse NCBITaxCitations from input string
parseNCBITaxCitations :: String -> Either ParseError [TaxCitation]
parseNCBITaxCitations = parse genParserNCBITaxCitations "parseTaxCitations"

-- | parse NCBITaxCitations from input filePath
readNCBITaxCitations :: String -> IO (Either ParseError [TaxCitation])
readNCBITaxCitations = parseFromFileEncISO88591 genParserNCBITaxCitations

-- | parse NCBITaxDelNodes from input string
parseNCBITaxDelNodes :: String -> Either ParseError [TaxDelNode]
parseNCBITaxDelNodes = parse genParserNCBITaxDelNodes "parseTaxDelNodes"

-- | parse NCBITaxDelNodes from input filePath
readNCBITaxDelNodes :: String -> IO (Either ParseError [TaxDelNode])
readNCBITaxDelNodes = parseFromFile genParserNCBITaxDelNodes

-- | parse NCBITaxDivisons from input string
parseNCBITaxDivisions :: String -> Either ParseError [TaxDivision]
parseNCBITaxDivisions = parse genParserNCBITaxDivisons "parseTaxDivisons"

-- | parse NCBITaxDivisons from input filePath
readNCBITaxDivisions :: String -> IO (Either ParseError [TaxDivision])
readNCBITaxDivisions = parseFromFile genParserNCBITaxDivisons

-- | parse NCBITaxGenCodes from input string
parseNCBITaxGenCodes :: String -> Either ParseError [TaxGenCode]
parseNCBITaxGenCodes = parse genParserNCBITaxGenCodes "parseTaxGenCodes"

-- | parse NCBITaxGenCodes from input filePath
readNCBITaxGenCodes :: String -> IO (Either ParseError [TaxGenCode])
readNCBITaxGenCodes = parseFromFile genParserNCBITaxGenCodes

-- | parse NCBITaxMergedNodes from input string
parseNCBITaxMergedNodes :: String -> Either ParseError [TaxMergedNode]
parseNCBITaxMergedNodes = parse genParserNCBITaxMergedNodes "parseTaxMergedNodes"

-- | parse NCBITaxMergedNodes from input filePath
readNCBITaxMergedNodes :: String -> IO (Either ParseError [TaxMergedNode])
readNCBITaxMergedNodes = parseFromFile genParserNCBITaxMergedNodes

-- | parse NCBITaxNames from input string
parseNCBITaxNames :: String -> Either ParseError [TaxName]
parseNCBITaxNames = parse genParserNCBITaxNames "parseTaxNames"

-- | parse NCBITaxNames from input filePath
readNCBITaxNames :: String -> IO (Either ParseError [TaxName])
readNCBITaxNames = parseFromFile genParserNCBITaxNames

-- | parse NCBITaxNames from input string
parseNCBITaxNodes :: String -> Either ParseError TaxNode
parseNCBITaxNodes = parse genParserNCBITaxNode "parseTaxNode"

-- | parse NCBITaxCitations from input filePath
readNCBITaxNodes :: String -> IO (Either ParseError [TaxNode])
readNCBITaxNodes = parseFromFile genParserNCBITaxNodes

-- | parse NCBISimpleTaxNames from input string
parseNCBISimpleTaxons :: String -> Either ParseError SimpleTaxon
parseNCBISimpleTaxons = parse genParserNCBISimpleTaxon "parseSimpleTaxon"

-- | parse NCBITaxCitations from input filePath
readNCBISimpleTaxons :: String -> IO (Either ParseError [SimpleTaxon])
readNCBISimpleTaxons = parseFromFile genParserNCBISimpleTaxons

-- | Parse the input as NCBITax datatype
readNCBITaxonomyDatabase :: String -> IO (Either [String] NCBITaxDump)
readNCBITaxonomyDatabase folder = do
  citations <- readNCBITaxCitations (folder ++ "citations.dmp")
  let citationsError = extractParseError citations
  taxdelNodes <- readNCBITaxDelNodes (folder ++ "delnodes.dmp")
  let delNodesError = extractParseError taxdelNodes
  divisons <- readNCBITaxDivisions (folder ++ "division.dmp")
  let divisonsError = extractParseError divisons
  genCodes <- readNCBITaxGenCodes (folder ++ "gencode.dmp")
  let genCodesError = extractParseError genCodes
  mergedNodes <- readNCBITaxMergedNodes (folder ++ "merged.dmp")
  let mergedNodesError = extractParseError mergedNodes
  names <- readNCBITaxNames (folder ++ "names.dmp")
  let namesError = extractParseError names
  taxnodes <- readNCBITaxNodes (folder ++ "nodes.dmp")
  let nodesError = extractParseError taxnodes
  let parseErrors =  [citationsError, delNodesError, divisonsError, genCodesError, mergedNodesError, namesError, nodesError]
  return (checkParsing parseErrors citations taxdelNodes divisons genCodes mergedNodes names taxnodes)

genParserNCBITaxCitations :: GenParser Char st [TaxCitation]
genParserNCBITaxCitations = many1 genParserNCBITaxCitation

genParserNCBITaxDelNodes :: GenParser Char st [TaxDelNode]
genParserNCBITaxDelNodes = many1 genParserNCBITaxDelNode

genParserNCBITaxDivisons :: GenParser Char st [TaxDivision]
genParserNCBITaxDivisons = many1 genParserNCBITaxDivision

genParserNCBITaxGenCodes :: GenParser Char st [TaxGenCode]
genParserNCBITaxGenCodes = many1 genParserNCBITaxGenCode


genParserNCBITaxMergedNodes :: GenParser Char st [TaxMergedNode]
genParserNCBITaxMergedNodes = many1 genParserNCBITaxMergedNode


genParserNCBITaxNames :: GenParser Char st [TaxName]
genParserNCBITaxNames = many1 genParserNCBITaxName

genParserNCBITaxNodes :: GenParser Char st [TaxNode]
genParserNCBITaxNodes = many1 genParserNCBITaxNode

genParserNCBISimpleTaxons :: GenParser Char st [SimpleTaxon]
genParserNCBISimpleTaxons = many1 genParserNCBISimpleTaxon


genParserNCBITaxCitation :: GenParser Char st TaxCitation
genParserNCBITaxCitation = do
  _citId <- many1 digit
  string "\t|\t"
  _citKey <- many (noneOf "\t")
  string "\t|\t"
  _pubmedId <- optionMaybe (many1 digit)
  string "\t|\t"
  _medlineId <- optionMaybe (many1 digit)
  tab
  char '|'
  _url <- genParserTaxURL
  char '|'
  tab
  _text <- (many (noneOf "\t"))
  string "\t|\t"
  _taxIdList <- (many genParserTaxIdList)
  string "\t|\n"
  return $ TaxCitation (readInt _citId) (B.pack _citKey) (liftM readInt _pubmedId) (liftM readInt _medlineId) _url (B.pack _text) _taxIdList

genParserNCBITaxDelNode :: GenParser Char st TaxDelNode
genParserNCBITaxDelNode = do
  taxdelNode <- many1 digit
  space
  char '|'
  char '\n'
  return $ TaxDelNode (readInt taxdelNode)

genParserNCBITaxDivision :: GenParser Char st TaxDivision
genParserNCBITaxDivision = do
  _divisionId <- many1 digit
  string "\t|\t"
  _divisionCDE <- many1 upper
  string "\t|\t"
  _divisionName <- many1 (noneOf "\t")
  string "\t|\t"
  _comments <- many1 (noneOf "\t")
  string "\t|\n"
  return $ TaxDivision (readInt _divisionId) (B.pack _divisionCDE) (B.pack _divisionName) (B.pack _comments)

genParserNCBITaxGenCode :: GenParser Char st TaxGenCode
genParserNCBITaxGenCode = do
  _geneticCodeId <- many1 digit
  string "\t|\t"
  _abbreviation <- (many1 (noneOf "\t"))
  string "\t|\t"
  _genCodeName <- many1 (noneOf "\t")
  string "\t|\t"
  _cde <- many1 (noneOf "\t")
  string "\t|\t"
  _starts <- many1 (noneOf "\t")
  string "\t|\n"
  return $ TaxGenCode (readInt _geneticCodeId) (B.pack _abbreviation) (B.pack _genCodeName) (B.pack _cde) (B.pack _starts)

genParserNCBITaxMergedNode :: GenParser Char st TaxMergedNode
genParserNCBITaxMergedNode = do
  _oldTaxId <- many1 digit
  string "\t|\t"
  _newTaxId <- many1 digit
  string "\t|\n"
  return $ TaxMergedNode (readInt _oldTaxId) (readInt _newTaxId)

genParserNCBITaxName :: GenParser Char st TaxName
genParserNCBITaxName = do
  _taxId <- many1 digit
  string "\t|\t"
  _nameTxt <- many1 (noneOf "\t\n")
  string "\t|\t"
  _uniqueName <- many (noneOf "\t\n")
  string "\t|\t"
  _nameClass <- many1 (noneOf "\t\n")
  tab
  char '|'
  newline
  return $! TaxName (readInt _taxId) (T.pack _nameTxt) (B.pack _uniqueName) (B.pack _nameClass)

genParserNCBISimpleTaxon :: GenParser Char st SimpleTaxon
genParserNCBISimpleTaxon = do
  _simpleTaxId <- many1 digit
  string "\t|\t"
  _simpleParentTaxId <- many1 digit
  string "\t|\t"
  _simpleRank <- many1 (noneOf "\t")
  many1 (noneOf "\n")
  char '\n'
  return $! SimpleTaxon (readInt _simpleTaxId) T.empty (readInt _simpleParentTaxId) (readRank _simpleRank)

genParserNCBITaxNode :: GenParser Char st TaxNode
genParserNCBITaxNode = do
  _taxId <- many1 digit
  string "\t|\t"
  _parentTaxId <- many1 digit
  string "\t|\t"
  _rank <- many1 (noneOf "\t")
  string "\t|\t"
  _emblCode <- (many (noneOf "\t"))
  string "\t|\t"
  _divisionId <- many1 digit
  string "\t|\t"
  _inheritedDivFlag <- many1 digit
  string "\t|\t"
  _geneticCodeId <- many1 digit
  string "\t|\t"
  _inheritedGCFlag <- many1 digit
  string "\t|\t"
  _mitochondrialGeneticCodeId <- many1 digit
  string "\t|\t"
  _inheritedMGCFlag <- many1 digit
  string "\t|\t"
  _genBankHiddenFlag <- many1 digit
  string "\t|\t"
  _hiddenSubtreeRootFlag <- many1 digit
  string "\t|\t"
  _comments <- many (noneOf "\t")
  tab
  char '|'
  char '\n'
  return $ TaxNode (readInt _taxId) (readInt _parentTaxId) (readRank _rank) (B.pack _emblCode) (read _divisionId :: Int) (readBool _inheritedDivFlag) (read _geneticCodeId ::Int) (readBool _inheritedGCFlag) (read _mitochondrialGeneticCodeId ::Int) (readBool _inheritedMGCFlag) (readBool _genBankHiddenFlag) (readBool _hiddenSubtreeRootFlag) (B.pack _comments)

---------------------------------------
-- Auxiliary functions
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
  return (readInt _taxId)

genParserTaxURL :: GenParser Char st B.ByteString
genParserTaxURL = do
  tab
  url1 <- many (noneOf "\t")
  tab
  url2 <- many (noneOf "|")
  return (B.pack (url1 ++ url2))
  --return (concatenateURLParts url1 url2)

concatenateURLParts :: Maybe String -> Maybe String -> Maybe String
concatenateURLParts url1 url2
  | isJust url1 && isJust url2 = maybeStringConcat url1 url2
  | isJust url1 && isNothing url2 = url1
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
checkParsing :: [String] -> Either ParseError [TaxCitation] -> Either ParseError [TaxDelNode] -> Either ParseError [TaxDivision] -> Either ParseError [TaxGenCode] -> Either ParseError [TaxMergedNode] -> Either ParseError [TaxName] -> Either ParseError [TaxNode]-> Either [String] NCBITaxDump
checkParsing parseErrors citations taxdelNodes divisons genCodes mergedNodes names taxnodes
  | join parseErrors == "" = Right (NCBITaxDump (E.fromRight citations) (E.fromRight taxdelNodes) (E.fromRight divisons) (E.fromRight genCodes) (E.fromRight mergedNodes) (E.fromRight names) (E.fromRight taxnodes))
  | otherwise = Left parseErrors

extractParseError :: Either ParseError a -> String
extractParseError _parse
  | E.isLeft _parse = show (E.fromLeft _parse)
  | otherwise = ""
