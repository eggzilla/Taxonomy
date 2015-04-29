-- | Parse and process taxonomy data

module Bio.Taxonomy (                      
                       module Bio.TaxonomyData,
                       extractTaxonomySubTree,
                       drawTaxonomy,    
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
                       readNCBITaxonomyDatabase,
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
import qualified Data.GraphViz as GV
import qualified Data.GraphViz.Printing as GVP
import qualified Data.Text.Lazy as TL

--------------------------------------------------------

--fgl graph representation

-- | draw Graph in dot format
drawTaxonomy :: Gr SimpleTaxon Double -> String
drawTaxonomy inputGraph = do
  let params = GV.nonClusteredParams {GV.isDirected       = True
                       , GV.globalAttributes = []
                       , GV.isDotCluster     = const True
                       , GV.fmtNode = \ (_,l) -> [GV.textLabel (TL.pack ((show (simpleRank l)) ++ "\n" ++ simpleScientificName l))]
                       , GV.fmtEdge          = const []
                       }
  let dotFormat = GV.graphToDot params inputGraph
  let dottext = GVP.renderDot $ GVP.toDot dotFormat
  TL.unpack dottext

-- | parse Taxonomy from input filePath                      
readNamedTaxonomy :: String -> IO (Either ParseError (Gr SimpleTaxon Double))  
readNamedTaxonomy directoryPath = do
  nodeNames <- readNCBITaxNames (directoryPath ++ "names.dmp")
  if (isLeft nodeNames)
     then do 
       return (Left (E.fromLeft nodeNames))
     else do
       let filteredNodeNames = filter (\a -> nameClass a == "scientific name") (E.fromRight nodeNames)
       taxonomyGraph <- parseFromFileEncISO88591 (genParserNamedTaxonomyGraph filteredNodeNames) (directoryPath ++ "nodes.dmp")
       return taxonomyGraph

-- | parse Taxonomy from file path
readTaxonomy :: String -> IO (Either ParseError (Gr SimpleTaxon Double))  
readTaxonomy filepath = parseFromFileEncISO88591 genParserTaxonomyGraph filepath

-- | parse Taxonomy from input string
parseTaxonomy :: [Char] -> Either ParseError (Gr SimpleTaxon Double)
parseTaxonomy input = parse genParserTaxonomyGraph "parseTaxonomy" input

genParserTaxonomyGraph :: GenParser Char st (Gr SimpleTaxon Double)
genParserTaxonomyGraph = do
  nodesEdges <- many1 (try (genParserGraphNodeEdge))
  optional eof
  let (nodesList,edgesList) =  unzip nodesEdges
  let taxedges = concat edgesList
  let taxnodes = concat nodesList
  return (mkGraph taxnodes taxedges)

genParserNamedTaxonomyGraph :: [TaxName] -> GenParser Char st (Gr SimpleTaxon Double)
genParserNamedTaxonomyGraph filteredNodeNames = do
  nodesEdges <- many1 (try (genParserGraphNodeEdge))
  optional eof
  let (nodesList,edgesList) =  unzip nodesEdges
  let taxedges = concat edgesList
  let taxnodes = concat nodesList
  let taxnamednodes = map (setNodeScientificName filteredNodeNames) taxnodes
  return (mkGraph taxnamednodes taxedges)

setNodeScientificName :: [TaxName] -> (t, SimpleTaxon) -> (t, SimpleTaxon)
setNodeScientificName inputTaxNames (inputNode,inputTaxon) = outputNode
  where maybeRetrievedName = find (\a -> nameTaxId a == simpleTaxId inputTaxon) inputTaxNames
        retrievedName = maybe "no name" nameTxt maybeRetrievedName
        outputNode = (inputNode,inputTaxon{simpleScientificName = retrievedName})

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
      
-- | Extract a subtree correpsonding to input node paths to root. If a Rank is provided, all node that are less or equal are omitted
extractTaxonomySubTree :: [Node] -> (Gr SimpleTaxon Double) -> Maybe Rank -> (Gr SimpleTaxon Double)
extractTaxonomySubTree inputNodes graph highestRank = taxonomySubTree
  where paths = nub (concatMap (\n -> (sp (n :: Node) (1 :: Node) graph)) inputNodes)
        contexts = map (context graph) paths
        lnodes = map labNode' contexts
        filteredLNodes = filterNodesByRank highestRank lnodes
        ledges = nub (concatMap (out graph) (map fst filteredLNodes))
        taxonomySubTree = (mkGraph filteredLNodes ledges) :: (Gr SimpleTaxon Double)

filterNodesByRank :: Maybe Rank -> [(t, SimpleTaxon)] -> [(t, SimpleTaxon)]
filterNodesByRank highestRank inputNodes
  | (isJust highestRank) = filter (\(_,t) -> simpleRank t >= (fromJust highestRank)) inputNodes
  | otherwise = inputNodes
        
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

constructTaxTree :: [TaxNode] -> Tree TaxNode
constructTaxTree (taxnode:taxnodes) = Node taxnode (concat (addChildElements (taxId taxnode) taxnodes))

addChildElements :: Int -> [TaxNode] -> [[Tree TaxNode]]
addChildElements currentTaxId taxnodes = do
  let (childElements, remainingElements) = partition (\x -> parentTaxId x == currentTaxId) taxnodes
  let subtreeLists = map (\x -> (x:remainingElements)) childElements
  let subtrees = constructSubTrees subtreeLists
  return subtrees
         
constructSubTrees :: [[TaxNode]] -> [Tree TaxNode]
constructSubTrees subtreeLists =  map constructTaxTree subtreeLists

-- | parse NCBITaxCitations from input string
parseNCBITaxCitations :: [Char] -> Either ParseError [TaxCitation]
parseNCBITaxCitations input = parse genParserNCBITaxCitations "parseTaxCitations" input

-- | parse NCBITaxCitations from input filePath                      
readNCBITaxCitations :: String -> IO (Either ParseError [TaxCitation])  
readNCBITaxCitations filePath = parseFromFileEncISO88591 genParserNCBITaxCitations filePath

-- | parse NCBITaxDelNodes from input string
parseNCBITaxDelNodes :: [Char] -> Either ParseError [TaxDelNode]
parseNCBITaxDelNodes input = parse genParserNCBITaxDelNodes "parseTaxDelNodes" input

-- | parse NCBITaxDelNodes from input filePath                      
readNCBITaxDelNodes :: String -> IO (Either ParseError [TaxDelNode])  
readNCBITaxDelNodes filePath = parseFromFile genParserNCBITaxDelNodes filePath

-- | parse NCBITaxDivisons from input string
parseNCBITaxDivisions :: [Char] -> Either ParseError [TaxDivision]
parseNCBITaxDivisions input = parse genParserNCBITaxDivisons "parseTaxDivisons" input

-- | parse NCBITaxDivisons from input filePath                      
readNCBITaxDivisions :: String -> IO (Either ParseError [TaxDivision])  
readNCBITaxDivisions filePath = parseFromFile genParserNCBITaxDivisons filePath

-- | parse NCBITaxGenCodes from input string
parseNCBITaxGenCodes :: [Char] -> Either ParseError [TaxGenCode]
parseNCBITaxGenCodes input = parse genParserNCBITaxGenCodes "parseTaxGenCodes" input

-- | parse NCBITaxGenCodes from input filePath                      
readNCBITaxGenCodes :: String -> IO (Either ParseError [TaxGenCode])  
readNCBITaxGenCodes filePath = parseFromFile genParserNCBITaxGenCodes filePath

-- | parse NCBITaxMergedNodes from input string
parseNCBITaxMergedNodes :: [Char] -> Either ParseError [TaxMergedNode]
parseNCBITaxMergedNodes input = parse genParserNCBITaxMergedNodes "parseTaxMergedNodes" input

-- | parse NCBITaxMergedNodes from input filePath                      
readNCBITaxMergedNodes :: String -> IO (Either ParseError [TaxMergedNode])  
readNCBITaxMergedNodes filePath = parseFromFile genParserNCBITaxMergedNodes filePath

-- | parse NCBITaxNames from input string
parseNCBITaxNames :: [Char] -> Either ParseError [TaxName]
parseNCBITaxNames input = parse genParserNCBITaxNames "parseTaxNames" input

-- | parse NCBITaxNames from input filePath                      
readNCBITaxNames :: String -> IO (Either ParseError [TaxName])  
readNCBITaxNames filePath = parseFromFile genParserNCBITaxNames filePath

-- | parse NCBITaxNames from input string
parseNCBITaxNodes :: [Char] -> Either ParseError TaxNode
parseNCBITaxNodes input = parse genParserNCBITaxNode "parseTaxNode" input

-- | parse NCBITaxCitations from input filePath                      
readNCBITaxNodes :: String -> IO (Either ParseError [TaxNode])  
readNCBITaxNodes filePath = parseFromFile genParserNCBITaxNodes filePath

-- | parse NCBISimpleTaxNames from input string
parseNCBISimpleTaxons :: [Char] -> Either ParseError SimpleTaxon
parseNCBISimpleTaxons input = parse genParserNCBISimpleTaxon "parseSimpleTaxon" input

-- | parse NCBITaxCitations from input filePath                      
readNCBISimpleTaxons :: String -> IO (Either ParseError [SimpleTaxon])  
readNCBISimpleTaxons filePath = parseFromFile genParserNCBISimpleTaxons filePath

-- | Parse the input as NCBITax datatype
readNCBITaxonomyDatabase :: String -> IO (Either [[Char]] NCBITaxDump)
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
  return $ (checkParsing parseErrors citations taxdelNodes divisons genCodes mergedNodes names taxnodes)

genParserNCBITaxCitations :: GenParser Char st [TaxCitation]
genParserNCBITaxCitations = do
  citations <- many1 genParserNCBITaxCitation
  return $ citations

genParserNCBITaxDelNodes :: GenParser Char st [TaxDelNode]
genParserNCBITaxDelNodes = do
  taxdelNodes <- many1 genParserNCBITaxDelNode
  return $ taxdelNodes
  
genParserNCBITaxDivisons :: GenParser Char st [TaxDivision]
genParserNCBITaxDivisons = do
  divisions <- many1 genParserNCBITaxDivision
  return $ divisions

genParserNCBITaxGenCodes :: GenParser Char st [TaxGenCode]
genParserNCBITaxGenCodes = do
  genCodes <- many1 genParserNCBITaxGenCode
  return $ genCodes

genParserNCBITaxMergedNodes :: GenParser Char st [TaxMergedNode]
genParserNCBITaxMergedNodes = do
  mergedNodes <- many1 genParserNCBITaxMergedNode
  return $ mergedNodes

genParserNCBITaxNames :: GenParser Char st [TaxName]
genParserNCBITaxNames = do
  names <- many1 genParserNCBITaxName
  return $ names

genParserNCBITaxNodes :: GenParser Char st [TaxNode]
genParserNCBITaxNodes = do
  taxnodes <- many1 genParserNCBITaxNode
  return $ taxnodes

genParserNCBISimpleTaxons :: GenParser Char st [SimpleTaxon]
genParserNCBISimpleTaxons = do
  taxnodes <- many1 genParserNCBISimpleTaxon
  return $ taxnodes
----------------------------

genParserNCBITaxCitation :: GenParser Char st TaxCitation
genParserNCBITaxCitation = do
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
  return $ TaxCitation (readInt _citId) _citKey (liftM readInt _pubmedId) (liftM readInt _medlineId) _url _text _taxIdList

genParserNCBITaxDelNode :: GenParser Char st TaxDelNode
genParserNCBITaxDelNode = do
  taxdelNode <- many1 digit
  space
  char ('|')
  char ('\n')
  return $ TaxDelNode (readInt taxdelNode)
  
genParserNCBITaxDivision :: GenParser Char st TaxDivision
genParserNCBITaxDivision = do
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
  return $ TaxDivision (readInt _divisionId) _divisionCDE _divisionName _comments 

genParserNCBITaxGenCode :: GenParser Char st TaxGenCode
genParserNCBITaxGenCode = do
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
  return $ TaxGenCode (readInt _geneticCodeId) _abbreviation _genCodeName _cde _starts

genParserNCBITaxMergedNode :: GenParser Char st TaxMergedNode
genParserNCBITaxMergedNode = do
  _oldTaxId <- many1 digit
  tab
  char ('|')
  tab  
  _newTaxId <- many1 digit
  tab
  char ('|')
  char ('\n')
  return $ TaxMergedNode (readInt _oldTaxId) (readInt _newTaxId)

genParserNCBITaxName :: GenParser Char st TaxName
genParserNCBITaxName = do
  _taxId <- many1 digit
  tab
  char ('|')
  tab
  _nameTxt <- many1 (noneOf ("\t\n"))
  tab
  char ('|')
  tab
  _uniqueName <- optionMaybe (many1 (noneOf "\t\n"))
  tab
  char ('|')
  tab
  _nameClass <- many1 (noneOf ("\t\n"))
  tab
  char ('|')
  newline
  return $ TaxName (readInt _taxId) _nameTxt _uniqueName _nameClass

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

genParserNCBITaxNode :: GenParser Char st TaxNode
genParserNCBITaxNode = do
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
  return $ TaxNode (readInt _taxId) (readInt _parentTaxId) (readRank _rank) _emblCode _divisionId (readBool _inheritedDivFlag) _geneticCodeId (readBool _inheritedGCFlag) _mitochondrialGeneticCodeId (readBool _inheritedMGCFlag) (readBool _genBankHiddenFlag) (readBool _hiddenSubtreeRootFlag) _comments

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
checkParsing :: [[Char]] -> Either ParseError [TaxCitation] -> Either ParseError [TaxDelNode] -> Either ParseError [TaxDivision] -> Either ParseError [TaxGenCode] -> Either ParseError [TaxMergedNode] -> Either ParseError [TaxName] -> Either ParseError [TaxNode]-> Either [[Char]] NCBITaxDump
checkParsing parseErrors citations taxdelNodes divisons genCodes mergedNodes names taxnodes
  | join (parseErrors) == "" = Right (NCBITaxDump (E.fromRight citations) (E.fromRight taxdelNodes) (E.fromRight divisons) (E.fromRight genCodes) (E.fromRight mergedNodes) (E.fromRight names) (E.fromRight taxnodes))
  | otherwise = Left (parseErrors)

extractParseError :: Either ParseError a -> String
extractParseError _parse
  | isLeft _parse = show (E.fromLeft _parse)
  | otherwise = ""
