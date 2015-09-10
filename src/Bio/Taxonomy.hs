-- | Functions for parsing, processing and visualization of taxonomy data.
--
-- === Usage example:
-- * Read in taxonomy data
--   
-- * Process data
--
-- * Visualize result
--
module Bio.Taxonomy (  -- * Datatypes
                       -- Datatypes used to represent taxonomy data
                       module Bio.TaxonomyData,
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
                       readNCBITaxonomyDatabase,
                       -- * Processing
                       compareSubTrees,    
                       extractTaxonomySubTreebyLevel,
                       extractTaxonomySubTreebyRank,
                       -- * Visualization
                       getParentbyRank,
                       drawTreeComparison,
                       drawTaxonomy
                      ) where
import Prelude 
import System.IO 
import Bio.TaxonomyData
import Text.Parsec.Prim (runP)
import Text.ParserCombinators.Parsec
import Control.Monad
import Data.List
import qualified Data.Vector as V
import Data.Maybe    
import Data.Either
import qualified Data.Either.Unwrap as E
import Data.Graph.Inductive
import qualified Data.GraphViz as GV
import qualified Data.GraphViz.Printing as GVP
import qualified Data.GraphViz.Attributes.Colors as GVAC
import qualified Data.GraphViz.Attributes.Complete as GVA
import qualified Data.Text.Lazy as TL
import qualified Data.ByteString.Char8 as B
--------------------------------------------------------
-- *  Visualization


-- | Draw graph in dot format
drawTaxonomy :: Gr SimpleTaxon Double -> String
drawTaxonomy inputGraph = do
  let params = GV.nonClusteredParams {GV.isDirected       = True
                       , GV.globalAttributes = [GV.GraphAttrs [GVA.Size (GVA.GSize (20 :: Double) (Just (20 :: Double)) False)]]
                       , GV.isDotCluster     = const True
                       , GV.fmtNode = \ (_,l) -> [GV.textLabel (TL.pack (show (simpleRank l) ++ "\n" ++ B.unpack (simpleScientificName l)))]
                       , GV.fmtEdge          = const []
                       }
  let dotFormat = GV.graphToDot params inputGraph
  let dottext = GVP.renderDot $ GVP.toDot dotFormat
  TL.unpack dottext

-- | Draw tree comparison graph in dot format
drawTreeComparison :: (Int,Gr CompareTaxon Double) -> String
drawTreeComparison (treeNumber,inputGraph) = do
  let cList = makeColorList treeNumber 
  let params = GV.nonClusteredParams {GV.isDirected = True
                      , GV.globalAttributes = []
                       , GV.isDotCluster = const True
                       , GV.fmtNode = \ (_,l) -> [GV.textLabel (TL.pack (show (compareRank l) ++ "\n" ++ B.unpack (compareScientificName l))), GV.style GV.wedged, GVA.Color (selectColors (inTree l) cList)]
                       , GV.fmtEdge = const []
                       }
  let dotFormat = GV.graphToDot params (grev inputGraph)
  let dottext = GVP.renderDot $ GVP.toDot dotFormat
  TL.unpack dottext

-- | Colors from color list are selected according to in which of the compared trees the node is contained.
selectColors :: [Int] -> [GVA.Color] -> GVAC.ColorList
selectColors inTrees currentColorList = GVAC.toColorList (map (\i -> currentColorList !! i) inTrees)

-- | A color list is sampled from the spectrum according to how many trees are compared.
makeColorList :: Int -> [GVA.Color]
makeColorList treeNumber = cList
  where cList = map (\i -> GVAC.HSV ((fromIntegral i/fromIntegral neededColors) * 0.708) 0.5 1.0)  [0..neededColors]
        neededColors = treeNumber - 1

-- | NCBI taxonomy dump nodes and names in the input directory path are parsed and a SimpleTaxon tree with scientific names for each node is generated.  
readNamedTaxonomy :: String -> IO (Either ParseError (Gr SimpleTaxon Double))  
readNamedTaxonomy directoryPath = do
  nodeNames <- readNCBITaxNames (directoryPath ++ "names.dmp")
  let scientificNameBS = B.pack "scientific name"
  if isLeft nodeNames
     then return (Left (E.fromLeft nodeNames))
     else do
       let nodeNamesVector = V.fromList (E.fromRight nodeNames)
       let filteredNodeNames = V.filter (\a -> nameClass a == scientificNameBS) nodeNamesVector
       parseFromFileEncISO88591 (genParserNamedTaxonomyGraph filteredNodeNames) (directoryPath ++ "nodes.dmp")

-- * Functions for parsing taxonomy

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
  let taxedges = filter (\(a,b,_) -> a /= b) edgesList
  --let taxnodes = concat nodesList
  --return (mkGraph taxnodes taxedges)
  return $! mkGraph nodesList taxedges

genParserNamedTaxonomyGraph :: V.Vector TaxName -> GenParser Char st (Gr SimpleTaxon Double)
genParserNamedTaxonomyGraph filteredNodeNames = do
  nodesEdges <- many1 (try genParserGraphNodeEdge)
  optional eof
  let (nodesList,edgesList) = unzip nodesEdges
  let taxedges = filter (\(a,b,_) -> a /= b) edgesList
  let taxnamednodes = map (setNodeScientificName filteredNodeNames) nodesList
  return $! mkGraph taxnamednodes taxedges

setNodeScientificName :: V.Vector TaxName -> (t, SimpleTaxon) -> (t, SimpleTaxon)
setNodeScientificName inputTaxNames (inputNode,inputTaxon) = outputNode
  where maybeRetrievedName = V.find (\a -> nameTaxId a == simpleTaxId inputTaxon) inputTaxNames
        retrievedName = maybe (B.pack "no name") nameTxt maybeRetrievedName
        outputNode = (inputNode,inputTaxon{simpleScientificName = retrievedName})

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
  return ((_simpleTaxIdInt,SimpleTaxon _simpleTaxIdInt B.empty _simpleParentTaxIdInt (readRank _simpleRank)),(_simpleTaxIdInt,_simpleParentTaxIdInt,1 :: Double))
      
-- | Extract a subtree correpsonding to input node paths to root. Only nodes in level number distance to root are included
compareSubTrees :: [Gr SimpleTaxon Double] -> (Int,Gr CompareTaxon Double)
compareSubTrees graphs = (length graphs,resultGraph)
  where treesLabNodes = map labNodes graphs
        treesLabEdges = map labEdges graphs
        mergedNodes = nub (concat treesLabNodes)
        mergedEdges = nub (concat treesLabEdges)
        --annotate node in which of the compared trees they are present
        comparedNodes = annotateTaxonsDifference treesLabNodes mergedNodes
        resultGraph = mkGraph comparedNodes mergedEdges :: Gr CompareTaxon Double

annotateTaxonsDifference  :: [[LNode SimpleTaxon]] -> [LNode SimpleTaxon] -> [LNode CompareTaxon]
annotateTaxonsDifference  treesNodes mergedtreeNodes = comparedNodes
  where comparedNodes = map (annotateTaxonDifference indexedTreesNodes) mergedtreeNodes
        indexedTreesNodes = zip [0..(length treesNodes)] treesNodes
        

annotateTaxonDifference :: [(Int,[LNode SimpleTaxon])] -> LNode SimpleTaxon -> LNode CompareTaxon
annotateTaxonDifference indexedTreesNodes mergedtreeNode = comparedNode
  where comparedNode = (simpleTaxId (snd mergedtreeNode),CompareTaxon (simpleScientificName (snd mergedtreeNode)) (simpleRank (snd mergedtreeNode)) currentInTree)
        currentInTree = concatMap (\(i,treeNodes) -> [i | mergedtreeNode `elem` treeNodes]) indexedTreesNodes
        
-- | Extract a subtree corresponding to input node paths to root. Only nodes in level number distance to root are included
extractTaxonomySubTreebyLevel :: [Node] -> Gr SimpleTaxon Double -> Maybe Int -> Gr SimpleTaxon Double
extractTaxonomySubTreebyLevel inputNodes graph levelNumber = taxonomySubTree
  where paths = nub (concatMap (getPath (1 :: Node) graph) inputNodes)
        contexts = map (context graph) paths
        lnodes = map labNode' contexts
        ledges = nub (concatMap (out graph . fst) lnodes)
        unfilteredTaxonomySubTree = mkGraph lnodes ledges :: Gr SimpleTaxon Double
        filteredLNodes = filterNodesByLevel levelNumber lnodes unfilteredTaxonomySubTree
        filteredledges = nub (concatMap (out graph . fst) filteredLNodes)
        taxonomySubTree = mkGraph filteredLNodes filteredledges :: Gr SimpleTaxon Double      

-- | Extract a subtree corresponding to input node paths to root. If a Rank is provided, all node that are less or equal are omitted
extractTaxonomySubTreebyRank :: [Node] -> Gr SimpleTaxon Double -> Maybe Rank -> Gr SimpleTaxon Double
extractTaxonomySubTreebyRank inputNodes graph highestRank = taxonomySubTree
  where paths = nub (concatMap (getPath (1 :: Node) graph) inputNodes)
        contexts = map (context graph) paths
        lnodes = map labNode' contexts
        filteredLNodes = filterNodesByRank highestRank lnodes
        filteredledges = nub (concatMap (out graph . fst) filteredLNodes)
        taxonomySubTree = mkGraph filteredLNodes filteredledges :: Gr SimpleTaxon Double

getPath :: Node -> Gr SimpleTaxon Double -> Node -> Path
getPath root graph node =  sp node root graph

-- | Extract parent node iwth specified Rank
getParentbyRank :: Node -> Gr SimpleTaxon Double -> Maybe Rank -> Maybe (Node, SimpleTaxon)
getParentbyRank inputNode graph requestedRank = filteredLNode
  where path = sp (inputNode :: Node) (1 :: Node) graph
        nodeContext = map (context graph) path
        lnode = map labNode' nodeContext
        filteredLNode = findNodeByRank requestedRank lnode
           
filterNodesByLevel :: Maybe Int -> [(Node, SimpleTaxon)] -> Gr SimpleTaxon Double -> [(Node, SimpleTaxon)]
filterNodesByLevel levelNumber inputNodes graph
  | isJust levelNumber = filteredNodes
  | otherwise = inputNodes
    --distances of all nodes to root
    where nodedistances = level (1::Node) (undir graph)
          sortedNodeDistances = sortBy sortByNodeID nodedistances
          sortedInputNodes = sortBy sortByNodeID inputNodes
          zippedNodeDistancesInputNodes = zip sortedNodeDistances sortedInputNodes
          zippedFilteredNodes = filter (\((_,distance),(_,_)) -> distance <= fromJust levelNumber) zippedNodeDistancesInputNodes
          filteredNodes = map snd zippedFilteredNodes

sortByNodeID :: (Node,a) -> (Node,a) -> Ordering
sortByNodeID (n1, _) (n2, _)
  | n1 < n2 = GT
  | n1 > n2 = LT
  | n1 == n2 = EQ
  | otherwise = EQ

findNodeByRank :: Maybe Rank -> [(t, SimpleTaxon)] -> Maybe (t, SimpleTaxon)
findNodeByRank requestedRank inputNodes
  | isJust requestedRank = filteredNodes
  | otherwise = Nothing
    where filteredNodes = find (\(_,t) -> simpleRank t == fromJust requestedRank) inputNodes


filterNodesByRank :: Maybe Rank -> [(t, SimpleTaxon)] -> [(t, SimpleTaxon)]
filterNodesByRank highestRank inputNodes
  | isJust highestRank = filteredNodes
  | otherwise = inputNodes
    where filteredNodes = filter (\(_,t) -> simpleRank t >= fromJust highestRank) inputNodes ++ noRankNodes
          noRankNodes = filter (\(_,t) -> simpleRank t == Norank) inputNodes

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
  _citKey <- optionMaybe (many1 (noneOf "\t"))
  string "\t|\t"
  _pubmedId <- optionMaybe (many1 digit)
  string "\t|\t"
  _medlineId <- optionMaybe (many1 digit)
  tab
  char '|' 
  _url <- genParserTaxURL
  char '|'
  tab
  _text <- optionMaybe (many1 (noneOf "\t"))
  string "\t|\t"
  _taxIdList <- optionMaybe (many1 genParserTaxIdList)
  string "\t|\n"
  return $ TaxCitation (readInt _citId) _citKey (liftM readInt _pubmedId) (liftM readInt _medlineId) _url _text _taxIdList

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
  _comments <- optionMaybe (many1 (noneOf "\t"))
  string "\t|\n"
  return $ TaxDivision (readInt _divisionId) _divisionCDE _divisionName _comments 

genParserNCBITaxGenCode :: GenParser Char st TaxGenCode
genParserNCBITaxGenCode = do
  _geneticCodeId <- many1 digit 
  string "\t|\t"
  _abbreviation <- optionMaybe (many1 (noneOf "\t"))
  string "\t|\t"
  _genCodeName <- many1 (noneOf "\t")
  string "\t|\t"
  _cde <- many1 (noneOf "\t")
  string "\t|\t"
  _starts <- many1 (noneOf "\t")
  string "\t|\n"
  return $ TaxGenCode (readInt _geneticCodeId) _abbreviation _genCodeName _cde _starts

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
  _uniqueName <- optionMaybe (many1 (noneOf "\t\n"))
  string "\t|\t"
  _nameClass <- many1 (noneOf "\t\n")
  tab
  char '|'
  newline
  return $! TaxName (readInt _taxId) (B.pack _nameTxt) (maybe B.empty B.pack _uniqueName) (B.pack _nameClass)

genParserNCBISimpleTaxon :: GenParser Char st SimpleTaxon
genParserNCBISimpleTaxon = do
  _simpleTaxId <- many1 digit
  string "\t|\t"
  _simpleParentTaxId <- many1 digit
  string "\t|\t"
  _simpleRank <- many1 (noneOf "\t")
  many1 (noneOf "\n")
  char '\n'
  return $! SimpleTaxon (readInt _simpleTaxId) B.empty (readInt _simpleParentTaxId) (readRank _simpleRank) 

genParserNCBITaxNode :: GenParser Char st TaxNode
genParserNCBITaxNode = do
  _taxId <- many1 digit
  string "\t|\t"
  _parentTaxId <- many1 digit
  string "\t|\t"
  _rank <- many1 (noneOf "\t")
  string "\t|\t"
  _emblCode <- optionMaybe (many1 (noneOf "\t"))
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
  _comments <- optionMaybe (many1 (noneOf "\t"))
  tab
  char '|'
  char '\n'
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
  return (readInt _taxId)

genParserTaxURL :: GenParser Char st (Maybe String)
genParserTaxURL = do
  tab 
  url1 <- optionMaybe (many1 (noneOf "\t"))
  tab
  url2 <- optionMaybe (many1 (noneOf "|"))
  return (concatenateURLParts url1 url2)

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
  | isLeft _parse = show (E.fromLeft _parse)
  | otherwise = ""
