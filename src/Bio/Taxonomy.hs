-- | Parse and process taxonomy data

module Bio.Taxonomy (                      
                       module Bio.TaxonomyData,
                       getParentbyRank,
---                       drawTreeComparison,
                       compareSubTrees,    
                       extractTaxonomySubTreebyLevel,
                       extractTaxonomySubTreebyRank,
---                      drawTaxonomy,    
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
---import qualified Data.GraphViz as GV
---import qualified Data.GraphViz.Printing as GVP
---import qualified Data.GraphViz.Attributes.Colors as GVAC
---import qualified Data.GraphViz.Attributes.Complete as GVA
import qualified Data.Text.Lazy as TL
import qualified Data.ByteString.Char8 as B
--------------------------------------------------------

--fgl graph representation

-- | draw Graph in dot format
---drawTaxonomy :: Gr SimpleTaxon Double -> String
---drawTaxonomy inputGraph = do
---  let params = GV.nonClusteredParams {GV.isDirected       = True
---                       , GV.globalAttributes = [GV.GraphAttrs [GVA.Size (GVA.GSize (20 :: Double) (Just (20 :: Double)) False)]]
---                       , GV.isDotCluster     = const True
---                       , GV.fmtNode = \ (_,l) -> [GV.textLabel (TL.pack ((show (simpleRank l)) ++ "\n" ++ (B.unpack (simpleScientificName l))))]
---                       , GV.fmtEdge          = const []
---                       }
---  let dotFormat = GV.graphToDot params inputGraph
---  let dottext = GVP.renderDot $ GVP.toDot dotFormat
---  TL.unpack dottext

-- | draw Comparison graph in dot format
---drawTreeComparison :: (Int,(Gr CompareTaxon Double)) -> String
---drawTreeComparison (treeNumber,inputGraph) = do
---  let cList = makeColorList treeNumber 
---  let params = GV.nonClusteredParams {GV.isDirected = True
---                      , GV.globalAttributes = []
---                       , GV.isDotCluster = const True
---                       , GV.fmtNode = \ (_,l) -> [GV.textLabel (TL.pack ((show (compareRank l)) ++ "\n" ++ (B.unpack (compareScientificName l)))), GV.style GV.wedged, GVA.Color (selectColors (inTree l) cList)]
---                       , GV.fmtEdge = const []
---                       }
---  let dotFormat = GV.graphToDot params (grev inputGraph)
---  let dottext = GVP.renderDot $ GVP.toDot dotFormat
---  TL.unpack dottext

---selectColors :: [Int] -> [GVA.Color] -> GVAC.ColorList
---selectColors inTrees currentColorList = GVAC.toColorList (map (\i -> currentColorList !! i) inTrees)

---makeColorList :: Int -> [GVA.Color]
---makeColorList treeNumber = cList
---  where cList = map (\i -> GVAC.HSV (((fromIntegral i)/(fromIntegral neededColors)) * 0.708) 0.5 1.0)  [0..neededColors]
---        neededColors = treeNumber - 1

-- | parse Taxonomy from input filePath                      
readNamedTaxonomy :: String -> IO (Either ParseError (Gr SimpleTaxon Double))  
readNamedTaxonomy directoryPath = do
  nodeNames <- readNCBITaxNames (directoryPath ++ "names.dmp")
  let scientificNameBS = B.pack ("scientific name")
  if (isLeft nodeNames)
     then do 
       return (Left (E.fromLeft nodeNames))
     else do
       let nodeNamesVector = V.fromList (E.fromRight nodeNames)
       let filteredNodeNames = V.filter (\a -> nameClass a == scientificNameBS) nodeNamesVector
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
  let taxedges = filter (\(a,b,_) -> a /= b) edgesList
  --let taxnodes = concat nodesList
  --return (mkGraph taxnodes taxedges)
  return (mkGraph nodesList taxedges)

genParserNamedTaxonomyGraph :: V.Vector TaxName -> GenParser Char st (Gr SimpleTaxon Double)
genParserNamedTaxonomyGraph filteredNodeNames = do
  nodesEdges <- many1 (try (genParserGraphNodeEdge))
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
  many1 (noneOf "\n")
  char ('\n')
  let _simpleTaxIdInt = readInt _simpleTaxId
  let _simpleParentTaxIdInt = readInt _simpleParentTaxId
  return $! ((_simpleTaxIdInt,SimpleTaxon _simpleTaxIdInt B.empty _simpleParentTaxIdInt (readRank _simpleRank)),(_simpleTaxIdInt,_simpleParentTaxIdInt,(1 :: Double)))
      
-- | Extract a subtree correpsonding to input node paths to root. Only nodes in level number distance to root are included
compareSubTrees :: [(Gr SimpleTaxon Double)] -> (Int,(Gr CompareTaxon Double))
compareSubTrees graphs = (length graphs,resultGraph)
  where treesLabNodes = map labNodes graphs
        treesLabEdges = map labEdges graphs
        mergedNodes = nub (concat treesLabNodes)
        mergedEdges = nub (concat treesLabEdges)
        --annotate node in which of the compared trees they are present
        comparedNodes = annotateTaxonsDifference treesLabNodes mergedNodes
        resultGraph = (mkGraph comparedNodes mergedEdges) :: (Gr CompareTaxon Double)

annotateTaxonsDifference  :: [[LNode SimpleTaxon]] -> [LNode SimpleTaxon] -> [LNode CompareTaxon]
annotateTaxonsDifference  treesNodes mergedtreeNodes = comparedNodes
  where comparedNodes = map (annotateTaxonDifference indexedTreesNodes) mergedtreeNodes
        indexedTreesNodes = zip [0..(length treesNodes)] treesNodes
        

annotateTaxonDifference :: [(Int,[LNode SimpleTaxon])] -> LNode SimpleTaxon -> LNode CompareTaxon
annotateTaxonDifference indexedTreesNodes mergedtreeNode = comparedNode
  where comparedNode = ((simpleTaxId (snd mergedtreeNode)),(CompareTaxon (simpleScientificName (snd mergedtreeNode)) (simpleRank (snd mergedtreeNode)) currentInTree))
        currentInTree = concatMap (\(i,treeNodes) -> if (elem mergedtreeNode treeNodes) then [i] else []) indexedTreesNodes
        
-- | Extract a subtree corresponding to input node paths to root. Only nodes in level number distance to root are included
extractTaxonomySubTreebyLevel :: [Node] -> (Gr SimpleTaxon Double) -> Maybe Int -> (Gr SimpleTaxon Double)
extractTaxonomySubTreebyLevel inputNodes graph levelNumber = taxonomySubTree
  where paths = nub (concatMap (\n -> (sp (n :: Node) (1 :: Node) graph)) inputNodes)
        contexts = map (context graph) paths
        lnodes = map labNode' contexts
        ledges = nub (concatMap (out graph) (map fst lnodes))
        unfilteredTaxonomySubTree = (mkGraph lnodes ledges) :: (Gr SimpleTaxon Double)
        filteredLNodes = filterNodesByLevel levelNumber lnodes unfilteredTaxonomySubTree
        filteredledges = nub (concatMap (out graph) (map fst filteredLNodes))
        taxonomySubTree = (mkGraph filteredLNodes  filteredledges) :: (Gr SimpleTaxon Double)                 

-- | Extract a subtree corresponding to input node paths to root. If a Rank is provided, all node that are less or equal are omitted
extractTaxonomySubTreebyRank :: [Node] -> (Gr SimpleTaxon Double) -> Maybe Rank -> (Gr SimpleTaxon Double)
extractTaxonomySubTreebyRank inputNodes graph highestRank = taxonomySubTree
  where paths = nub (concatMap (\n -> (sp (n :: Node) (1 :: Node) graph)) inputNodes)
        contexts = map (context graph) paths
        lnodes = map labNode' contexts
        filteredLNodes = filterNodesByRank highestRank lnodes
        ledges = nub (concatMap (out graph) (map fst lnodes))
        taxonomySubTree = (mkGraph filteredLNodes ledges) :: (Gr SimpleTaxon Double)

-- | Extract parent node iwth specified Rank
getParentbyRank :: Node -> (Gr SimpleTaxon Double) -> Maybe Rank -> Maybe (Node, SimpleTaxon)
getParentbyRank inputNode graph requestedRank = filteredLNode
  where path = sp (inputNode :: Node) (1 :: Node) graph
        nodeContext = map (context graph) path
        lnode = map labNode' nodeContext
        filteredLNode = findNodeByRank requestedRank lnode
           
filterNodesByLevel :: Maybe Int -> [(Node, SimpleTaxon)] -> (Gr SimpleTaxon Double) -> [(Node, SimpleTaxon)]
filterNodesByLevel levelNumber inputNodes graph
  | (isJust levelNumber) = filteredNodes
  | otherwise = inputNodes
    --distances of all nodes to root
    where nodedistances = level (1::Node) (undir graph)
          sortedNodeDistances = sortBy sortByNodeID nodedistances
          sortedInputNodes = sortBy sortByNodeID inputNodes
          zippedNodeDistancesInputNodes = zip sortedNodeDistances sortedInputNodes
          zippedFilteredNodes = filter (\((_,distance),(_,_)) -> distance <= (fromJust levelNumber)) zippedNodeDistancesInputNodes
          filteredNodes = map snd zippedFilteredNodes

sortByNodeID :: (Node,a) -> (Node,a) -> Ordering
sortByNodeID (n1, _) (n2, _)
  | n1 < n2 = GT
  | n1 > n2 = LT
  | n1 == n2 = EQ
  | otherwise = EQ

findNodeByRank :: Maybe Rank -> [(t, SimpleTaxon)] -> Maybe (t, SimpleTaxon)
findNodeByRank requestedRank inputNodes
  | (isJust requestedRank) = filteredNodes
  | otherwise = Nothing
    where filteredNodes = find (\(_,t) -> simpleRank t == (fromJust requestedRank)) inputNodes


filterNodesByRank :: Maybe Rank -> [(t, SimpleTaxon)] -> [(t, SimpleTaxon)]
filterNodesByRank highestRank inputNodes
  | (isJust highestRank) = filteredNodes
  | otherwise = inputNodes
    where filteredNodes = filter (\(_,t) -> simpleRank t >= (fromJust highestRank)) inputNodes ++ noRankNodes
          noRankNodes = filter (\(_,t) -> simpleRank t == Norank) inputNodes

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
  return $! names

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
  return $! TaxName (readInt _taxId) (B.pack _nameTxt) (maybe B.empty B.pack _uniqueName) (B.pack _nameClass)

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
  return $! SimpleTaxon (readInt _simpleTaxId) B.empty (readInt _simpleParentTaxId) (readRank _simpleRank) 

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
