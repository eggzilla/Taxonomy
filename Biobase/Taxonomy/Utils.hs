-- | Functions for processing of taxonomy data.
--
module Biobase.Taxonomy.Utils (  -- * Datatypes
                       -- Datatypes used to represent taxonomy data
                       module Biobase.Taxonomy.Types,
                       -- * Processing
                       compareSubTrees,
                       extractTaxonomySubTreebyLevel,
                       extractTaxonomySubTreebyLevelNew,
                       extractTaxonomySubTreebyRank,
                       safeNodePath,
                       getParentbyRank,
                      ) where
import Prelude
import Biobase.Taxonomy.Types
import Data.List
import qualified Data.Vector as V
import Data.Maybe
import Data.Graph.Inductive.Graph
import Data.Graph.Inductive.Query.SP (sp)
import Data.Graph.Inductive.Query.BFS (level)
import Data.Graph.Inductive.Tree
import Data.Graph.Inductive.Basic

---------------------------------------
-- Processing functions

-- | Extract a subtree correpsonding to input node paths to root. Only nodes in level number distance to root are included. Used in Ids2TreeCompare tool.
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

-- | Extract a subtree corresponding to input node paths to root. Only nodes in level number distance to root are included. Used in Ids2Tree tool.
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

-- | Extract a subtree corresponding to input node paths to root. Only nodes in level number distance to root are included. Used in Ids2Tree tool.
extractTaxonomySubTreebyLevelNew :: [Node] -> Gr SimpleTaxon Double -> Maybe Int -> Gr SimpleTaxon Double
extractTaxonomySubTreebyLevelNew inputNodes graph levelNumber = taxonomySubTree
  where inputNodeVector = V.fromList inputNodes
        paths = V.concatMap (getVectorPath (1 :: Node) graph) inputNodeVector
        contexts = V.map (context graph) paths
        vlnodes = V.map labNode' contexts
        ledges = concatMap (out graph . fst) lnodes
        lnodes = V.toList vlnodes
        --ledges = V.toList vledges
        unfilteredTaxonomySubTree = mkGraph lnodes ledges :: Gr SimpleTaxon Double
        filteredLNodes = filterNodesByLevel levelNumber lnodes unfilteredTaxonomySubTree
        --filteredLNodesVector = V.fromList filteredLNodes
        filteredledges = concatMap (out graph . fst) filteredLNodes
        --filteredledges = V.toList filteredledgesVector
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

getVectorPath :: Node -> Gr SimpleTaxon Double -> Node -> V.Vector Node
getVectorPath root graph node =  maybe V.empty V.fromList (sp node root graph)

getPath :: Node -> Gr SimpleTaxon Double -> Node -> Path
getPath root graph node =  maybe [] id (sp node root graph)

-- | Extract parent node with specified Rank
getParentbyRank :: Node -> Gr SimpleTaxon Double -> Maybe Rank -> Maybe (Node, SimpleTaxon)
getParentbyRank inputNode graph requestedRank = filteredLNode
  where path =  maybe [] id (sp (inputNode :: Node) (1 :: Node) graph)
        nodeContext = map (context graph) path
        lnode = map labNode' nodeContext
        filteredLNode = findNodeByRank requestedRank lnode

-- | Filter nodes by distance from root
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

-- | Find only taxons of a specific rank in a list of input taxons
findNodeByRank :: Maybe Rank -> [(t, SimpleTaxon)] -> Maybe (t, SimpleTaxon)
findNodeByRank requestedRank inputNodes
  | isJust requestedRank = filteredNodes
  | otherwise = Nothing
    where filteredNodes = find (\(_,t) -> simpleRank t == fromJust requestedRank) inputNodes

-- | Filter a list of input taxons for a minimal provided rank
filterNodesByRank :: Maybe Rank -> [(t, SimpleTaxon)] -> [(t, SimpleTaxon)]
filterNodesByRank highestRank inputNodes
  | isJust highestRank = filteredNodes
  | otherwise = inputNodes
    where filteredNodes = filter (\(_,t) -> simpleRank t >= fromJust highestRank) inputNodes ++ noRankNodes
          noRankNodes = filter (\(_,t) -> simpleRank t == Norank) inputNodes

-- | Returns path between 2 maybe nodes. Used in TreeDistance tool.
safeNodePath :: Maybe Node -> Gr SimpleTaxon Double -> Maybe Node -> Either String Path
safeNodePath nodeid1 graphOutput nodeid2
  | isJust nodeid1 && isJust nodeid2 = Right  (maybe [] id (sp (fromJust nodeid1) (fromJust nodeid2) (undir graphOutput)))
  | otherwise = Left "Both taxonomy ids must be provided for distance computation"
