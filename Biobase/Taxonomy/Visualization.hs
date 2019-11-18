-- | Functions for  visualization of taxonomy data.
module Biobase.Taxonomy.Visualization (  -- * Datatypes
                       -- Datatypes used to represent taxonomy data
                       module Biobase.Taxonomy.Types,
                       -- * Visualization
                       drawTaxonomyComparison,
                       drawTaxonomy,
                       writeTree,
                       writeDotTree,
                       writeJsonTree
                      ) where
import Prelude
import Biobase.Taxonomy.Types
import Data.Graph.Inductive.Tree
import Data.Graph.Inductive.Basic
import qualified Data.GraphViz as GV
import qualified Data.GraphViz.Printing as GVP
import qualified Data.GraphViz.Attributes.Colors as GVAC
import qualified Data.GraphViz.Attributes.Complete as GVA
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Aeson as AE
import qualified Data.Text.Lazy as T

---------------------------------------
-- Visualisation functions

-- | Draw graph in dot format. Used in Ids2Tree tool.
drawTaxonomy :: Bool -> Gr SimpleTaxon Double -> String
drawTaxonomy withRank inputGraph = do
  let nodeFormating = if withRank then nodeFormatWithRank else nodeFormatWithoutRank
  let params = GV.nonClusteredParams {GV.isDirected       = True
                       , GV.globalAttributes = [GV.GraphAttrs [GVA.Size (GVA.GSize (20 :: Double) (Just (20 :: Double)) False)]]
                       , GV.isDotCluster     = const True
                       --, GV.fmtNode = \ (_,l) -> [GV.textLabel (TL.pack (show (simpleRank l) ++ "\n" ++ T.unpack (simpleScientificName l)))]
                       , GV.fmtNode = nodeFormating
                       , GV.fmtEdge          = const []
                       }
  let dotFormat = GV.graphToDot params inputGraph
  let dottext = GVP.renderDot $ GVP.toDot dotFormat
  T.unpack dottext

nodeFormatWithRank :: (t, SimpleTaxon) -> [GVA.Attribute]
nodeFormatWithRank (_,l) = [GV.textLabel (T.concat [T.pack (show (simpleRank l)), T.pack ("\n") , simpleScientificName l])]

nodeFormatWithoutRank :: (t, SimpleTaxon) -> [GVA.Attribute]
nodeFormatWithoutRank (_,l) = [GV.textLabel (simpleScientificName l)]

-- | Draw tree comparison graph in dot format. Used in Ids2TreeCompare tool.
drawTaxonomyComparison :: Bool -> (Int,Gr CompareTaxon Double) -> String
drawTaxonomyComparison withRank (treeNumber,inputGraph) = do
  let cList = makeColorList treeNumber
  let nodeFormating = if withRank then (compareNodeFormatWithRank cList) else (compareNodeFormatWithoutRank cList)
  let params = GV.nonClusteredParams {GV.isDirected = True
                       , GV.globalAttributes = []
                       , GV.isDotCluster = const True
                       --, GV.fmtNode = \ (_,l) -> [GV.textLabel (TL.pack (show (compareRank l) ++ "\n" ++ B.unpack (compareScientificName l))), GV.style GV.wedged, GVA.Color (selectColors (inTree l) cList)]
                       , GV.fmtNode = nodeFormating
                       , GV.fmtEdge = const []
                       }
  let dotFormat = GV.graphToDot params (grev inputGraph)
  let dottext = GVP.renderDot $ GVP.toDot dotFormat
  T.unpack dottext

compareNodeFormatWithRank :: [GVA.Color] -> (t, CompareTaxon) -> [GVA.Attribute]
compareNodeFormatWithRank cList (_,l) = [GV.textLabel (T.concat [T.pack (show (compareRank l) ++ "\n"),compareScientificName l]), GV.style GV.wedged, GVA.Color (selectColors (inTree l) cList)]

compareNodeFormatWithoutRank :: [GVA.Color] -> (t, CompareTaxon) -> [GVA.Attribute]
compareNodeFormatWithoutRank cList (_,l) = [GV.textLabel (compareScientificName l), GV.style GV.wedged, GVA.Color (selectColors (inTree l) cList)]

-- | Colors from color list are selected according to in which of the compared trees the node is contained.
selectColors :: [Int] -> [GVA.Color] -> GVAC.ColorList
selectColors inTrees currentColorList = GVAC.toColorList (map (\i -> currentColorList !! i) inTrees)

-- | A color list is sampled from the spectrum according to how many trees are compared.
makeColorList :: Int -> [GVA.Color]
makeColorList treeNumber = cList
  where cList = map (\i -> GVAC.HSV ((fromIntegral i/fromIntegral neededColors) * 0.708) 0.5 1.0) [0..neededColors]
        neededColors = treeNumber - 1

-- | Write tree representation either as dot or json to provided file path
writeTree :: String -> String -> Bool -> Gr SimpleTaxon Double -> IO ()
writeTree requestedFormat outputDirectoryPath withRank inputGraph = do
  case requestedFormat of
    "dot" -> writeDotTree outputDirectoryPath withRank inputGraph
    "json"-> writeJsonTree outputDirectoryPath inputGraph
    _ -> writeDotTree outputDirectoryPath withRank inputGraph

-- | Write tree representation as dot to provided file path.
-- Graphviz tools like dot can be applied to the written .dot file to generate e.g. svg-format images.
writeDotTree :: String -> Bool -> Gr SimpleTaxon Double -> IO ()
writeDotTree outputDirectoryPath withRank inputGraph = do
  let diagram = drawTaxonomy withRank (grev inputGraph)
  writeFile (outputDirectoryPath ++ "taxonomy.dot") diagram

-- | Write tree representation as json to provided file path.
-- You can visualize the result for example with 3Djs.
writeJsonTree :: String -> Gr SimpleTaxon Double -> IO ()
writeJsonTree outputDirectoryPath inputGraph = do
  let jsonOutput = AE.encode (grev inputGraph)
  L.writeFile (outputDirectoryPath ++ "taxonomy.json") jsonOutput
