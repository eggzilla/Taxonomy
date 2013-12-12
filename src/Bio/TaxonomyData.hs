-- | This module contains a hierarchical data structure for
--   taxonomy data


module Bio.TaxonomyData where
    
-- | RNAplex output consists of a set of interactions 
data Taxonomy = Taxonomy
  {
    taxData :: String
  }
  deriving (Show, Eq)

