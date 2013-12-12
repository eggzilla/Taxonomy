-- | Parse and process taxonomy data

module Bio.Taxonomy (                      
                       module Bio.TaxonomyData
                      ) where

import Bio.TaxonomyData
import Text.ParserCombinators.Parsec
import Text.ParserCombinators.Parsec.Token
import Text.ParserCombinators.Parsec.Language (emptyDef)    
import Control.Monad

