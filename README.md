Taxonomy [![Hackage](https://img.shields.io/hackage/v/Taxonomy.svg)](https://hackage.haskell.org/package/Taxonomy) [![Build Status](https://travis-ci.org/eggzilla/Taxonomy.svg?branch=master)](https://travis-ci.org/eggzilla/Taxonomy)
=============

Haskell cabal Taxonomy libary contains tools, parsers, datastructures and visualisation
for the NCBI (National Center for Biotechnology Information) Taxonomy datasources.

It can utilize information from the Entrez REST interface (via [EntrezHTTP](https://github.com/eggzilla/EntrezHTTP),
as well as from the files of the Taxonomy database dump.

Input data is parsed into a FGL based datastructure, which enables a wealth of processing
steps like node distances, retrieval of parent nodes or extraction of
subtrees.

Trees can be visualised via dot-format (graphviz) or
via json-format (d3js).
