## Welcome ##

Welcome to the Lemon-Tree project. Lemon-Tree is a "one-stop shop" software suite for module network inference that aims to make module network methods available to a broader community of users while simultaneously
facilitating and encouraging collaborations between developers. Module networks are probabilistic graphical models which consist of modules of coregulated genes and their regulatory programs, and they can be inferred from integrated genomic data compendia (including gene expression, microRNA expression, copy number variation, and more) measured in a large number of individuals or experimental conditions.

## Download ##

Please note that due to a change in Google code downloads policy, all the files available for **download** can now be found on the home page, under **"External links"** section.

## Release info ##


### Version 3.0.`*` ###

For this first public release of Lemon-Tree we have ported all of the data structures and algorithms from our previous tool [LeMoNe](http://bioinformatics.psb.ugent.be/software/details/LeMoNe) and developed a number of extensions and improvements, including:

  * A complete Java implementation of all steps in the model learning process with uniform input/output format.
  * A novel algorithm based on rigorous algebraic results to reconstruct overlapping consensus clusters from an ensemble of clustering results, each of which may be overlapping or not.
  * Novel methods to incorporate _discrete_ regulators (e.g. eQTLs, CNVs, clinical grades) into the model, thereby greatly expanding the number and types of data that can be integrated.
  * Novel methods to visualise co-expression modules and their predicted regulators allowing for better interpretation and hypothesis generation.



## Where to go next? ##

  * Start with the step-by-step [tutorial](https://code.google.com/p/lemon-tree/wiki/Tutorial).
  * Read more about module networks and their applications in one of the many publications compiled in the [module network bibliography](https://code.google.com/p/lemon-tree/wiki/Bibliography).
  * Compare features of Lemon-Tree to existing module network software tools in our [comparison table](https://code.google.com/p/lemon-tree/wiki/ModuleNetworkTools).
  * Join the [discussion group](http://groups.google.com/group/lemon-tree-discuss) if the above does not answer all your questions.
  * And last but not least, if you are a developer and you are interested in porting your existing module network code to the Lemon-Tree project or you have ideas for further project extensions, or want to contribute for testing and documentation, please   [join the project](http://code.google.com/p/support/wiki/HowToJoinAProject)!

## Please cite ##

The Lemon-Tree software (version 3.0.`*`) is described in the following publication, part of the [PLOS Computational Biology Software collection](http://www.ploscollections.org/article/browse/issue/info%3Adoi%2F10.1371%2Fissue.pcol.v03.i10):

  * Bonnet E, Calzone L, Michoel T. (2015) [Integrative multi-omics module network inference with Lemon-Tree](http://dx.doi.org/10.1371/journal.pcbi.1003983). PLoS Comput Biol 11(2): e1003983.

Please cite this manuscript if you used Lemon-Tree in your research, in addition to any other relevant publications from the [module network bibliography](https://code.google.com/p/lemon-tree/wiki/Bibliography).

## Who are we? ##

Current project members/developers are:

  * [Eric Bonnet](http://eric.d.bonnet.free.fr) is a senior research engineer at the Institut Curie.
  * [Tom Michoel](http://lab.michoel.info) is a group leader at the Roslin Institute.

Parts of the code ported from [LeMoNe](http://bioinformatics.psb.ugent.be/software/details/LeMoNe) were developed by:

  * [Anagha Joshi](http://www.roslin.ed.ac.uk/anagha-joshi)
  * [Steven Maere](http://www.psb.ugent.be/esb)
