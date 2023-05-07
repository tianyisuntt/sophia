# sophia

DESCRIPTION
--------------------------------

A _tower_ is a sequence of simplicial complexes connected by simplicial maps. A simplicial map can be decomposed in a composition of two different elementary operations: the _inclusion_ of a simplex and the _contraction_ of two vertices.
Sophia computes the persistent diagram of a tower by converting it first into an equivalent filtration, which is a sequence of simplicial complexes connected by inclusions. From there, it uses classical algorithms (see [1] and [2]) to compute the persistence pairs.

The convertion is a variation of the coning strategy by Dey et al. [3], yielding a filtration that is asymptotically only marginally larger than the tower. For more details about the construction, we refer to [4]. It is implemented as a streaming algorithm.  
To be able to compute the persistence of really long towers, when they have relatively small complexes, the persistence algorithm was also implemented as a streaming algorithm. The previous converting algorithm, which can be used simultaneously, signals simplices as "inactive", which allows us to clean out the boundary matrix. This way, the space complexity of the algorithm does not depend on the length of the tower, but the maximal size of any subcomplex within the tower.

Sophia is published under the GNU Lesser General Public License.  
**Author:** Hannah Schreiber (hannah.schreiber.k [at] gmail.com)  
**Contributor:** Michael Kerber

#### References
[1] C. Chen and M. Kerber. Persistent homology computation with a twist. In _European Workshop on Computational Geometry (EuroCG)_, pages 197-200, 2011.  
[2] U. Bauer, M. Kerber, and J. Reininghaus. Clear and compress: Computing persistent homology in chunks. In _Topological Methods in Data Analysis and Visualization III_, Mathematics and Visualization, pages 103-117. Springer, 2014.  
[3] T. Dey, F. Fan, and Y. Wang. Computing Topological Persistence for Simplicial Maps. In _ACM Symposium on Computational Geometry (SoCG)_, pages 345-354, 2014.  
[4] M. Kerber, H. Schreiber. Barcodes of Towers and a Streaming Algorithm for Persistent Homology. <https://arxiv.org/abs/1701.02208>, TU Graz, 2016

INSTALLATION
--------------------------------

To compile Sophia easily, we recommand to use CMake (<https://cmake.org/>).

Sophia has two options related to two external libraries: PHAT (<https://bitbucket.org/phat-code/phat>) and Gudhi (<http://gudhi.gforge.inria.fr/>), see USAGE section for more details. They also require to have installed the libraries OpenMP and Boost (version 1.48.0 or more recent) respectively. To activate these options, please specify the path to the include directory of PHAT and/or Gudhi in the "CMakeLists.txt" file at line 8 and 9:  
> \#set(PHAT /path-to-PHAT-include-directory/)  
> \#set(GUDHI /path-to-Gudhi-include-directory/)  

and decomment them by removing the '#' symbol.

To compile, run in a terminal:
> $	cd /path-to-software/  
> $	mkdir build  
> $	cd build  
> $	cmake ..  
> $	make

And then to use the software:
> $	./sophia [args]

USAGE
--------------------------------

For Information about the different file formats mentionned in this section, refer to the section FILE FORMATS.  
To each of the following options can be added the argument "--silent" to avoid output in the terminal.

By default, there are 4 options:

1.	`./sophia -cv input-file-name output-file-name`  
Builds the equivalent filtration of the tower, whose description is in _input-file-name_. Stores the filtration in _output-file-name_ with a VERTICES Format.

2.	`./sophia -cf input-file-name output-file-name`  
Builds the equivalent filtration of the tower, whose description is in _input-file-name_. Stores the filtration in _output-file-name_ with a FACES Format.

3.	`./sophia -cp input-file-name output-file-name chunk-size`  
Builds the equivalent filtration of the tower, whose description is in _input-file-name_ and simultaneously computes its persistence diagram, using a chunk of size _chunk-size_. Stores the persistence pairs in _output-file-name_.

4.	`./sophia -p input-file-name output-file-name chunk-size`  
Computes the persistence diagram of the filtration, whose description is in _input-file-name_, using a chunk of size _chunk-size_. The format of _input-file-name_ needs to be the FACES Format. Stores the persistence pairs in _output-file-name_.

  Two additional options can be activated (see section INSTALLATION). Note that they are thought to be more convenient and not efficient; they are often slower than using PHAT or Gudhi separately because the filtration is kept in memory:

5.	`./sophia -cphat input-file-name output-file-name`  
Builds the equivalent filtration of the tower, whose description is in _input-file-name_. Then uses the PHAT library, with its default parameters, to compute its persistent diagram. Stores the persistence pairs in _output-file-name_, in the same output format PHAT uses, but without the first line.

6.	`./sophia -cgudhi input-file-name output-file-name`  
Builds the equivalent filtration of the tower, whose description is in _input-file-name_. Then uses the Gudhi library, with its default parameters, to compute its persistent diagram. Stores the persistence pairs in _output-file-name_, in the same output format Gudhi uses.

FILE FORMATS
--------------------------------

For each file format, every simplex has an unique integer _identifier_. The identifiers are assigned to the simplices in increasing order by appearence, starting with 0.

In this section, "f^i" will always correspond to an identifier of a simplex and "v^i" more specifically to the identifier of a vertex.

The example files are in the "example" directory.

#### Tower Format
Options 1 to 3 and 5 to 6 use the same format to describe the input Tower: one line corresponds to one action.
In a case of an inclusion of a _d_-simplex _s_:
> [ts] i v^1 ... v^{d+1}

where _v^i_ in _s_ and for each _i_ < _j_ in {1,...,d+1}, _v^i_ < _v^j_,
and _ts_ is an optional time indicator. If there is no time indication, time will start at 0 and increase by one at each insertion or contraction.
In a case of a contraction:
> [ts] c v^d v^k

where _v^d_ and _v^k_ are the vertices to be contracted. From here on, the remaining vertex needs to be refered by _v^k_ and **NOT** _v^d_,
and _ts_ is the optional time indicator.  
It is possible to insert comments in the file by letting the comment line begin with '#'.  
An example is file "example\_tower.txt".

#### VERTICES Format
Each line corresponds to an inclusion of a _d_-simplex _s_:
> d v^1 ... v^{d+1} ts

where _v^i_ in _s_ and for each _i_ < _j_ in {1,...,d+1}, _v^i_ < _v^j_,
and _ts_ corresponds to the time _s_ was included.  
It corresponds to the input format of Gudhi.  
An example is file "example\_vertices\_filtration.txt".

#### FACES Format
Each line corresponds to either an inclusion of a _d_-simplex _s_:
> ts i f^1 ... f^{d+1} f^s

where _f^i_ is a face of _s_ and for each _i_ < _j_ in {1,...,d+1}, _f^i_ < _f^j_,
_f^s_ is the identifier of _s_,
and _ts_ corresponds to the time s was included.  
Or signals the inactivity of a simplex _s_:
> a f^s

where _f^s_ is the identifier of _s_.  
An example is file "example\_faces\_filtration.txt".

#### Persistence Pairs Format
Each line corresponds to one pair:
> dim b d

where _b_ is the birth time, _d_ the death time and _dim_ the dimension of the homology class.  
An example is file "example\_res.txt".
