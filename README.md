# SeGraM: A Universal Genomic Mapping Accelerator for both Sequence-to-Graph Mapping and Sequence-to-Sequence Mapping

SeGraM is a universal genomic mapping accelerator that supports both sequence-to-graph mapping and sequence-to sequence mapping, for both short and long reads. SeGraM consists of two main components: (1) MinSeed, the first minimizer-based seeding accelerator, which finds the candidate mapping locations (i.e., subgraphs) in a given genome graph; and (2) BitAlign, the first bitvector-based sequence-to-graph alignment accelerator, which performs alignment between a given read and the subgraph identified by MinSeed. MinSeed is built upon a memory-efficient minimizer-based seeding algorithm, and BitAlign is built upon our novel bitvector-based, highly-parallel sequence-to-graph alignment algorithm. 

In MinSeed, the minimizer-based seeding approach decreases the memory footprint of the index and provides speedup during seed queries. MinSeed logic requires only basic operations (e.g., comparisons, simple arithmetic operations, scratchpad read-write operations) that are implemented with simple logic. Due to frequent memory accesses required for fetching the seeds, we couple MinSeed with High-Bandwidth Memory (HBM) to enable low-latency and highly-parallel memory access, which alleviates the memory latency bottleneck. 

In BitAlign, we design a new bitvector-based alignment approach, which is amenable to efficient hardware acceleration. BitAlign employs a systolic-array-based design to circulate the internal data (i.e., bitvectors) generated by different processing elements, which provides scalability and reduces both memory bandwidth and memory footprint. In order to handle hops (i.e., non-neighbor nodes in the graph-based reference), BitAlign provides a simple design that contains queue structures between each processing element, which store the most recently generated bitvectors.

## BitAlign

After MinSeed or any seeding tool determines the subgraphs to perform alignment for each query read, for each (read, subgraph) pair, BitAlign calculates the edit distance and the corresponding alignment between the two. In order to provide an efficient, hardware-friendly, and low-cost solution, we modify the sequence alignment algorithm of [GenASM](https://github.com/CMU-SAFARI/GenASM/), which is bitvector-based, to support sequence-to-graph alignment, and we exploit the bit-parallelism that the GenASM algorithm provides.

#### GenASM
GenASM makes the bitvector-based Bitap algorithm suitable for efficient hardware implementation. GenASM shares a common characteristic with the well-known DP-based algorithms: both algorithms operate on tables. The key difference between GenASM-based alignment and DP-based alignment is that cell values are bitvectors in GenASM, whereas cell values are numerical values in DP-based algorithms. In GenASM, the rules for computing new cell values can be formulated as simple bitwise operations, which are particularly easy and cheap to implement in hardware. Unfortunately, GenASM is limited to sequence-to-sequence align- ment. We build on GenASM’s bitvector-based algorithm to develop our new sequence-to-graph alignment algorithm, BitAlign.

#### BitAlign
There is a major difference between sequence-to-sequence alignment and sequence-to-graph alignment: for the current character, sequence-to-sequence alignment needs to know about only the neighboring (i.e., previous/adjacent) text character, whereas sequence-to-graph alignment must incorporate non-neighboring characters as well whenever there is an edge (i.e., hop) from the current character to the non-neighboring character. To ensure that each of these data depen- dencies can be resolved as quickly as possible, we topologically sort the input graph during pre-processing.

BitAlign starts with a linearized and topologically sorted input subgraph. This ensures that (1) we can iterate over each node of the input graph sequentially, and (2) all of the required bitvectors for the current iteration have already been generated in previous iterations. Besides the subgraph, BitAlign also takes the query read and the edit distance threshold (i.e., maximum number of edits to tolerate when performing approximate string matching) as inputs.

Similar to GenASM, as a pre-processing step, we generate four pattern bitmasks for the query read (one for each character in the alphabet: A, C, G, T). Unlike in GenASM, which stores only the most recent status bitvectors (i.e., R[d] bitvectors, where 0 ≤ 𝑑 ≤ 𝑘) that hold the partial match information, BitAlign needs to store all of the status bitvectors for all of the text iterations (i.e., allR[n][d], where 𝑛 is the length of the linearized reference subgraph). These allR[n][d] bitvectors will be later used by the traceback step of BitAlign (BitAlign-TB).

Next, BitAlign iterates over each node of the linearized graph and retrieves the pattern bitmask for each node, based on the character stored in the current node. Unlike in GenASM, when computing three of the intermediate bitvectors (i.e., match, substitution, and deletion), BitAlign incorporates the hops as well by examining all successor nodes that the current node has. When calculating the deletion (𝐷), substitution (𝑆), and match (𝑀) bitvectors, we take the hops into consideration, whereas when calculating the insertion (𝐼) bitvector, we do not need to, since an insertion does not consume a character from the reference subgraph, but does so from the query read only. After computing all of these intermediate bitvectors, we store only the R[d] bitvector (i.e., ANDed version of the intermediate bitvectors) in memory. After completing all iterations, we perform traceback by traversing the stored bitvectors in the opposite direction to find the optimal alignment (based on the user-supplied alignment scoring function).

#### Implementation Details

- To perform GenASM-like traceback, BitAlign stores 3(𝑘 + 1) bitvectors per graph edge (similar to how GenASM stores three out of the four intermediate vectors), where 𝑘 is the edit distance threshold. Since the number of edges in the graph can only be bounded very loosely, the potential memory footprint increases significantly, which is expensive to implement both in software and hardware. We solve this problem by storing only 𝑘 + 1 bitvectors per node (i.e., R[d] bitvectors), from which the 3(𝑘 + 1) bitvectors per edge can be regenerated on-demand during traceback. While this modification incurs small additional computational overhead, it decreases the memory footprint of the algorithm by at least 3×. Since the main area and power cost of the alignment hardware comes from memory, we find this trade-off favorable.

- Similar to GenASM, BitAlign also follows the divide-and-conquer approach, where we divide the linearized subgraph and the query read into overlapping windows and execute BitAlign for each window. After all windows’ traceback outputs are found, we merge them to find the final traceback output. This approach helps us to decrease the memory footprint and lower the complexity of the algorithm.

- We implement the hop information between nodes of the graph as an adjacency matrix called HopBits. Based on the HopBits entry of the current text character, either the actual hop bitvector (if the HopBits entry is 1), or a bitvector containing all ones such that it will not have any effect on the bitwise operations (if the HopBits entry is 0), is used when calculating the match, deletion, and substitution bitvectors.

- In order to decrease the size of the HopBits matrix, we also provide the option to limit the distance between the farthest node to take into account and the current node (i.e., hop limit).

## Running SeGraM

Call the following two functions in `src/graph.c` and `src/align.c` files in your C code, respectively, or update the existing main() function in `src/main.c` file:

```bash
struct SeqNode* generateGraphFromGFA(char *filename, int *numNodes, int *numEdges)
bitalign_aligner(struct SeqNode *nodes, char *pattern, int startNode, int startOffset, int endNode, int endOffset, int k, int hopLimit, int W, int O, int scoreM, int scoreS, int scoreOpen, int scoreExtend)
```
For example:

```bash
    int numNodes, numEdges;
    
    struct SeqNode *nodes = generateGraphFromGFA("test.gfa", &numNodes, &numEdges);
    bitalign_aligner(nodes, "ACGTCATGCAGTCGTAACGTAGTCGTCACAGTCAGTCGTAGCTAGTA", 0, 0, 3, 41, 1, 5, 200, 0, 1, 1, 1, 1);
    
    deleteGraph(nodes, numNodes);
```

#### Limitations

- The current implementation follows the divide-and-conquer approach, and does not have the ability to disable it.
- The current implementation expects the start node and the start offset within the start node, along with the possible end node and the end offset within the end node to perform the alignment. In order to run BitAlign as a standalone tool without the need of an initial seeding approach, this requirement should be removed.
- The current implementation expects the graph representation of the reference text to be generated using the `struct SeqNode* generateGraphFromGFA(char *filename, int *numNodes, int *numEdges)` function in `graph.c` file under `src/`.

## Datasets

We evaluate SeGraM using the latest major release of the human genome assembly, GRCh38, as the starting reference genome. To incorporate known genetic variations and thus form a genome graph, we use 7 VCF files for HG001-007 from the GIABproject (v3.3.2). Across the 24 graphs generated (one for each chromosome; 1–22, X, Y), in total, we have 20.4M nodes, 27.9 M edges, 3.1B sequence characters, and 7.1M variations. 

For the read datasets, we generate four sets of long reads (i.e., PacBio and ONT datasets) using PBSIM2 and three sets of short reads (i.e., Illumina datasets) using Mason. For the PacBio and ONT datasets, we have reads of length 10kbp, each simulated with 5% and 10% error rates. The Illumina datasets have reads of length 100bp, 150bp, and 250bp, each simulated with a 1% error rate. Each dataset has 10,000 reads.

All our prepared datasets can be downloaded from [this link](https://drive.google.com/file/d/18eSrcC1mCRCy9TUI2xvq8MH3mWykEunV/view?usp=sharing). The unzipped directory has the following structure:

```
└── datasets
    └── graphs
        └── gfa_files : our graph files (1 for each chromosome: 1-22, X, and Y) in GFA format
        └── vg_files : our graph files (1 for each chromosome: 1-22, X, and Y) in VG format
        └── gfa_files : our graph files (1 for each chromosome: 1-22, X, and Y) in FASTA format (i.e., each node of the entry has its own row so that we can use these files to generate the index files (.mmi) with minimap2)
    └── reads
        └── illumina_reads : our simulated short reads with 1% error rate and 100/150/250bp length
        └── pacbio_ont_reads : our simulated 10k-length long reads with 5% and 10% error rates and different error profiles for PacBio and ONT
```

## Citing SeGraM
If you use SeGraM in your work, please cite:
- >Damla Senol Cali, Konstantinos Kanellopoulos, Joel Lindegger, Zulal Bingol, Gurpreet S. Kalsi, Ziyi Zuo, Can Firtina, Meryem Banu Cavlak, Jeremie S. Kim, Nika Mansouri Ghiasi, Gagandeep Singh, Juan Gomez-Luna, Nour Almadhoun Alserr, Mohammed Alser, Sreenivas Subramoney, Can Alkan, Saugata Ghose, and Onur Mutlu.
[**"SeGraM: A Universal Hardware Accelerator for Genomic Sequence-to-Graph and Sequence-to-Sequence Mapping."**](https://people.inf.ethz.ch/omutlu/pub/SeGraM_genomic-sequence-mapping-universal-accelerator_isca22.pdf)
In _Proceedings of the 49th International Symposium on Computer Architecture (ISCA),_ New York City, NY, USA, June 2022.

- >Damla Senol Cali, Gurpreet S. Kalsi, Zülal Bingöl, Can Firtina, Lavanya Subramanian, Jeremie S. Kim, Rachata Ausavarungnirun, Mohammed Alser, Juan Gomez-Luna, Amirali Boroumand, Anant Nori, Allison Scibisz, Sreenivas Subramoney, Can Alkan, Saugata Ghose, and Onur Mutlu.
[**"GenASM: A High-Performance, Low-Power Approximate String Matching Acceleration Framework for Genome Sequence Analysis."**](https://people.inf.ethz.ch/omutlu/pub/GenASM-approximate-string-matching-framework-for-genome-analysis_micro20.pdf)
In _Proceedings of the 53rd International Symposium on Microarchitecture (MICRO),_ Virtual, October 2020.

## Contact
Damla Senol Cali (damlasenolcali@gmail.com)
