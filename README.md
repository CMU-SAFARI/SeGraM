# SeGraM: A Universal Genomic Mapping Accelerator for both Sequence-to-Graph Mapping and Sequence-to-Sequence Mapping

SeGraM is a universal genomic mapping accelerator that supports both sequence-to-graph mapping and sequence-to sequence mapping, for both short and long reads. SeGraM consists of two main components: (1) MinSeed, the first minimizer-based seeding accelerator, which finds the candidate mapping locations (i.e., subgraphs) in a given genome graph; and (2) BitAlign, the first bitvector-based sequence-to-graph alignment accelerator, which performs alignment between a given read and the subgraph identified by MinSeed. MinSeed is built upon a memory-efficient minimizer-based seeding algorithm, and BitAlign is built upon our novel bitvector-based, highly-parallel sequence-to-graph alignment algorithm. 

In MinSeed, the minimizer-based seeding approach decreases the memory footprint of the index and provides speedup during seed queries. MinSeed logic requires only basic operations (e.g., comparisons, simple arithmetic operations, scratchpad read-write operations) that are implemented with simple logic. Due to frequent memory accesses required for fetching the seeds, we couple MinSeed with High-Bandwidth Memory (HBM) to enable low-latency and highly-parallel memory access, which alleviates the memory latency bottleneck. 

In BitAlign, we design a new bitvector-based alignment approach, which is amenable to efficient hardware acceleration. BitAlign employs a systolic-array-based design to circulate the internal data (i.e., bitvectors) generated by different processing elements, which provides scalability and reduces both memory bandwidth and memory footprint. In order to handle hops (i.e., non-neighbor nodes in the graph-based reference), BitAlign provides a simple design that contains queue structures between each processing element, which store the most recently generated bitvectors.

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

## Citation
>Damla Senol Cali, Konstantinos Kanellopoulos, Joel Lindegger, Zulal Bingol, Gurpreet S. Kalsi, Ziyi Zuo, Can Firtina, Meryem Banu Cavlak, Jeremie S. Kim, Nika Mansouri Ghiasi, Gagandeep Singh, Juan Gomez-Luna, Nour Almadhoun Alserr, Mohammed Alser, Sreenivas Subramoney, Can Alkan, Saugata Ghose, and Onur Mutlu.
[**"SeGraM: A Universal Hardware Accelerator for GenomicSequence-to-Graph and Sequence-to-Sequence Mapping."**](https://people.inf.ethz.ch/omutlu/pub/SeGraM_genomic-sequence-mapping-universal-accelerator_isca22.pdf)
In _Proceedings of the 49th International Symposium on Computer Architecture (ISCA),_ New York City, NY, USA, June 2022.

## Contact
Damla Senol Cali (damlasenolcali@gmail.com)
