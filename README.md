# U1snRNA_GraphGenome_RF

This page shares the script we used for our paper entitled "Diversity of U1 small nuclear RNAs and Diagnostic Methods for their Mutations", published in Cancer Science, 2025. Detailed methods can be found in our manuscript and supplementary methods.

1. Mutation calling with human pangenome reference (Related to Figure 4)

2. Splicing-based machine learning (Random Forest Classifier; Related to Figure 5)

   Step1. Intron-centric alternative splicing analysis
     We implemented LeafCutter to perform this analysis.
     Intron clustering was run with the minimum required reads = 50 and max_intron = 500,000, and differentially spliced introns were annotated using LeafViz with GENCODE v.38 .gtf file. According to the annotation and ES calculation described above, we defined three differentially spliced introns: 1) all-DSIs, 2) Cryptic-DSIs, and 3) C6-DSIs. The cryptically spliced introns were defined in LeafCutter by intersecting junctions with GENCODE v.38 transcripts.
     For detailed usage of LeafCutter, refer to the repository (https://github.com/davidaknowles/leafcutter) and the paper (Li YI, Knowles DA, Humphrey J, et al. Annotation-free quantification of RNA splicing
using LeafCutter. Nat Genet. 2018; 50: 151-158.).

   Step2. Cross-validation and consideration of the best model

         Script: U1analysis.random.forest_modeleval_list.py
   
   Step3. Estimate U1 snRNA mutation status

         Script: U1analysis.random.forest.py

      
