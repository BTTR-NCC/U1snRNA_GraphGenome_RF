# U1snRNA_GraphGenome_RF

This page shares the script we used for our paper entitled "Diversity of U1 small nuclear RNAs and Diagnostic Methods for their Mutations", published in Cancer Science, 2025. Detailed methods can be found in our manuscript and supplementary methods.

1. Mutation calling with human pangenome reference (Related to Figure 4)

   Human pangenome references can be downloaded from the Human Pangenome Reference Consortium (https://humanpangenome.org/) and the mutation calling is based on vg algorithms (https://github.com/vgteam/vg).

   Step1. Extract sequencing reads aligned on snRNA genes and variant genes from the BAM file as in the analysis with the linear genome, then align to the reference graph

      Script: MutationCall_Pangenome/bam2fq_vg_giraff_mapping.sh (bamtofastq subscript locates under MutationCall_Pangenome/subscript/)
      Sample region bed file: sample_files/gencode.v38.chr_patch_hapl_scaff.annotation.snRNA.transcript.2k.bed

   Step2. Call mutations using EB Call and bcftools.

      Script: MutationCall_Pangenome/call_bcftools.sh , MutationCall_Pangenome/call_Genomon.sh

3. Splicing-based machine learning (Random Forest Classifier; Related to Figure 5)

   Step1. Intron-centric alternative splicing analysis
     We implemented LeafCutter to perform this analysis.
     Intron clustering was run with the minimum required reads = 50 and max_intron = 500,000, and differentially spliced introns were annotated using LeafViz with GENCODE v.38 .gtf file. According to the annotation and ES calculation described above, we defined three differentially spliced introns: 1) all-DSIs, 2) Cryptic-DSIs, and 3) C6-DSIs. The cryptically spliced introns were defined in LeafCutter by intersecting junctions with GENCODE v.38 transcripts.
     For detailed usage of LeafCutter, refer to the repository (https://github.com/davidaknowles/leafcutter) and the paper (Li YI, Knowles DA, Humphrey J, et al. Annotation-free quantification of RNA splicing
using LeafCutter. Nat Genet. 2018; 50: 151-158.).

   Step2. Cross-validation and define the best differentially spliced introns (DSIs) for the estimation 

      Script: Splicing-based_ML/U1analysis.random.forest_modeleval_list.py
   
   Step3. Estimate U1 snRNA mutation status with defined DSIs

      Script: Splicing-based_ML/U1analysis.random.forest.py
      Example input:
         psitrainF = Splicing-based_ML/sample_files/Example_RFinput.tsv
         psitestF = Splicing-based_ML/Example_psitable_corrected.txt
 
      
