#!/bin/bash

detectedits="/PHShome/mi825/Desktop/CRISPR_WGS_Detection/edits_detection_comparison/src/revision/detect_edits.py"
genome="/data/pinello/COMMON_DATA/REFERENCE_GENOMES/Broad/hg38/Homo_sapiens_assembly38.fasta"
dnmt1site3_gm12878="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/DNMT1Site3.bam"
dnmt1site3_k562="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/K562_DNMT1Site3.cram"
outdir="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/detectWithOtherTools/manuel_experiments/VCFs_revision_casoffinder"

# EMX1 - GM12878 - Mutect2
emx1_targets="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/offtargetDetection/casoffinder/EMX1.GAGTCCGAGCAGAAGAAGAA.gap1.offby4.NNN.CasOFFinder.out"
emx1_bam="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/EMX1.cram"
echo "python $detectedits --tool mutect2 --targets $emx1_targets --genome $genome --bam1 $dnmt1site3_gm12878 --bam2 $emx1_bam --normal-sample DNMT1Site3 --out "${outdir}/EMX1_GM12878_MUTECT2" --casoffinder --threads 16"

# # EMX1 - K562 - Mutect2
# emx1_bam="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/K562_EMX1.cram"
# python $detectedits --tool mutect2 --targets $emx1_targets --genome $genome --bam1 $dnmt1site3_k562 --bam2 $emx1_bam --normal-sample DNMT1Site3 --out "${outdir}/EMX1_K562_MUTECT2" --casoffinder --threads 16

# # HEK4 - GM12878 - Mutect2
# hek4_targets="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/offtargetDetection/casoffinder/HEK4.GGCACTGCGGCTGGAGGTGG.gap1.offby4.NNN.CasOFFinder.out"
# hek4_bam="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/HEKSite4.bam"
# python $detectedits --tool mutect2 --targets $hek4_targets --genome $genome --bam1 $dnmt1site3_gm12878 --bam2 $hek4_bam --normal-sample DNMT1Site3 --out "${outdir}/HEK4_GM12878_MUTECT2" --casoffinder --threads 16

# # HEK4 - K562 - Mutect2
# hek4_bam="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/K562_HEKSite4.cram"
# python $detectedits --tool mutect2 --targets $hek4_targets --genome $genome --bam1 $dnmt1site3_k562 --bam2 $hek4_bam --normal-sample DNMT1Site3 --out "${outdir}/HEK4_K562_MUTECT2" --casoffinder --threads 16

# # RNF2 - GM12878 - Mutect2
# rnf2_targets="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/offtargetDetection/casoffinder/RNF2.GTCATCTTAGTCATTACCTG.gap1.offby4.NNN.CasOFFinder.out"
# rnf2_bam="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/RNF2.bam"
# python $detectedits --tool mutect2 --targets $rnf2_targets --genome $genome --bam1 $dnmt1site3_gm12878 --bam2 $rnf2_bam --normal-sample DNMT1Site3 --out "${outdir}/RNF2_GM12878_MUTECT2" --casoffinder --threads 16

# # RNF2 - K562 - Mutect2
# rnf2_bam="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/K562_RNF2.cram"
# python $detectedits --tool mutect2 --targets $rnf2_targets --genome $genome --bam1 $dnmt1site3_k562 --bam2 $rnf2_bam --normal-sample DNMT1Site3 --out "${outdir}/RNF2_K562_MUTECT2" --casoffinder --threads 16

# # VEGFA3 - GM12878 - Mutect2
# vegfa3_targets="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/offtargetDetection/casoffinder/VEGFA3.GGTGAGTGAGTGTGTGCGTG.gap1.offby4.NNN.CasOFFinder.out"
# vegfa3_bam="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/VEGFASite3.bam"
# python $detectedits --tool mutect2 --targets $vegfa3_targets --genome $genome --bam1 $dnmt1site3_gm12878 --bam2 $vegfa3_bam --normal-sample DNMT1Site3 --out "${outdir}/VEGFA3_GM12878_MUTECT2" --casoffinder --threads 16

# # VEGFA3 - K562 - Mutect2
# vegfa3_bam="/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/wgs/GM12878-Cas9/WGS1000/data/K562_VEGFASite3.cram"
# python $detectedits --tool mutect2 --targets $vegfa3_targets --genome $genome --bam1 $dnmt1site3_k562 --bam2 $vegfa3_bam --normal-sample DNMT1Site3 --out "${outdir}/VEGFA3_K562_MUTECT2" --casoffinder --threads 16
