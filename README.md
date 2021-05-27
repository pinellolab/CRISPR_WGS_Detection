# CRISPR WGS Detection

Whole-genome-sequencing (WGS) porvides the ability to investigate genome edits at every position in the genome. Unlike other methods that enrich for sites of potental genome editing using biochemical properties, WGS reads the DNA sequence at all sites in the genome. It is noted that WGS may be a very costly way to identify on- and off-target genome editing, but it is not subject to the same assumptions that are present in other methods.

We are primarily concerned with two use cases -- one in which the CRISPR sgRNA used to perform the editing is known (*Supervised Method*), and one in which the sgRNA is unknown (*Unsupervised Method*). 

## Supervised Method

If the sgRNA used to perform genome editing is known, we can use WGS to identify editing at on- and off- targets in an unbiased manner. Briefly, at target regions of interest (e.g. suggested by in-silico off-target prediction) we extract reads from the WGS sequencing and analyze reads. Our overall goal is to distinguish true editing from other sources of noise. 

We first use Casoffinder (Bae et al., Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases, Bioinformatics 2014) to identify all putative cut sites with up to 4 mismatches and 1 bulge with respect to the guide sequence. Next, at each putative site, we define the *cutting rate* and the *correlation* between the predicted and observed alleles. 

The *cutting rate* is simply the number of reads with indels divided by the total number of reads. High cutting rate indicates that the site contains frequent indels. 

The *correlation* is the correlation between the observed frequency of each allele containing an allele and the predicted frequency of that allele. High correlation means that the observed frequencies of indel alleles is very similar to the predicted frequencies of those alleles. We predicted editing outcomes using InDelphi (Shen et al., Predictable and precise template-free CRISPR editing of pathogenic variants, Nature 2018). 

At sites of real editing, we expect to see high cutting rate and high correlation in the treatment sample relative to the control sample. At noisy sites with no editing, even if there is a high cutting rate, we expect that the frequency of each allele does not correlate with the predicted frequency. We then define two scores, delta correlation and delta cutting rate which are simply the difference between those values in the treated and control samples, that is,

1)	delta_correlation = correlation(Treated_indel_reads, inDelphi_score) - correlation(Control_indel_reads, inDelphi_score) 

2)	delta_cutting_rate = Treated_cutting_rate – Control_cutting_rate

<img src=https://github.com/pinellolab/CRISPR_WGS_Detection/blob/main/supervised/wgs.png alt="Supervised Detection" width="70%" height="70%" align="bottom" />
<img src=https://github.com/pinellolab/CRISPR_WGS_Detection/blob/main/supervised/supervised_approach.png alt="Supervised Detection" width="65%" height="65%" align="bottom" />


## Unsupervised Method

For this task, we have leveraged characteristics of genome editing to detect and characterize genome editing events. In the case for which the on-target is known, we can computationally enumerate sites with sequence similarity to the on-target where we expect that off-target editing may occur. Using our whole-genome sequencing data, we can examine these computationally-predicted sites for editing. In addition, genome edits produced by double-strand cleavage can also be predicted based on sequence homology at the cut site and other characteristics. We can compare the predicted editing outcomes with the observed editing outcomes in our WGS data to determine which sites may evidence true nuclease editing. 

Our current approach to detect genome editing genome-wide depends on the calculation of the Indel Ratio, computed as follows:
1.	The entire genome is scanned for insertions or deletions (Indels), and the number of indels in overlapping 40bp windows is tracked.
2.	The ratio of unique indels in the edited WGS sample to the control WGS sample is calculated, weighted by the number of reads in each sample. This is the Indel Ratio.
3.	40bp windows are sorted by the Indel Ratio, and windows with the highest Indel Ratio indicate sites of possible editing.
