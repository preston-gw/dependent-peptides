# News
News relating to the research article 'Visualisation tools for dependent peptide searches to support the exploration of *in vitro* protein modifications' by Preston and co-authors (https://doi.org/10.1371/journal.pone.0235263)
## Update regarding dependent-peptide matching
*Posted to the [article's comments section](https://journals.plos.org/plosone/article/comments?id=10.1371/journal.pone.0235263) on 14 Feb 2021*

BACKGROUND. Scripts II, III and IV perform pairwise comparisons of dependent peptides (DPs; handled as combinations of sequence and mass shift). If two DPs share the same sequence, and their mass shifts are close enough, they are considered to match.

PROBLEM. Script III is inappropriately matching DPs with the annotations ‘Comp.:13C-12C’ and ‘Deamidation’. Other problematic pairs are ‘di-Oxidation/Sulfide’ and ‘Comp.:S-NH3O/Comp.:O-NH3’. Further pairs with very close shifts could also be problematic. Script II matches DPs using the same basic method as Script III, and will be similarly affected. Script IV uses a different method, and may not be affected to the same extent.

SOLUTIONS. For Script III, I have identified three potential solutions: (i) excluding interfering modifications (e.g., ‘Comp.:13C-12C’); (ii) requiring annotations to match; and (iii) halving the tolerance. I tested these solutions by making variants of Script III and reprocessing the BSA data reported in the paper. A match was considered inappropriate if the DPs’ annotations were different. The matching of unannotated DPs with different mass shifts was not tested for presently (however, see CONCLUSIONS).

RESULTS. While none of the solutions is perfect, each was effective in reducing the frequency with which differentially annotated DPs were matched (from 32 observations to 0-6 observations). For each solution, only a modest (<8%) change in the overall number of DPs post enrichment was observed. There was no change in the number of marker DPs (mass shift = 125 or 143 Da) observed pre or post enrichment. There was only a modest decrease in the calculated enrichment factor for the marker DPs (from 6.0 to 5.7-5.8).

CONCLUSIONS. Users could implement a solution that is appropriate to their experiment. A good general solution seems to be reducing the tolerance, because this does not require affected DPs to have been annotated.
