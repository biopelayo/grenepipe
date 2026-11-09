# lncRNA population analysis: summary

## lncRNA annotation

- Final lncRNA transcripts: **1500** (1500 novel from assembly, 0 from external annotations)
- By type: lincRNA: 858, antisense: 485, intronic: 157
- TE-like (>= 50% exonic overlap with TAIR10 TEs): 477

## Expression matrices

| set    |   features_total |   features_expressed |
|:-------|-----------------:|---------------------:|
| lncRNA |             1500 |                 1500 |
| coding |             4000 |                 4000 |

## Consensus clustering (lncRNA expression)

| set    |   k |       PAC |   mean_silhouette | cluster_sizes          | selected   |
|:-------|----:|----------:|------------------:|:-----------------------|:-----------|
| lncRNA |   2 | 0.589539  |         0.469359  | 19;29                  | False      |
| lncRNA |   3 | 0         |         0.675094  | 19;15;14               | True       |
| lncRNA |   4 | 0.101064  |         0.437465  | 18;15;1;14             | False      |
| lncRNA |   5 | 0.105496  |         0.437855  | 1;15;17;1;14           | False      |
| lncRNA |   6 | 0.100177  |         0.433457  | 2;15;15;1;14;1         | False      |
| lncRNA |   7 | 0.0930851 |         0.433741  | 2;15;14;1;14;1;1       | False      |
| lncRNA |   8 | 0.0718085 |         0.437315  | 1;15;1;14;14;1;1;1     | False      |
| lncRNA |   9 | 0.0602837 |         0.433464  | 1;15;1;13;1;14;1;1;1   | False      |
| lncRNA |  10 | 0.0434397 |         0.435696  | 1;15;1;12;14;1;1;1;1;1 | False      |
| coding |   2 | 0.160461  |         0.185723  | 19;29                  | False      |
| coding |   3 | 0         |         0.237424  | 19;15;14               | True       |
| coding |   4 | 0.190603  |         0.150041  | 11;15;8;14             | False      |
| coding |   5 | 0.190603  |         0.0771823 | 14;5;10;5;14           | False      |
| coding |   6 | 0.171099  |         0.0732389 | 13;5;10;4;14;2         | False      |
| coding |   7 | 0.170213  |         0.0724804 | 12;5;4;4;14;6;3        | False      |
| coding |   8 | 0.177305  |         0.0723288 | 5;5;7;4;4;14;6;3       | False      |
| coding |   9 | 0.166667  |         0.0736409 | 4;4;7;4;4;14;7;3;1     | False      |
| coding |  10 | 0.161348  |         0.0741475 | 5;4;7;4;4;14;5;2;2;1   | False      |

### lncRNA vs coding clustering agreement

| comparison               |   k_lncrna |   k_coding |   ARI |
|:-------------------------|-----------:|-----------:|------:|
| lncRNA_vs_coding         |          3 |          3 |     1 |
| lncRNA_vs_external_label |          3 |        nan |     1 |
| coding_vs_external_label |        nan |          3 |     1 |

## lncRNA signature (elastic-net, nested CV)

- Signature size: **49** lncRNAs
|   fold |   n_test |   accuracy |   n_features |
|-------:|---------:|-----------:|-------------:|
|      1 |       10 |          1 |         56   |
|      2 |       10 |          1 |         30   |
|      3 |       10 |          1 |         58   |
|      4 |        9 |          1 |         48   |
|      5 |        9 |          1 |         52   |
|    nan |       48 |          1 |         48.8 |

Top signature members:

| feature   |   max_abs_coef |   classes |   coefs | lncrna_type   | te_like   |
|:----------|---------------:|----------:|--------:|:--------------|:----------|
| LNC00663  |     0.0160784  |         1 |   0.016 | lincRNA       | False     |
| LNC00877  |     0.0154494  |         1 |  -0.015 | lincRNA       | True      |
| LNC00521  |     0.0146284  |         1 |  -0.015 | antisense     | False     |
| LNC00048  |     0.0140447  |         1 |   0.014 | lincRNA       | False     |
| LNC00404  |     0.0138024  |         1 |   0.014 | lincRNA       | False     |
| LNC01122  |     0.0133689  |         3 |   0.013 | intronic      | False     |
| LNC00656  |     0.0131254  |         3 |  -0.013 | lincRNA       | True      |
| LNC00692  |     0.0130923  |         2 |  -0.013 | antisense     | False     |
| LNC00206  |     0.0127341  |         2 |  -0.013 | lincRNA       | True      |
| LNC00033  |     0.0122474  |         2 |  -0.012 | lincRNA       | False     |
| LNC01235  |     0.0122299  |         1 |   0.012 | antisense     | False     |
| LNC01292  |     0.0116122  |         2 |  -0.012 | lincRNA       | True      |
| LNC00117  |     0.0114845  |         1 |   0.011 | intronic      | True      |
| LNC00157  |     0.0114173  |         3 |   0.011 | lincRNA       | True      |
| LNC00963  |     0.0112953  |         2 |  -0.011 | lincRNA       | False     |
| LNC00068  |     0.0110834  |         3 |   0.011 | lincRNA       | False     |
| LNC00893  |     0.010674   |         3 |   0.011 | lincRNA       | False     |
| LNC00574  |     0.0102951  |         1 |   0.01  | intronic      | False     |
| LNC01457  |     0.00947402 |         3 |   0.009 | antisense     | True      |
| LNC00081  |     0.00920525 |         2 |   0.009 | lincRNA       | False     |
| LNC00599  |     0.00907193 |         2 |   0.009 | lincRNA       | False     |
| LNC01434  |     0.00889431 |         2 |  -0.009 | lincRNA       | False     |
| LNC00517  |     0.00817176 |         1 |  -0.008 | lincRNA       | False     |
| LNC00241  |     0.0081263  |         3 |   0.008 | lincRNA       | False     |
| LNC01301  |     0.00726621 |         1 |   0.007 | lincRNA       | False     |

## Genotype x condition (DESeq2 LRT)

- Features with significant G x E (padj < 0.05): **203** of which lncRNAs: 77

## cis-eQTL (MatrixEQTL)

- lncRNA eGenes at FDR 0.05: **63**
