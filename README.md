# uncertainty

Code and data to accompany the paper:

    Shravan Vasishth and Andrew Gelman. How to embrace variation and accept uncertainty in linguistic and psycholinguistic data analysis. Linguistics, 59:1311--1342, 2021
    doi: https://doi.org/10.1515/ling-2019-0051

Please contact Shravan Vasishth if there are any problems with this repository.

The files:

- The directory R has some functions that I used in this work.
- code contains reproducible code and a compiled pdf file that shows the results of the code.
- model has some Stan modeling results that were precomputed (see code).
- Paper contains the original paper. Note that it may not compile on your particular machine.

To extract R code from the Rmd file (in the code directory), use purl() in the knitr library.

.
├── R
│   ├── createStanDat.R
│   ├── createStanDatAcc.R
│   ├── gen_fake_lnorm.R
│   ├── gen_fake_lnorm2x2.R
│   ├── gen_fake_norm.R
│   ├── magnifytext.R
│   ├── multiplot.R
│   ├── plotpredictions.R
│   ├── plotresults.R
│   └── stan_results.R
├── README.md
├── code
│   ├── Uncertainty.Rmd
│   ├── Uncertainty.pdf
│   data
│   ├── DillonE1.txt
│   ├── JMVV2019replication.Rda
│   ├── Lago.csv
│   ├── Tucker.RData
│   ├── Wagers.Rdata
│   ├── agrmt_mismatch.Rda
│   ├── bayesfactors.txt
│   ├── data_model.Rda
│   ├── data_model_dillonrep.Rda
│   ├── gibsonwu2012data.txt
│   ├── lmer_estimates2.Rda
│   ├── lmer_estimates3.Rda
│   ├── posteriorsTargetMismatch.Rda
│   ├── public_article_data.txt
│   ├── remafit.Rda
│   ├── smallsamplesestimates.Rda
│   └── tvals.Rda
├── model
│   └── au_predicted_meansD13rep.Rda
└── paper
    ├── UncertaintyVasishthGelman2021Preprint.Rnw
    ├── bibio.bib
    ├── dgruyter.sty


You can extract the R code from the Rnw file by typing

    knitr::purl("")
