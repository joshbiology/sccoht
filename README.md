# The ATPase module of mammalian SWI/SNF family complexes mediates subcomplex identity and catalytic activity-independent genomic targeting

## Repository for Pan *et al*., Nature Genetics, 2019

Code repository for reproducing analyses and figures in this report. Credit to Robin Meyers for the ProjectTemplate setup from our [previous repository](https://github.com/robinmeyers/pan-meyers-et-al).

## Setup

### Clone this repository

From the command line:

```
$ git clone https://github.com/joshbiology/sccoht
```

### Open an R session from this directory

This can be done a variety of ways.

```
$ cd sccoht
$ R
```

or if you use RStudio

```
$ cd sccoht
$ open sccoht.Rproj
```

### Execute the SetupProject script

From the R console, run the SetupProject.R script. This will install any missing packages and download data files from the following figshare records:


1. [All raw data](https://doi.org/10.6084/m9.figshare.7801718) -- all raw and cached data used in the paper. 

This repo is a superset that contains the two peak repositories reference in the supplemental info of the paper:

2. [Raw peaks](https://doi.org/10.6084/m9.figshare.6965498) -- raw MACS 2.0 derived peaks from the ChIP-seq data used in the paper.

3. [Derived peaksets](https://doi.org/10.6084/m9.figshare.6965510) -- overlapped peaksets that are used in the paper to define complex groups across conditions. This includes residual complex positioning, activity-(in)dependent localization, BAF/PBAF etc. 

SetupProject.R also executes the scripts in the "munge" folder, which will save munged datasets in the "cache" folder as Rdata objects. 

```
> source("./R/SetupProject.R")
```

This will take a while to complete. Once it does, we recommend restarting your R session before continuing.

*Note: downloading the data files from figshare for the first time requires a web browser to authenticate.* If running on a remote machine without this capability, execute the following commands in an R console locally.

```
> library(rfigshare)
> figshare_article <- fs_details(6005297)
```

This will open a web browser session, log in through figshare, and generate a `.httr-oauth` file in your working directory. Copy this file to the project directory on the remote machine. Running the `SetupProject.R` script as above should no longer require the authentication in a web browser.
