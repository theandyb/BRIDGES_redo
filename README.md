# Introduction

The point of this repository is to re-do my analyses of the germline mutation rate using singletons as a proxy. I want to ensure that my results using the BRIDGES dataset are both correct and reproducible. To that end I am making a few changes:

1. I'd like to stream-line the pipeline as much as possible. While a tool like snakemake would ensure this, for now I'm going to maintain a connected set of "steps" encoded in scripts to generate my results. (I really should learn snakemake though).
2. If I get run over by a bus, someone else should be able to read this document, and based on that (and access to the data) they should be able to reproduce my results.

# Step 1: Generate singleton files
