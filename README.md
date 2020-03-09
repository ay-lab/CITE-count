[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

CITE-count (Totalseq A module compatible with 10X data)
============

* Motivation: As a popular antibody labeling technic, totalseq A has been applied in CITE-Seq and LEAP-Seq based studies, yet the 10X software has very limited support for totalseq A data preprocessing. Here we provided a pipeline that preprosses raw totalseq A fastq files and conducts QC of the results. The pipeline can also import other totalseq data, merge with RNA-seq result and export seamlessly to Scanpy pipeline.

Fastq preprocessing background

![alt text](./img/Pre_processing.png "Totalseq A 3' capture preprocessing")