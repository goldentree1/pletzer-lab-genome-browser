#!/bin/bash

# EXAMPLE for PA14 (I versus Un)::

bigwigAverage --binSize 1 --bigwigs \
public/data/GCF_000014625/PA14_Un_1.reheadered.bw \
public/data/GCF_000014625/PA14_Un_2.reheadered.bw \
public/data/GCF_000014625/PA14_Un_3.reheadered.bw \
-o public/data/GCF_000014625/PA14_Un_1-to-3-average__binsize1.reheadered.bw


bigwigAverage --binSize 1 --bigwigs \
public/data/GCF_000014625/PA14_I_1.reheadered.bw \
public/data/GCF_000014625/PA14_I_2.reheadered.bw \
public/data/GCF_000014625/PA14_I_3.reheadered.bw \
-o public/data/GCF_000014625/PA14_I_1-to-3-average__binsize1.reheadered.bw
