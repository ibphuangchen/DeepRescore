#!/usr/bin/env sh

set -e

snakemake -j 500 --cluster-config ./cluster.jason \
					--rerun-incomplete \
					--use-singularity --singularity-args '-B /rsrch3,/home/chuang8' \
		          --cluster 'bsub -n {cluster.n} -q {cluster.queue} -W {cluster.time} -M {cluster.memory} -R {cluster.resources} -o {cluster.output} -e {cluster.error}' $@

