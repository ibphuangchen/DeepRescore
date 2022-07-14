
configfile:"/home/chuang8/deepRescore_config.yaml"

rule calc_basic_features:
	input:
		mzid="mzid/{sample}.mzid"
		spectrum="mzML/{sample}.mzML"
	output:
		"features/{sample}/features.txt"
	shell:
		"""
		java -jar {config["deepRescorePath"]}/bin/PDV-1.6.1.beta.features/PDV-1.6.1.beta.features-jar-with-dependencies.jar \
				-r {input.mzid} \
				-rt $result_type \
				-s {input.spetrum} \
				-st 1 \
				-i * \
				-k s \
				-o features/{wildcards.sample} \
				-a 0.05 \
				-c 0 \
				-decoy ${decoy_prefix} \
				-ft pdf \
				--features
		"""


rule pga_fdr_control:
	input:
		"features/{sample}/features.txt"
	output:
		rawPSM="rawPSM/{sample}_rawPSMs.txt",
		pgaPep="pga/{sample}_peptide_level/pga-peptideSummary.txt",
		pgaPsm="pga/{sample}_psm_level/pga-peptideSummary.txt"

	shell:
		"""
		Rscript {config["deepRescorePath"]}/bin/got_pga_input.R {input} {config['searchBy']} {output.rawPSM}
    Rscript {config["deepRescorePath"]}/bin/calculate_fdr.R {output.rawPSM} pga/{wildcards.sample} {config["deepRescorePath"]}/bin/protein.pro-ref.fasta rev_
		"""

rule generate_train_prediction_data:
	input:
		pgaPep="{sample}_peptide_level/pga-peptideSummary.txt",
		pgaPsm="{sample}_psm_level/pga-peptideSummary.txt",
		features="features/{sample}/features.txt",
		rawPSM="rawPSM/{sample}_rawPSMs"

	output:
		autoRT_tr="autoRT_train/{sample}_autoRT_tr.tsv",
		autoRT_pd="autoRT_prediction/{sample}_autoRT_pd.tsv",
		pdeep2_pd="pdeep2_prediction/{sample}_pdeep2_pd.tsv",
		pdeep2_pd_u="pdeep2_prediction/{sample}_pdeep2_pd_unique.tsv"
	shell:
		"""
		Rscript {config["deepRescorePath"]}/bin/got_train_prediction.R {input.pgaPep} \
        {input.pgaPsm} {input.features} {output.autoRT_tr} {output.autoRT_pd} \
				pdeep2_prediction/{wildcards.sample}_pdeep2 {input.rawPSM}
		"""

rule train_autoRT:
	input:
		"autoRT_train/{sample}_autoRT_tr.tsv",
	output:
		directory("autoRT_models/{sample}")
	shell:
		"""
		python /opt/AutoRT/autort.py train \
        -i {input} \
        -o {output} \
        -e 40 \
        -b 64 \
        -u m \
        -m /opt/AutoRT/models/base_models_PXD006109/model.json \
        -rlr \
        -n 10
    """

rule predict_autoRT:
	input:
		models="autoRT_models/{sample}/model.json",
		autoRT_pd="autoRT_prediction/{sample}_autoRT_pd.tsv"
	output:
		"autoRT_prediction/{sample}_autoRT_results.txt"
	shell:
		"""
		python /opt/AutoRT/autort.py predict \
        -t {input.autoRT_pd} \
        -s {input.models} \
        -o autoRT_prediction \
        -p {wildcards.sample}_autoRT_results
    
		"""


rule run_pdeep2:
	input:
		"pdeep2_prediction/{sample}_pdeep2_pd_unique.tsv"
	output:
		"pdeep2_prediction/{sample}_pdeep2_pd_results.txt"
	shell:
    """
    python /opt/pDeep2/predict.py -e $energy -i $instrument \
				-in {input} \
				-out {output} 
    """

rule process_pDeep2_results:
	input:
		pdeep2_pd="pdeep2_prediction/{sample}_pdeep2_pd.tsv"
		pdeep2_pd_res="pdeep2_prediction/{sample}_pdeep2_pd_results.txt",
		spectrum="mzML/{sample}.mzML",

	output:
		format_title="pdeep2_prediction/{sample}_format_titles.txt",
		pairs="pdeep2_prediction/{sample}_spectrum_pairs.txt",
		similarity="pdeep2_prediction/{sample}_similarity_SA.txt",

	threads:
		12
	shell:
		"""
    mv {input.pdeep2_pd_res} {input.pdeep2_pd_res}.mgf

    Rscript {config["deepRescorePath"]}/bin/format_pDeep2_titile.R {input.pdeep2_pd} \
				{input.spectrum} {output.format_title}

    java  -cp {config["deepRescorePath"]}/bin/PDV-1.6.1.beta.features/PDV-1.6.1.beta.features-jar-with-dependencies.jar PDVGUI.GenerateSpectrumTable \
				{output.format_title} {input.spectrum} {input.pdeep2_pd_res}.mgf {output.pairs} {config['searchBy']}

    mkdir -p pdeep2_prediction/{wildcards.sample}_sections 
		mkdir -p pdeep2_prediction/{wildcards.sample}.sections_results

    Rscript {config["deepRescorePath"]}/bin/similarity/devide_file.R {output.pairs} {threads} pdeep2_prediction/{wildcards.sample}_sections/
    
		for file in pdeep2_prediction/{wildcards.sample}_sections/*
    do
        name=`basename $file`
        Rscript {config["deepRescorePath"]}/bin/similarity/calculate_similarity_SA.R $file pdeep2_prediction/{wildcards.sample}_sections_results/{name}_results.txt &
    done
    wait
    awk 'NR==1 {header=\$_} FNR==1 && NR!=1 { \$_ ~ \$header getline; } {print}' pdeep2_prediction/{wildcards.sample}_sections_results/*_results.txt > {output.similarity}
		"""

rule generate_percolator_input:
	input:
		rawPSM="rawPSM/{sample}_rawPSMs.txt",
		features="features/{sample}/features.txt",
		autort_res="autoRT_prediction/{sample}_autoRT_results.txt",
		similarity="pdeep2_prediction/{sample}_similarity_SA.txt"
	output:
		"percolator_in/{sample}_format.pin"
	shell:
		"""
		Rscript {config["deepRescorePath"]}/bin/percolator/format_percolator_input.R {input.features} \
				{input.rawPSM} {input.autort_res} {input.similarity} {output} {config['searchBy']}
		"""

rule run_percolator:
	input:
		"percolator_in/{sample}_format.pin"
	output:
		pep="percolator_results/{sample}_pep.tsv",
		psm="percolator_results/{sample}_psm.tsv"
	shell:
		"""
		percolator -r {output.pep} -m {output.psm} {input}
		"""

rule generate_pdv_input:
	input:
		percolator_in="percolator_in/{sample}_format.pin",
		features="features/{sample}/features.txt"
	output:
		"DeepRescore_results/{sample}_pep.final.tsv",
		"DeepRescore_results/{sample}_psm.final.tsv"
	shell:
		"""
		Rscript {config["deepRescorePath"]}/bin/PDV/generate_pdv_input.R ./ {wildcards.sample} \
				{input.features} {input.percolator_in}
		"""
