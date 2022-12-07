import os
import re


path = "data/h5ad"
files = os.listdir(path)

h5adfilesbool = [bool(re.search(".*\.h5ad", file)) for file in files]
h5adfiles = []
for i, bool in enumerate(h5adfilesbool):
    if bool:
        h5adfiles.append(files[i])

prefixes = [h5adfile.split(".")[0] for h5adfile in h5adfiles]



rule all:
    input:
       'data/figs/auc.jpeg'


        
rule subsample:
    input:
        data = "data/h5ad/{h5adfile}.h5ad"
    output:
        "data/subsampled/{h5adfile}.h5ad"
    params:
        script = "bin/subsample_cells.py"
    shell:
        """
        python {params.script} {input.data}
        """
        
        
rule normalize_depth:
    input:
        data = "data/subsampled/{h5adfile}"
    output:
        adata = "data/subsampled_normalized/{h5adfile}"
    params:
        script = "bin/normalize_depth.py"

    shell:
        """
        python {params.script} {input.data}
        """
        
rule average_UMIs:
    input:
        data = "data/subsampled_normalized/{h5adfile}.h5ad"
    output:
        organism_part_averages = "data/average_counts/{h5adfile}.csv"
    params:
        script = "bin/average_UMI.py"

    shell:
        """
        python {params.script} {input.data}
        """
        
        
rule create_pseudobulk:
    input:
        data = expand("data/average_counts/{h5adfile}.csv", h5adfile = prefixes)
    output:
        pseudobulk ="data/pseudobulk/pseudobulk.csv"

    shell:
        """
        array=($(ls data/average_counts))
        echo ${{array[0]}}
        
        touch {output.pseudobulk}
        
        head -n 1 {input.data}/${{array[0]}} > {output.pseudobulk}
        
        
        for file in ${{array[@]}}; do
            tail -n+2  {input.data}/${{file}} >> {output.pseudobulk}
        done
        """
        
rule run_EGAD:
    input:
        data="data/pseudobulk/pseudobulk.csv"
    output:
        egad="data/EAGD/EGAD.csv"

    script:
        "bin/EGAD.R"

rule vis_auc:
    input:
        data = "data/EAGD/EGAD.csv"
    output:
        'data/figs/auc.jpeg'
    params:
        script = 'bin/vis_auruc.py'
    shell:
        """
        python {params.script} {input.data}
        """