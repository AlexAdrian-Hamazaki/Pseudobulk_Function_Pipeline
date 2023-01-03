configfile: "config.yaml"


import os
import re


files = os.listdir(config["path"])
print(files)

h5adfilesbool = [bool(re.search(".*\.h5ad", file)) for file in files]
h5adfiles = []
for i, bool in enumerate(h5adfilesbool):
    if bool:
        h5adfiles.append(files[i])

prefixes = [h5adfile.split(".")[0] for h5adfile in h5adfiles]


rule all:
    input:
       'data/figs/auc.jpeg',
       "data/average_counts/TS_Liver.csv"


        

        
rule normalize_depth:
    input:
        data = expand("{path}/{{h5adfile}}.h5ad", path=config["path"])
    output:
        adata = expand("{path}/normalized/{{h5adfile}}.h5ad", path=config["path"])
    params:
        script = "bin/normalize_depth.py"

    shell:
        """
        python {params.script} {input.data}
        """
        
rule average_UMIs:
    input:
        data = expand("{path}/normalized/{{h5adfile}}.h5ad", path=config["path"])
    output:
        organism_part_averages = "data/average_counts/{h5adfile}.csv"
    params:
        script = "bin/average_UMI.py",
        cell_type_column = config["UserInput_CellTypeColumn"]

    shell:
        """
        python {params.script} {input.data} {params.cell_type_column}
        """
        
        
rule create_pseudobulk:
    input:
        data = expand("data/average_counts/{h5adfile}.csv", h5adfile = prefixes)
    output:
        pseudobulk ="data/pseudobulk/pseudobulk.csv"

    shell:
        """


        touch {output.pseudobulk}
        
        head -n 1 {input.data[0]} > {output.pseudobulk}
        
        array2=({input.data})
        
        for file in ${{array2}}; do
            echo ${{file}}
            tail -n+2  ${{file}} >> {output.pseudobulk}
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