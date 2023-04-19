version 1.0

struct RuntimeAttr {
  Int? cpu
  Int? memory_gb
  String? docker
}

workflow analysisFlow {
    input {
        Object params
        RuntimeAttr runtime_analysisflow
    }

    call calculate_memory {
        input:
            params = params,
            runtime_attr = runtime_analysisflow
    }

    call analysis_flow {
        input:
            params = params,
            runtime_attr_override = runtime_analysisflow,
            mem_on_mtx = calculate_memory.mem_on_mtx
    }

    output {
        File tar_gz_output = analysis_flow.tar_gz_output
    }
}

task calculate_memory {
    input {
        Object params
        RuntimeAttr runtime_attr
    }

    command {
        set -euo pipefail

        path_list="$(echo ~{select_first([params.path, ""])} | tr ',' ' ')"
        total_cells=0
        total_genes=0
        for path in $path_list; do
            if [ -f "$path/barcodes.tsv" ]; then
                cell_num=$(wc -l < "$path/barcodes.tsv")
            elif [ -f "$path/barcodes.tsv.gz" ]; then
                cell_num=$(zcat "$path/barcodes.tsv.gz" | wc -l)
            else
                echo "Error: barcodes.tsv or barcodes.tsv.gz not found in $path"
                exit 1
            fi

            if [ -f "$path/genes.tsv" ]; then
                gene_num=$(wc -l < "$path/genes.tsv")
            elif [ -f "$path/genes.tsv.gz" ]; then
                gene_num=$(zcat "$path/genes.tsv.gz" | wc -l)
            elif [ -f "$path/features.tsv.gz" ]; then
                gene_num=$(zcat "$path/features.tsv.gz" | wc -l)
            else
                echo "Error: genes.tsv, genes.tsv.gz or features.tsv.gz not found in $path"
                exit 1
            fi

            total_cells=$((total_cells + cell_num))
            total_genes=$((total_genes + gene_num))
        done
        echo "total_cells: $total_cells"
        echo "total_genes: $total_genes"
        mem_on_mtx=$((total_cells * total_genes))
        echo

    }

    output {
        Int mem_on_mtx = ceil(mem_on_mtx * 0.00000003 + 2)
    }

    RuntimeAttr runtime_attr_default = object {
        cpu: 2,
        memory_gb: 2
    }
    Int cpu = select_first([runtime_attr.cpu, runtime_attr_default.cpu])
    Int memory_gb = select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])
    String docker = select_first([runtime_attr.docker, runtime_attr_default.docker])

    runtime {
        cpu: cpu
        memory: memory_gb
        docker: docker
    }
}

task analysis_flow {
    input {
        Object params
        RuntimeAttr runtime_attr_override
        Int mem_on_mtx
    }

    command {
        set -euo pipefail
        python /tmp/SPXcanpy.py  \
            --path "~{params.path}" \
            --spname "~{params.spname}" \
            --gname "~{params.gname}" \
            --datatype "~{params.datatype}" \
            --mingene "~{params.mingene}" \
            --maxgene "~{params.maxgene}" \
            --minumi "~{params.minumi}" \
            --maxumi "~{params.maxumi}" \
            --mtfiler "~{params.mtfiler}" \
            --mincell "~{params.mincell}" \
            --subscale "~{params.subscale}" \
            --clust_method "~{params.clust_method}" \
            --dim "~{params.dim}" \
            --resolution "~{params.resolution}" \
            --rmbatch "~{params.rmbatch}" \
            --batchname "~{params.batchname}" \
            --logfc "~{params.logfc}" \
            --minpct  "~{params.minpct}" \
            --color  "~{params.color}" \
            --dendrogram  "~{params.dendrogram}" \
            --outdir  "~{params.outdir}" \
            --prefix  "~{params.prefix}" \
        tar -czvf output.tar.gz ~{params.outdir}
        echo 'scanpyqc done'
    }

    RuntimeAttr runtime_attr_default = object {
        cpu: 4,
        memory_gb: mem_on_mtx
    }

    Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
    String docker = select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    runtime {
        cpu: cpu
        memory: memory_gb
        docker: docker
    }
    output {
        File tar_gz_output = "scanpyqc1_out/scanpyqcv1.tar.gz"
    }
}
