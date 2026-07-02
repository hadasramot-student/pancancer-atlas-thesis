# my_scimilarity.py
#
# Ch.2 (Annotation) - Scimilarity cell-type annotation for one object,
# called between 04_object_annotation.R and 05_consensus_annotation.R.
# Embeds cells with the Scimilarity foundation model and assigns a cell
# type by weighted kNN vote against the reference atlas.
#
# Inputs:
#   --input  : h5ad written by 04_object_annotation.R (log-normalized counts
#              expected in adata.layers["logcounts"], raw counts in
#              adata.layers["counts"] if present)
#   --model-path : local path to a downloaded Scimilarity model checkpoint
# Outputs:
#   <temp_objects_dir>/<output> - h5ad with scimilarity_celltype,
#   scimilarity_vs2nd, scimilarity_celltype_filtered obs columns
#
# Usage:
#   conda run -n scimilarity python my_scimilarity.py \
#     --input step2.h5ad --project-name Azizi2018_Breast --output step3.h5ad \
#     --model-path /path/to/scimilarity/model_v1.1

import os
import scanpy as sc
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts
from optparse import OptionParser

def main():
    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-i", "--input",
                      dest="input",
                      help="input file name(required)"
    )
    parser.add_option("-o", "--output",
                    dest="output",
                    type="string",
                    help="output file name(required)"
    )
    parser.add_option("-m", "--model-path",
                      dest="model_path",
                      default=os.environ.get("SCIMILARITY_MODEL_PATH", "models/scimilarity/model_v1.1"),
                      help="Path to the scimilarity model (default: $SCIMILARITY_MODEL_PATH or models/scimilarity/model_v1.1)")
    parser.add_option("-p", "--project-name",
                      dest="project_name",
                      help="Project name for output file (required)")
    parser.add_option("-d", "--temp_objects_dir",
                      dest="temp_objects_dir",
                      default=os.environ.get("ATLAS_DATA_DIR", "data/temp_objects/"),
                      help="Output directory for results (default: $ATLAS_DATA_DIR or data/temp_objects/)")
    parser.add_option("-t", "--threshold",
                      dest="threshold",
                      type="float",
                      default=0.60,
                      help="Confidence threshold for filtering (default: 0.60)")
    parser.add_option("-k", "--knn",
                      dest="k",
                      type="int",
                      default=50,
                      help="Number of nearest neighbors (default: 50)")
    parser.add_option("-g", "--gpu",
                      dest="use_gpu",
                      action="store_true",
                      default=False,
                      help="Use GPU for computation")

    (options, args) = parser.parse_args()
    print(options)
    # Validate required arguments
    if not options.input or not options.project_name:
        parser.error("Both --input and --project-name are required")

    adata = sc.read_h5ad(options.temp_objects_dir + options.input)

    # (Strongly recommended) Ensure adata.X contains raw-ish counts before SCimilarity normalization.
    # If your counts are stored in a layer, use it:
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    ca = CellAnnotation(model_path=options.model_path, use_gpu=options.use_gpu)

    # Use your existing log-normalized layer as the main matrix
    adata.X = adata.layers["logcounts"].copy()

    # Align genes to the model's expected order (adds missing genes as zeros, reorders)
    adata = align_dataset(adata, ca.gene_order)

    # Embed + predict
    adata.obsm["X_scimilarity"] = ca.get_embeddings(adata.X)

    pred, nn_idxs, nn_dists, nn_stats = ca.get_predictions_knn(
        adata.obsm["X_scimilarity"],
        k=options.k,
        weighting=True,
    )

    adata.obs["scimilarity_celltype"] = pred.values
    adata.obs["scimilarity_vs2nd"] = nn_stats["vs2nd_weighted"]

    # Optional: confidence filter
    adata.obs["scimilarity_celltype_filtered"] = adata.obs["scimilarity_celltype"].where(
        adata.obs["scimilarity_vs2nd"] >= options.threshold, "Unknown"
    )
    save_path = options.temp_objects_dir + options.output
    adata.write_h5ad(save_path)
    print("Saved: " + save_path)

    os.remove((options.temp_objects_dir + options.input))
    print("removed: " + options.input)


if __name__ == "__main__":
    main()
