import sys, os
import pickle
import time
import numpy as np
import pandas as pd

if __name__=="__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Run BHtsne on a given compressed Pandas DataFrame.')
    parser.add_argument('input', type=str, 
                       help='Input dataframe filename. Loaded with load_df.')
    parser.add_argument('output', type=str, 
                       help='Output dataframe filename. Saved with save_df.')
    parser.add_argument('-s', '--seed', type=int, default=-1, dest='seed',
                       help='Random seed used ')
    parser.add_argument('-p', '--perplexity', type=int, default=50, dest='perplexity',
                       help='Perplexity setting')
    parser.add_argument('-d', '--bhtsne_dir', type=str, default="/n/home15/adrianveres/software/bhtsne", dest='bhtsne_dir',
                       help='Package directory')
    parser.add_argument('-i', '--max_iter', type=int, default=1000, dest='max_iter',
                       help='Maximum iterations')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose')
    parser.add_argument('-r', '--dims', type=int, default=-1, dest='input_dims')

    args = parser.parse_args()

    sys.path.append(args.bhtsne_dir)
    from bhtsne import run_bh_tsne

    with np.load(args.input) as f:
        input_df = pd.DataFrame(**f)

    if args.input_dims > 0:
        input_df = input_df.iloc[:, :args.input_dims]

    tsne_proj = run_bh_tsne(input_df.values, perplexity=args.perplexity, randseed=args.seed,
        verbose=args.verbose,
        max_iter=args.max_iter,
        use_pca=False,
        tmp_dir_prefix=os.path.dirname(args.input)+'/')

    tsne_proj = pd.DataFrame(tsne_proj, index=input_df.index, columns=np.arange(2))
    np.savez_compressed(args.output, data=tsne_proj.values, index=tsne_proj.index.values, columns=tsne_proj.columns.values)
