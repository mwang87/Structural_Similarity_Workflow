import pandas as pd
import argparse
import smilesComparison_draft
import random
import json
import copy
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description='Run comparisons')
    parser.add_argument('input_smiles', type=str, help='Input file')
    parser.add_argument('output_filename', type=str, help='output_filename')
    parser.add_argument('--node_current', type=int, default=0, help='Current Node')
    parser.add_argument('--node_total', type=int, default=1, help='Total Nodes')
    
    random.seed(10)
    SUBSAMPLE_RATIO = 1

    args = parser.parse_args()
    smiles_df = pd.read_csv(args.input_smiles, sep='\t')

    print(smiles_df)

    all_smiles_list = smiles_df.to_dict(orient="records")
    all_smiles_list_subset = copy.deepcopy(all_smiles_list)

    output_list = []

    current_node = args.node_current
    total_nodes = min(args.node_total, len(all_smiles_list))

    if current_node >= total_nodes:
        return

    # Subsample the list
    all_smiles_list_subset = all_smiles_list_subset[current_node::total_nodes]

    for i, smiles1 in tqdm(enumerate(all_smiles_list_subset)):
        for j, smiles2 in enumerate(all_smiles_list):
            if i >= j:
                continue

            # Subsampling flag
            if random.randint(0, SUBSAMPLE_RATIO) != 0:
                continue

            # Doing the calculation
            comparison_result = smilesComparison_draft.smiles_comparisson(
                smiles1["SMILES"], smiles2["SMILES"], smiles1["ID"], smiles2["ID"], MCStimeout=10
            )

            try:
                json.dumps(comparison_result)
            except:
                print("ERROR SERIALIZING")
                print(comparison_result)
                continue

            output_list.append(comparison_result)
            
    with open(args.output_filename, 'w') as f:
        json.dump(output_list, f)


if __name__ == "__main__":
    main()
