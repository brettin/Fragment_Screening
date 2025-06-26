import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.ML.Cluster import Butina
from rdkit.Chem import AllChem
from rdkit.Chem import QED
import numpy as np
import argparse
import random
from rdkit import DataStructs

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Cluster molecular fragments and create screening pools',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cluster_fragments.py
  python cluster_fragments.py -i my_fragments.csv -o my_results
  python cluster_fragments.py -b 512 -r 1 -c 0.3 --props MW,logP
  python cluster_fragments.py -b 2048 -r 3 -c 0.5 --props all -p 2 -s 42
        """
    )
    
    # Input/Output options
    parser.add_argument('-i', '--input-file', 
                       default='fragments.csv',
                       help='Input CSV file path (default: fragments.csv)')
    parser.add_argument('-o', '--output-prefix', 
                       default='fragment',
                       help='Output file prefix (default: fragment)')
    parser.add_argument('-v', '--verbose', 
                       action='store_true',
                       help='Enable verbose output')
    
    # Fingerprint parameters
    parser.add_argument('-b', '--fingerprint-bits', 
                       type=int, default=1024,
                       help='Morgan fingerprint bit count (default: 1024)')
    parser.add_argument('-r', '--fingerprint-radius', 
                       type=int, default=2,
                       help='Morgan fingerprint radius (default: 2)')
    
    # Clustering parameters
    parser.add_argument('-c', '--cutoff', 
                       type=float, default=0.4,
                       help='Clustering similarity cutoff (default: 0.4)')
    
    # Pooling parameters
    parser.add_argument('-p', '--pool-size', 
                       type=int, default=3,
                       help='Number of fragments per pool (default: 3)')
    parser.add_argument('-s', '--random-seed', 
                       type=int, default=None,
                       help='Random seed for reproducible results')
    
    # Property calculation
    parser.add_argument('--props', '--properties', 
                       default='all',
                       help='Molecular properties to calculate (comma-separated or "all"/"none")')
    
    return parser.parse_args()

def compute_properties(smiles, properties_to_calculate):
    """Compute molecular properties based on requested properties."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [None] * len(properties_to_calculate)
    
    results = []
    for prop in properties_to_calculate:
        if prop == 'MW':
            results.append(Descriptors.MolWt(mol))
        elif prop == 'logP':
            results.append(Crippen.MolLogP(mol))
        elif prop == 'HBD':
            results.append(Lipinski.NumHDonors(mol))
        elif prop == 'HBA':
            results.append(Lipinski.NumHAcceptors(mol))
        elif prop == 'TPSA':
            results.append(rdMolDescriptors.CalcTPSA(mol))
        elif prop == 'RotBonds':
            results.append(Lipinski.NumRotatableBonds(mol))
        elif prop == 'Solubility':
            # ESOL-like solubility estimate
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            rotb = Lipinski.NumRotatableBonds(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            logS = 0.16 - 0.63*logp - 0.0062*mw + 0.066*hbd - 0.74*rotb - 0.006*tpsa
            solubility_mg_per_mL = 10**logS * mw
            results.append(solubility_mg_per_mL)
        else:
            results.append(None)
    
    return results

def main():
    """Main function to process fragments."""
    args = parse_arguments()
    
    # Set random seed if provided
    if args.random_seed is not None:
        random.seed(args.random_seed)
        np.random.seed(args.random_seed)
    
    if args.verbose:
        print(f"Processing input file: {args.input_file}")
        print(f"Output prefix: {args.output_prefix}")
        print(f"Fingerprint bits: {args.fingerprint_bits}")
        print(f"Fingerprint radius: {args.fingerprint_radius}")
        print(f"Clustering cutoff: {args.cutoff}")
        print(f"Pool size: {args.pool_size}")
        print(f"Properties to calculate: {args.props}")
    
    # ---------- Step 1: Load input file ----------
    try:
        df = pd.read_csv(args.input_file)
        if args.verbose:
            print(f"Loaded {len(df)} fragments from {args.input_file}")
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return
    
    # Validate input format
    required_columns = ['FragmentID', 'SMILES']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: Missing required columns: {missing_columns}")
        print("Input file must contain: FragmentID, SMILES")
        return
    
    # ---------- Step 2: Calculate properties ----------
    if args.props.lower() != 'none':
        if args.props.lower() == 'all':
            properties_to_calculate = ['MW', 'logP', 'HBD', 'HBA', 'TPSA', 'RotBonds', 'Solubility']
        else:
            properties_to_calculate = [p.strip() for p in args.props.split(',')]
        
        if args.verbose:
            print(f"Calculating properties: {properties_to_calculate}")
        
        props = df['SMILES'].apply(lambda x: compute_properties(x, properties_to_calculate))
        df_props = pd.DataFrame(props.tolist(), columns=properties_to_calculate)
        df = pd.concat([df, df_props], axis=1)
        
        # Save enriched fragment file
        output_props_file = f"{args.output_prefix}_with_properties.csv"
        df.to_csv(output_props_file, index=False)
        if args.verbose:
            print(f"Saved enriched data to: {output_props_file}")
    
    # ---------- Step 3: Generate fingerprints for clustering ----------
    if args.verbose:
        print("Generating molecular fingerprints...")
    
    mols = [Chem.MolFromSmiles(s) for s in df['SMILES']]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, args.fingerprint_radius, nBits=args.fingerprint_bits) 
           for m in mols]
    
    # ---------- Step 4: Compute distance matrix ----------
    if args.verbose:
        print("Computing similarity matrix...")
    
    dists = []
    for i in range(len(fps)):
        for j in range(i):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            dists.append(1 - sim)
    
    # ---------- Step 5: Perform Butina clustering ----------
    if args.verbose:
        print(f"Clustering with cutoff {args.cutoff}...")
    
    clusters = Butina.ClusterData(dists, len(fps), args.cutoff, isDistData=True)
    
    if args.verbose:
        print(f"Created {len(clusters)} clusters")
    
    # ---------- Step 6: Generate pools ----------
    if args.verbose:
        print(f"Creating pools of size {args.pool_size}...")
    
    pool_list = []
    pool_counter = 1
    
    # Flatten clusters into list of fragments
    clustered_fragments = [list(cluster) for cluster in clusters]
    all_indices = sum(clustered_fragments, [])
    
    # Random shuffle to increase mixing
    random.shuffle(all_indices)
    
    for i in range(0, len(all_indices), args.pool_size):
        pool_fragments = all_indices[i:i+args.pool_size]
        for frag in pool_fragments:
            pool_list.append({'PoolID': f'Pool{pool_counter}', 'FragmentID': df.iloc[frag]['FragmentID']})
        pool_counter += 1
    
    # Save pool assignment file
    output_pools_file = f"{args.output_prefix}_pools.csv"
    df_pools = pd.DataFrame(pool_list)
    df_pools.to_csv(output_pools_file, index=False)
    
    if args.verbose:
        print(f"Created {pool_counter-1} pools")
        print(f"Saved pool assignments to: {output_pools_file}")
    
    print("Processing complete!")
    if args.props.lower() != 'none':
        print(f"Enriched file saved as: {output_props_file}")
    print(f"Pools file saved as: {output_pools_file}")

if __name__ == "__main__":
    main()

