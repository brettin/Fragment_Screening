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
                       help='Clustering similarity cutoff (0.0-1.0, default: 0.4). Lower values create tighter clusters.')
    
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
    
    # Validate cutoff parameter
    if args.cutoff < 0.0 or args.cutoff > 1.0:
        print(f"Error: Cutoff must be between 0.0 and 1.0, got {args.cutoff}")
        return
    
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
        print(f"Clustering with similarity cutoff {args.cutoff} (distance threshold {1-args.cutoff:.3f})...")
    
    # Convert similarity cutoff to distance threshold for Butina clustering
    distance_cutoff = 1 - args.cutoff
    clusters = Butina.ClusterData(dists, len(fps), distance_cutoff, isDistData=True)
    
    if args.verbose:
        print(f"Created {len(clusters)} clusters")
    
    # ---------- Step 6: Generate pools ----------
    if args.verbose:
        print(f"Creating pools with maximal diversity (max size {args.pool_size})...")
    
    pool_list = []
    
    # Strategy: Create pools with maximally diverse compounds
    # Respect maximum pool size while ensuring ALL fragments are assigned
    
    # Get cluster sizes and create a distribution plan
    cluster_sizes = [len(cluster) for cluster in clusters]
    total_fragments = sum(cluster_sizes)
    num_clusters = len(clusters)
    
    if args.verbose:
        print(f"Distributing {total_fragments} fragments from {num_clusters} clusters")
    
    # Initialize pools as empty lists
    pools = []
    current_pool = []
    
    # Distribute fragments to maximize diversity within each pool
    # Take one fragment from each cluster before starting a new pool
    cluster_indices = [0] * num_clusters  # Track position in each cluster
    
    while True:
        # Try to add one fragment from each cluster to current pool
        added_fragment = False
        
        for cluster_idx in range(num_clusters):
            cluster = clusters[cluster_idx]
            
            # If we haven't used all fragments from this cluster
            if cluster_indices[cluster_idx] < len(cluster):
                # Check if adding this fragment would exceed max pool size
                if len(current_pool) >= args.pool_size:
                    # Current pool is full, finalize it and start a new one
                    pools.append(current_pool)
                    current_pool = []
                
                # Add fragment to current pool
                current_pool.append(cluster[cluster_indices[cluster_idx]])
                cluster_indices[cluster_idx] += 1
                added_fragment = True
        
        # If we couldn't add any fragments, we're done
        if not added_fragment:
            break
    
    # Add any remaining fragments to the last pool
    if current_pool:
        pools.append(current_pool)
    
    # Ensure ALL fragments are assigned to pools
    # Check if any fragments are missing
    assigned_fragments = set()
    for pool in pools:
        for frag in pool:
            assigned_fragments.add(frag)
    
    # Find any unassigned fragments
    all_fragments = set(range(len(df)))
    unassigned_fragments = all_fragments - assigned_fragments
    
    if unassigned_fragments:
        if args.verbose:
            print(f"Found {len(unassigned_fragments)} unassigned fragments, distributing them...")
        
        # Distribute unassigned fragments to existing pools or create new ones
        unassigned_list = list(unassigned_fragments)
        
        # For each unassigned fragment, find which cluster it belongs to
        for frag in unassigned_list:
            # Find which cluster this fragment belongs to
            frag_cluster = None
            for cluster_idx, cluster in enumerate(clusters):
                if frag in cluster:
                    frag_cluster = cluster_idx
                    break
            
            if frag_cluster is None:
                if args.verbose:
                    print(f"Warning: Fragment {frag} not found in any cluster")
                continue
            
            # Find pools that don't contain any fragments from the same cluster
            compatible_pools = []
            for pool_idx, pool in enumerate(pools):
                pool_has_similar = False
                for pool_frag in pool:
                    # Check if this pool fragment is from the same cluster
                    for cluster_idx, cluster in enumerate(clusters):
                        if pool_frag in cluster and cluster_idx == frag_cluster:
                            pool_has_similar = True
                            break
                    if pool_has_similar:
                        break
                
                if not pool_has_similar and len(pool) < args.pool_size:
                    compatible_pools.append(pool_idx)
            
            # Try to add to a compatible pool
            added = False
            for pool_idx in compatible_pools:
                if len(pools[pool_idx]) < args.pool_size:
                    pools[pool_idx].append(frag)
                    added = True
                    break
            
            # If couldn't add to existing pools, create a new pool
            if not added:
                pools.append([frag])
    
    # Create pool assignments
    for pool_idx, pool_fragments in enumerate(pools):
        for frag in pool_fragments:
            pool_list.append({'PoolID': f'Pool{pool_idx + 1}', 'FragmentID': df.iloc[frag]['FragmentID']})
    
    # Save pool assignment file
    output_pools_file = f"{args.output_prefix}_pools.csv"
    df_pools = pd.DataFrame(pool_list)
    df_pools.to_csv(output_pools_file, index=False)
    
    if args.verbose:
        print(f"Created {len(pools)} pools with variable sizes")
        pool_sizes = [len(pool) for pool in pools]
        print(f"Pool sizes: min={min(pool_sizes)}, max={max(pool_sizes)}, avg={sum(pool_sizes)/len(pool_sizes):.1f}")
        print(f"All pools respect maximum size of {args.pool_size}")
        print(f"Total fragments assigned: {len(pool_list)} out of {total_fragments}")
        if len(pool_list) == total_fragments:
            print("✅ All fragments successfully assigned to pools!")
        else:
            print(f"⚠️  Warning: {total_fragments - len(pool_list)} fragments not assigned!")
        print(f"Saved pool assignments to: {output_pools_file}")
    
    print("Processing complete!")
    if args.props.lower() != 'none':
        print(f"Enriched file saved as: {output_props_file}")
    print(f"Pools file saved as: {output_pools_file}")

if __name__ == "__main__":
    main()

