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
from itertools import combinations
import heapq
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Advanced fragment clustering and diversity-optimized pooling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cluster_fragments_v2.py -i fragments.csv -v
  python cluster_fragments_v2.py -i fragments.csv -p 5 -c 0.3 --diversity-method greedy
  python cluster_fragments_v2.py -i fragments.csv -p 10 --diversity-method kmeans --iterations 100
        """
    )
    
    # Input/Output options
    parser.add_argument('-i', '--input-file', 
                       default='fragments.csv',
                       help='Input CSV file path (default: fragments.csv)')
    parser.add_argument('-o', '--output-prefix', 
                       default='fragment_v2',
                       help='Output file prefix (default: fragment_v2)')
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
                       help='Clustering similarity cutoff (0.0-1.0, default: 0.4)')
    
    # Pooling parameters
    parser.add_argument('-p', '--pool-size', 
                       type=int, default=5,
                       help='Target number of fragments per pool (default: 5)')
    parser.add_argument('-s', '--random-seed', 
                       type=int, default=None,
                       help='Random seed for reproducible results')
    
    # Advanced diversity optimization
    parser.add_argument('--diversity-method', 
                       choices=['greedy', 'kmeans', 'genetic'],
                       default='greedy',
                       help='Diversity optimization method (default: greedy)')
    parser.add_argument('--iterations', 
                       type=int, default=50,
                       help='Number of optimization iterations (default: 50)')
    parser.add_argument('--diversity-weight', 
                       type=float, default=1.0,
                       help='Weight for diversity vs size balance (default: 1.0)')
    parser.add_argument('--min-pool-size', 
                       type=int, default=3,
                       help='Minimum fragments per pool (default: 3)')
    
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

def calculate_diversity_score(fragments, similarity_matrix):
    """Calculate diversity score for a set of fragments."""
    if len(fragments) < 2:
        return 0.0
    
    # Calculate average pairwise dissimilarity (1 - similarity)
    total_dissimilarity = 0.0
    count = 0
    
    for i, j in combinations(fragments, 2):
        total_dissimilarity += (1.0 - similarity_matrix[i][j])
        count += 1
    
    return total_dissimilarity / count if count > 0 else 0.0

def greedy_diversity_pooling(fragments, similarity_matrix, target_pool_size, min_pool_size, diversity_weight=1.0):
    """
    Greedy algorithm to create pools with maximal diversity.
    
    Strategy:
    1. Start with the most diverse pair of fragments
    2. Iteratively add the fragment that maximizes pool diversity
    3. Create new pools when current pool reaches target size
    4. Ensure minimum pool size constraint
    """
    n_fragments = len(fragments)
    unassigned = set(range(n_fragments))
    pools = []
    
    while unassigned:
        if len(unassigned) < min_pool_size:
            # Create a small final pool
            current_pool = list(unassigned)
            pools.append(current_pool)
            break
        
        # Start new pool with the most diverse pair
        if len(unassigned) >= 2:
            best_pair = None
            best_diversity = -1
            
            for i, j in combinations(unassigned, 2):
                diversity = 1.0 - similarity_matrix[i][j]
                if diversity > best_diversity:
                    best_diversity = diversity
                    best_pair = (i, j)
            
            current_pool = list(best_pair)
            unassigned -= set(best_pair)
        else:
            current_pool = [unassigned.pop()]
        
        # Greedily add fragments to maximize diversity
        while len(current_pool) < target_pool_size and unassigned:
            best_fragment = None
            best_score = -1
            
            for frag in unassigned:
                # Calculate diversity score if we add this fragment
                test_pool = current_pool + [frag]
                diversity_score = calculate_diversity_score(test_pool, similarity_matrix)
                
                # Balance diversity with pool size preference
                size_penalty = len(test_pool) / target_pool_size
                combined_score = diversity_score * diversity_weight - size_penalty
                
                if combined_score > best_score:
                    best_score = combined_score
                    best_fragment = frag
            
            if best_fragment is not None:
                current_pool.append(best_fragment)
                unassigned.remove(best_fragment)
            else:
                break
        
        pools.append(current_pool)
    
    return pools

def kmeans_diversity_pooling(fragments, similarity_matrix, target_pool_size, min_pool_size, n_iterations=50):
    """
    K-means inspired approach for diversity-based pooling.
    
    Strategy:
    1. Initialize pool centers randomly
    2. Assign fragments to nearest (most similar) center
    3. Recalculate centers as most diverse representatives
    4. Iterate until convergence or max iterations
    """
    n_fragments = len(fragments)
    n_pools = max(1, n_fragments // target_pool_size)
    
    # Initialize pool centers randomly
    centers = random.sample(range(n_fragments), min(n_pools, n_fragments))
    
    for iteration in range(n_iterations):
        # Assign fragments to nearest center
        pool_assignments = defaultdict(list)
        for frag in range(n_fragments):
            if frag in centers:
                continue
            
            # Find the center with highest similarity
            best_center = None
            best_similarity = -1
            for center in centers:
                if center < frag:
                    similarity = similarity_matrix[center][frag]
                else:
                    similarity = similarity_matrix[frag][center]
                
                if similarity > best_similarity:
                    best_similarity = similarity
                    best_center = center
            
            pool_assignments[best_center].append(frag)
        
        # Recalculate centers as most diverse representatives
        new_centers = []
        for center in centers:
            pool_members = pool_assignments[center] + [center]
            
            # Find the fragment that is most similar to all others (best representative)
            best_representative = center
            best_avg_similarity = 0
            
            for member in pool_members:
                avg_similarity = 0
                count = 0
                for other in pool_members:
                    if member != other:
                        if member < other:
                            avg_similarity += similarity_matrix[member][other]
                        else:
                            avg_similarity += similarity_matrix[other][member]
                        count += 1
                
                if count > 0:
                    avg_similarity /= count
                    if avg_similarity > best_avg_similarity:
                        best_avg_similarity = avg_similarity
                        best_representative = member
            
            new_centers.append(best_representative)
        
        # Check for convergence
        if set(centers) == set(new_centers):
            break
        
        centers = new_centers
    
    # Convert assignments to pools
    pools = []
    for center in centers:
        pool = pool_assignments[center] + [center]
        if len(pool) >= min_pool_size:
            pools.append(pool)
        else:
            # Merge small pools
            for existing_pool in pools:
                if len(existing_pool) + len(pool) <= target_pool_size:
                    existing_pool.extend(pool)
                    break
            else:
                pools.append(pool)
    
    return pools

def genetic_diversity_pooling(fragments, similarity_matrix, target_pool_size, min_pool_size, n_iterations=50):
    """
    Genetic algorithm approach for diversity optimization.
    
    Strategy:
    1. Generate random pool assignments
    2. Evaluate fitness based on diversity and size constraints
    3. Crossover and mutate to create new generations
    4. Select best solutions
    """
    n_fragments = len(fragments)
    population_size = 20
    
    def create_random_assignment():
        """Create a random pool assignment."""
        assignment = list(range(n_fragments))
        random.shuffle(assignment)
        
        # Split into pools
        pools = []
        for i in range(0, n_fragments, target_pool_size):
            pool = assignment[i:i+target_pool_size]
            if len(pool) >= min_pool_size:
                pools.append(pool)
        
        return pools
    
    def evaluate_fitness(pools):
        """Evaluate fitness of a pool assignment."""
        total_diversity = 0
        size_penalty = 0
        
        for pool in pools:
            if len(pool) >= min_pool_size:
                diversity = calculate_diversity_score(pool, similarity_matrix)
                total_diversity += diversity * len(pool)
                
                # Penalize deviation from target size
                size_diff = abs(len(pool) - target_pool_size)
                size_penalty += size_diff
        
        return total_diversity - size_penalty
    
    # Initialize population
    population = [create_random_assignment() for _ in range(population_size)]
    
    for generation in range(n_iterations):
        # Evaluate fitness
        fitness_scores = [(evaluate_fitness(pools), pools) for pools in population]
        fitness_scores.sort(reverse=True)
        
        # Keep best half
        best_pools = [pools for _, pools in fitness_scores[:population_size//2]]
        
        # Create new generation through crossover and mutation
        new_population = best_pools.copy()
        
        while len(new_population) < population_size:
            # Crossover
            parent1, parent2 = random.sample(best_pools, 2)
            child = parent1[:len(parent1)//2] + parent2[len(parent2)//2:]
            
            # Mutation: randomly reassign some fragments
            if random.random() < 0.1:
                for pool in child:
                    if len(pool) > 1:
                        # Move a random fragment to a different pool
                        frag = random.choice(pool)
                        pool.remove(frag)
                        target_pool = random.choice(child)
                        target_pool.append(frag)
            
            new_population.append(child)
        
        population = new_population
    
    # Return best solution
    best_fitness, best_pools = max([(evaluate_fitness(pools), pools) for pools in population])
    return best_pools

def main():
    """Main function to process fragments with advanced diversity optimization."""
    args = parse_arguments()
    
    # Set random seed if provided
    if args.random_seed is not None:
        random.seed(args.random_seed)
        np.random.seed(args.random_seed)
    
    # Validate parameters
    if args.cutoff < 0.0 or args.cutoff > 1.0:
        print(f"Error: Cutoff must be between 0.0 and 1.0, got {args.cutoff}")
        return
    
    if args.verbose:
        print(f"Advanced Fragment Clustering and Diversity Optimization")
        print(f"Processing input file: {args.input_file}")
        print(f"Output prefix: {args.output_prefix}")
        print(f"Diversity method: {args.diversity_method}")
        print(f"Target pool size: {args.pool_size}")
        print(f"Min pool size: {args.min_pool_size}")
        print(f"Optimization iterations: {args.iterations}")
    
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
    
    # ---------- Step 3: Generate fingerprints ----------
    if args.verbose:
        print("Generating molecular fingerprints...")
    
    mols = [Chem.MolFromSmiles(s) for s in df['SMILES']]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, args.fingerprint_radius, nBits=args.fingerprint_bits) 
           for m in mols]
    
    # ---------- Step 4: Compute similarity matrix ----------
    if args.verbose:
        print("Computing similarity matrix...")
    
    n_fragments = len(fps)
    similarity_matrix = np.zeros((n_fragments, n_fragments))
    
    for i in range(n_fragments):
        for j in range(i+1):
            if i == j:
                similarity_matrix[i][j] = 1.0
            else:
                sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                similarity_matrix[i][j] = sim
                similarity_matrix[j][i] = sim
    
    # ---------- Step 5: Optional clustering for initial grouping ----------
    if args.cutoff < 1.0:
        if args.verbose:
            print(f"Performing initial clustering with cutoff {args.cutoff}...")
        
        dists = []
        for i in range(n_fragments):
            for j in range(i):
                dists.append(1 - similarity_matrix[i][j])
        
        distance_cutoff = 1 - args.cutoff
        clusters = Butina.ClusterData(dists, n_fragments, distance_cutoff, isDistData=True)
        
        if args.verbose:
            print(f"Created {len(clusters)} initial clusters")
    else:
        # No clustering, treat each fragment as its own cluster
        clusters = [[i] for i in range(n_fragments)]
    
    # ---------- Step 6: Advanced diversity-based pooling ----------
    if args.verbose:
        print(f"Creating pools using {args.diversity_method} optimization...")
    
    # Flatten clusters for diversity optimization
    all_fragments = []
    for cluster in clusters:
        all_fragments.extend(cluster)
    
    # Apply selected diversity optimization method
    if args.diversity_method == 'greedy':
        pools = greedy_diversity_pooling(all_fragments, similarity_matrix, 
                                       args.pool_size, args.min_pool_size, args.diversity_weight)
    elif args.diversity_method == 'kmeans':
        pools = kmeans_diversity_pooling(all_fragments, similarity_matrix, 
                                       args.pool_size, args.min_pool_size, args.iterations)
    elif args.diversity_method == 'genetic':
        pools = genetic_diversity_pooling(all_fragments, similarity_matrix, 
                                        args.pool_size, args.min_pool_size, args.iterations)
    else:
        print(f"Error: Unknown diversity method '{args.diversity_method}'")
        return
    
    # ---------- Step 7: Create pool assignments ----------
    pool_list = []
    total_fragments_assigned = 0
    
    for pool_idx, pool_fragments in enumerate(pools):
        pool_diversity = calculate_diversity_score(pool_fragments, similarity_matrix)
        
        for frag in pool_fragments:
            pool_list.append({
                'PoolID': f'Pool{pool_idx + 1}', 
                'FragmentID': df.iloc[frag]['FragmentID'],
                'PoolDiversity': pool_diversity
            })
            total_fragments_assigned += 1
    
    # Save pool assignment file
    output_pools_file = f"{args.output_prefix}_pools.csv"
    df_pools = pd.DataFrame(pool_list)
    df_pools.to_csv(output_pools_file, index=False)
    
    # ---------- Step 8: Report results ----------
    if args.verbose:
        print(f"\n=== DIVERSITY OPTIMIZATION RESULTS ===")
        print(f"Method: {args.diversity_method}")
        print(f"Created {len(pools)} pools")
        
        pool_sizes = [len(pool) for pool in pools]
        print(f"Pool sizes: min={min(pool_sizes)}, max={max(pool_sizes)}, avg={sum(pool_sizes)/len(pool_sizes):.1f}")
        
        # Calculate overall diversity metrics
        total_diversity = 0
        for pool in pools:
            diversity = calculate_diversity_score(pool, similarity_matrix)
            total_diversity += diversity * len(pool)
        
        avg_diversity = total_diversity / total_fragments_assigned if total_fragments_assigned > 0 else 0
        print(f"Average pool diversity: {avg_diversity:.3f}")
        print(f"Total fragments assigned: {total_fragments_assigned}")
        
        # Check for similar compounds in same pool
        similar_in_same_pool = 0
        for pool in pools:
            for i, j in combinations(pool, 2):
                if similarity_matrix[i][j] > args.cutoff:
                    similar_in_same_pool += 1
        
        if similar_in_same_pool == 0:
            print("✅ No similar compounds found in same pool!")
        else:
            print(f"⚠️  Warning: {similar_in_same_pool} similar compound pairs in same pool")
        
        print(f"Saved pool assignments to: {output_pools_file}")
    
    print("Advanced diversity optimization complete!")
    if args.props.lower() != 'none':
        print(f"Enriched file saved as: {output_props_file}")
    print(f"Pools file saved as: {output_pools_file}")

if __name__ == "__main__":
    main() 