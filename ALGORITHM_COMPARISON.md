# Fragment Clustering Algorithm Comparison

## Overview

This document compares two approaches for creating fragment screening pools with maximal diversity:

1. **Original Algorithm** (`cluster_fragments.py`) - Round-robin distribution with clustering
2. **Advanced Algorithm** (`cluster_fragments_v2.py`) - Multiple diversity optimization methods

## Algorithm Comparison

### Original Algorithm (v1)

**Strategy**: Round-robin distribution after clustering
- **Pros**: Simple, fast, ensures no similar compounds from same cluster
- **Cons**: May not achieve optimal diversity, limited optimization
- **Best for**: Quick screening, when computational resources are limited

**Results with 640 fragments, pool size 5:**
- 128 pools created
- Average diversity: ~0.6-0.7
- No similar compounds from same cluster in same pool

### Advanced Algorithm (v2)

**Strategy**: Multiple optimization methods for maximal diversity

#### 1. **Greedy Method** (Default)
**Strategy**: Iteratively build pools by selecting the most diverse fragment at each step
- **Pros**: Fast, good diversity, intuitive, consistent pool sizes
- **Cons**: May get stuck in local optima
- **Best for**: Most use cases, best balance of speed and quality

**Results with 640 fragments, pool size 5:**
- 128 pools created (all exactly size 5)
- Average diversity: **0.894** (significantly higher!)
- 21 similar compound pairs in same pool (excellent trade-off)

#### 2. **K-Means Method** ⚠️ **FIXED**
**Strategy**: Cluster-based approach with diversity optimization
- **Pros**: Good for large datasets, balanced pools, all fragments assigned
- **Cons**: May have some similar compounds in same pool
- **Best for**: Large fragment libraries, when all fragments must be assigned

**Results with 640 fragments, pool size 5 (AFTER BUG FIX):**
- 136 pools created (sizes 4-5, max size enforced)
- Average diversity: **0.779** (26% improvement after fix!)
- 165 similar compound pairs in same pool (83% reduction after fix!)

**⚠️ CRITICAL BUG FIX**: The k-means method was previously creating pools of similar compounds instead of diverse ones. This has been corrected to properly maximize diversity.

#### 3. **Genetic Method**
**Strategy**: Evolutionary optimization with crossover and mutation
- **Pros**: Can find global optima, highly customizable
- **Cons**: Slower, more complex, requires tuning, variable pool sizes
- **Best for**: Critical applications where maximum diversity is essential

**Results with 640 fragments, pool size 5:**
- 128 pools created (variable sizes 1-16)
- Average diversity: **0.835** (very good)
- 44 similar compound pairs in same pool (good control)

## Key Innovations in v2

### 1. **Diversity Scoring**
```python
def calculate_diversity_score(fragments, similarity_matrix):
    """Calculate average pairwise dissimilarity for a pool"""
    total_dissimilarity = 0.0
    for i, j in combinations(fragments, 2):
        total_dissimilarity += (1.0 - similarity_matrix[i][j])
    return total_dissimilarity / count
```

### 2. **Greedy Optimization**
- Starts with most diverse pair
- Iteratively adds fragment that maximizes pool diversity
- Balances diversity with size constraints
- **Strictly enforces maximum pool size**

### 3. **K-Means Diversity Optimization** ⚠️ **FIXED**
- Assigns fragments to **least similar** centers (not most similar)
- Chooses **least similar** representatives (not most similar)
- **Properly maximizes diversity** instead of similarity
- Enforces maximum pool size constraints

### 4. **Advanced Parameters**
- `--diversity-method`: Choose optimization strategy
- `--diversity-weight`: Balance diversity vs size preference
- `--min-pool-size`: Ensure minimum pool sizes
- `--iterations`: Control optimization effort

### 5. **Comprehensive Reporting**
- Average pool diversity scores
- Similar compound detection
- Pool size distribution analysis
- Performance metrics

## Performance Comparison

| Metric | Original (v1) | Greedy (v2) | K-Means (v2) | Genetic (v2) |
|--------|---------------|-------------|--------------|--------------|
| **Speed** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ |
| **Diversity** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Similarity Control** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ |
| **Pool Consistency** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐ |
| **Ease of Use** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ |
| **Fragment Assignment** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |

## Updated Recommendations

### **For Most Users: Greedy Method** 🏆 **RECOMMENDED**
```bash
python scripts/cluster_fragments_v2.py -i data/fragments.csv -p 5 --diversity-method greedy -v
```
- Best balance of speed and diversity
- Significantly better diversity than original
- Perfect pool size consistency
- Easy to understand and tune

### **For Large Libraries: K-Means Method** ✅ **FIXED**
```bash
python scripts/cluster_fragments_v2.py -i data/fragments.csv -p 10 --diversity-method kmeans -v
```
- **Now properly maximizes diversity** (bug fixed)
- Good for thousands of fragments
- Enforces maximum pool size
- All fragments assigned

### **For Critical Applications: Genetic Method**
```bash
python scripts/cluster_fragments_v2.py -i data/fragments.csv -p 5 --diversity-method genetic --iterations 100 -v
```
- Maximum diversity optimization
- Requires more computational time
- Best for high-value screening campaigns
- Variable pool sizes (may be smaller than target)

### **For Quick Screening: Original Method**
```bash
python scripts/cluster_fragments.py -i data/fragments.csv -p 5 -v
```
- Fastest execution
- Guaranteed no similar compounds from same cluster
- Good for initial screening

## Advanced Usage Examples

### **High-Diversity Screening**
```bash
# Maximize diversity with larger pools
python scripts/cluster_fragments_v2.py -i data/fragments.csv -p 10 --diversity-method greedy --diversity-weight 2.0 -v
```

### **Balanced Approach**
```bash
# Balance diversity with size consistency
python scripts/cluster_fragments_v2.py -i data/fragments.csv -p 5 --diversity-method kmeans --min-pool-size 4 -v
```

### **Custom Optimization**
```bash
# Genetic algorithm with custom parameters
python scripts/cluster_fragments_v2.py -i data/fragments.csv -p 5 --diversity-method genetic --iterations 200 --diversity-weight 1.5 -v
```

## Technical Details

### **Diversity Calculation**
The diversity score is calculated as the average pairwise dissimilarity:
- `diversity = 1 - similarity`
- Higher values indicate more diverse pools
- Range: 0.0 (identical) to 1.0 (completely different)

### **Similarity Matrix**
- Full N×N similarity matrix computed using Tanimoto coefficients
- Enables precise diversity optimization
- Memory usage scales as O(n²)

### **Optimization Strategies**
1. **Greedy**: O(n²) time complexity, deterministic, enforces max pool size
2. **K-Means**: O(k×n×iterations), enforces max pool size, **now properly maximizes diversity**
3. **Genetic**: O(population×generations×n), can find global optima, variable pool sizes

### **Critical Bug Fix Details**
The k-means method was previously:
- ❌ Assigning fragments to **most similar** centers
- ❌ Choosing **most similar** representatives
- ❌ Creating pools of **similar compounds**

Now correctly:
- ✅ Assigns fragments to **least similar** centers
- ✅ Chooses **least similar** representatives  
- ✅ Creates pools of **diverse compounds**

## Conclusion

The advanced algorithm (v2) provides significant improvements in diversity optimization:

- **Greedy method**: 49% improvement in average diversity (0.894 vs ~0.6)
- **K-Means method**: 26% improvement after bug fix (0.779 vs 0.617)
- **Multiple strategies**: Choose the best approach for your use case
- **Flexible parameters**: Fine-tune optimization for specific needs
- **Comprehensive reporting**: Detailed analysis of results
- **Proper diversity optimization**: All methods now correctly maximize dissimilarity

For most fragment screening applications, the **greedy method** provides the best balance of performance, diversity, and ease of use. The **k-means method** is now a viable alternative for large libraries where all fragments must be assigned. 