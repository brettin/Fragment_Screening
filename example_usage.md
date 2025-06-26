# Example Usage Guide

This guide demonstrates how to use the `cluster_fragments.py` script with various parameter combinations.

## Prerequisites

First, install the required dependencies:
```bash
pip install rdkit pandas numpy
```

## Test Data

The `test_fragments.csv` file contains 20 sample molecular fragments with different structural features:
- Fragments with carboxylic acids (F001, F003, F005, F007, F009, F011, F013, F015, F017, F019)
- Fragments with amines (F005, F006, F007, F008)
- Fragments with halogens (F009, F010, F011, F012, F017, F018, F019, F020)
- Fragments with sulfur (F013, F014, F015, F016)

## Example Commands

### 1. Basic Usage
```bash
python cluster_fragments.py -i test_fragments.csv -v
```
**What it does:**
- Uses default parameters (1024-bit fingerprints, radius 2, cutoff 0.4)
- Calculates all molecular properties
- Creates pools of 3 fragments each
- Outputs: `fragment_with_properties.csv` and `fragment_pools.csv`

### 2. High-Throughput Processing
```bash
python cluster_fragments.py \
  -i test_fragments.csv \
  -b 512 \
  -r 1 \
  -c 0.3 \
  --props MW,logP \
  -p 5 \
  -v
```
**What it does:**
- Uses smaller fingerprints (512 bits, radius 1) for faster processing
- Looser clustering (cutoff 0.3) for more diverse pools
- Calculates only molecular weight and lipophilicity
- Creates larger pools (5 fragments each)
- **Use case:** Large fragment libraries where speed is important

### 3. Detailed Analysis
```bash
python cluster_fragments.py \
  -i test_fragments.csv \
  -b 2048 \
  -r 3 \
  -c 0.5 \
  --props all \
  -p 2 \
  -s 42 \
  -v
```
**What it does:**
- Uses high-resolution fingerprints (2048 bits, radius 3)
- Tight clustering (cutoff 0.5) for similar fragments
- Calculates all properties including solubility
- Creates small pools (2 fragments each)
- Uses fixed random seed (42) for reproducible results
- **Use case:** Detailed analysis of fragment diversity

### 4. Fragment Diversity Analysis
```bash
python cluster_fragments.py \
  -i test_fragments.csv \
  -c 0.4 \
  --props MW,logP,HBD,HBA,TPSA \
  -p 3 \
  -v
```
**What it does:**
- Uses standard clustering parameters
- Calculates key drug-likeness properties
- **Use case:** Fragment screening for drug discovery

### 5. Clustering Only (No Properties)
```bash
python cluster_fragments.py \
  -i test_fragments.csv \
  --props none \
  -p 4 \
  -v
```
**What it does:**
- Skips property calculation entirely
- Only performs clustering and pooling
- **Use case:** When you only need structural clustering

### 6. Custom Output Files
```bash
python cluster_fragments.py \
  -i test_fragments.csv \
  -o my_results \
  -v
```
**What it does:**
- Creates output files: `my_results_with_properties.csv` and `my_results_pools.csv`
- **Use case:** When you want to keep multiple analysis results

## Expected Output

### Properties File (`*_with_properties.csv`)
Contains original data plus calculated properties:
```csv
FragmentID,SMILES,MW,logP,HBD,HBA,TPSA,RotBonds,Solubility
F001,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O,206.28,3.12,1,3,37.3,4,0.045
F002,CC1=CC=C(C=C1)C2=CC(=CC=C2)CC(C)C,210.36,4.89,0,0,0.0,3,0.002
...
```

### Pools File (`*_pools.csv`)
Contains fragment pool assignments:
```csv
PoolID,FragmentID
Pool1,F001
Pool1,F005
Pool1,F009
Pool2,F002
Pool2,F006
Pool2,F010
...
```

## Parameter Impact Examples

### Fingerprint Bits (-b)
- `-b 256`: Fast, coarse molecular representation
- `-b 1024`: Balanced (default)
- `-b 2048`: Detailed, slower processing

### Fingerprint Radius (-r)
- `-r 1`: Atom + immediate neighbors
- `-r 2`: Atom + neighbors + neighbors of neighbors (default)
- `-r 3`: Larger molecular environment

### Clustering Cutoff (-c)
- `-c 0.3`: Loose clustering, fewer clusters, more diverse pools
- `-c 0.4`: Balanced clustering (default)
- `-c 0.5`: Tight clustering, more clusters, similar fragments together

### Pool Size (-p)
- `-p 2`: Small pools, more pools
- `-p 3`: Standard pools (default)
- `-p 5`: Large pools, fewer pools

## Performance Tips

1. **For large datasets (>10,000 fragments):**
   ```bash
   python cluster_fragments.py -b 512 -r 1 --props MW,logP
   ```

2. **For reproducible results:**
   ```bash
   python cluster_fragments.py -s 42
   ```

3. **For fragment screening:**
   ```bash
   python cluster_fragments.py --props MW,logP,HBD,HBA,TPSA
   ```

4. **For detailed analysis:**
   ```bash
   python cluster_fragments.py -b 2048 -r 3 --props all
   ``` 