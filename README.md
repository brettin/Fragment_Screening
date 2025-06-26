# Fragment Clustering and Pooling Tool

A Python script for clustering molecular fragments based on structural similarity and organizing them into screening pools.

## Overview

This tool processes a CSV file containing molecular fragments (identified by SMILES strings) and:

1. **Calculates molecular properties** (molecular weight, lipophilicity, hydrogen bonding, etc.)
2. **Generates molecular fingerprints** for similarity calculations
3. **Clusters fragments** using the Butina algorithm based on structural similarity
4. **Creates screening pools** by grouping fragments into manageable sets

## Input Format

The input CSV file should contain:
- `FragmentID`: Unique identifier for each fragment
- `SMILES`: SMILES string representation of the molecular structure

Example:
```csv
FragmentID,SMILES
F001,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
F002,CC1=CC=C(C=C1)C2=CC(=CC=C2)CC(C)C
F003,CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O
```

## Installation

### Prerequisites
- Python 3.7+
- RDKit
- pandas
- numpy

### Install Dependencies
```bash
pip install rdkit pandas numpy
```

## Usage

### Basic Usage
```bash
python cluster_fragments.py
```

### Command Line Options

#### Input/Output Options
- `-i, --input-file`: Input CSV file path (default: `fragments.csv`)
- `-o, --output-prefix`: Output file prefix (default: `fragment`)
- `-v, --verbose`: Enable verbose output

#### Fingerprint Parameters
- `-b, --fingerprint-bits`: Morgan fingerprint bit count (default: 1024)
  - Range: 256-4096
  - Higher values provide more detailed molecular representation
  - Recommended: 1024 for fragments, 2048 for detailed analysis
- `-r, --fingerprint-radius`: Morgan fingerprint radius (default: 2)
  - Range: 0-7
  - Radius 0: atom types only
  - Radius 1: atom + immediate neighbors
  - Radius 2: atom + neighbors + neighbors of neighbors (recommended for fragments)
  - Radius 3+: larger molecular environment

#### Clustering Parameters
- `-c, --cutoff`: Clustering similarity cutoff (default: 0.4)
  - Range: 0.0-1.0
  - Lower values create tighter clusters
  - Higher values create looser clusters
  - Recommended: 0.4 for fragment diversity, 0.3 for high-throughput

#### Pooling Parameters
- `-p, --pool-size`: Number of fragments per pool (default: 3)
- `-s, --random-seed`: Random seed for reproducible results

#### Property Calculation
- `--props, --properties`: Molecular properties to calculate
  - Options: `MW`, `logP`, `HBD`, `HBA`, `TPSA`, `RotBonds`, `Solubility`
  - Use comma-separated values: `--props MW,logP,HBD`
  - Use `all` for all properties (default)
  - Use `none` to skip property calculation

## Output Files

1. **`{prefix}_with_properties.csv`**: Input data enriched with calculated molecular properties
2. **`{prefix}_pools.csv`**: Fragment pool assignments with columns:
   - `PoolID`: Pool identifier (e.g., "Pool1", "Pool2")
   - `FragmentID`: Original fragment identifier

## Molecular Properties

| Property | Description | Units | Relevance |
|----------|-------------|-------|-----------|
| `MW` | Molecular Weight | g/mol | Size filtering |
| `logP` | Octanol-water partition coefficient | - | Lipophilicity |
| `HBD` | Hydrogen Bond Donors | count | Drug-likeness |
| `HBA` | Hydrogen Bond Acceptors | count | Drug-likeness |
| `TPSA` | Topological Polar Surface Area | Å² | Membrane permeability |
| `RotBonds` | Rotatable Bonds | count | Flexibility |
| `Solubility` | Estimated aqueous solubility | mg/mL | Solubility |

## Examples

### Basic Fragment Screening
```bash
python cluster_fragments.py -i my_fragments.csv -o my_results
```

### High-Throughput Processing
```bash
python cluster_fragments.py \
  -i fragments.csv \
  -b 512 \
  -r 1 \
  -c 0.3 \
  --props MW,logP \
  -p 5
```

### Detailed Analysis
```bash
python cluster_fragments.py \
  -i fragments.csv \
  -b 2048 \
  -r 3 \
  -c 0.5 \
  --props all \
  -p 2 \
  -s 42
```

### Fragment Diversity Analysis
```bash
python cluster_fragments.py \
  -i fragments.csv \
  -c 0.4 \
  --props MW,logP,HBD,HBA,TPSA \
  -p 3
```

## Performance Considerations

| Parameter | Speed Impact | Memory Impact | Quality Impact |
|-----------|-------------|---------------|----------------|
| `-b` (bits) | High | High | Medium |
| `-r` (radius) | Medium | Low | High |
| `--props` | High | Low | High |
| `-c` (cutoff) | Low | Low | High |

## Troubleshooting

### Common Issues

1. **"Bad SMILES" warnings**: Some SMILES strings may be invalid or unsupported
2. **Memory errors**: Reduce fingerprint bits (`-b 512`) or process smaller datasets
3. **Slow processing**: Use `-r 1` and `--props MW,logP` for faster processing
4. **Too many/few clusters**: Adjust cutoff value (`-c`)

### Performance Tips

- For large datasets (>10,000 fragments): Use `-b 512 -r 1`
- For fragment screening: Use `--props MW,logP,HBD,HBA`
- For reproducible results: Always set `-s` (random seed)

## Algorithm Details

1. **Property Calculation**: Uses RDKit descriptors and ESOL-like solubility estimation
2. **Fingerprinting**: Morgan circular fingerprints for structural similarity
3. **Clustering**: Butina algorithm with Tanimoto similarity
4. **Pooling**: Random assignment within clusters to maximize diversity

## License

This tool is provided as-is for research and development purposes. 