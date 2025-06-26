#!/usr/bin/env python3
"""
Test script for cluster_fragments.py
Demonstrates various usage scenarios with different parameters.
"""

import subprocess
import sys
import os

def run_test(command, description):
    """Run a test command and print results."""
    print(f"\n{'='*60}")
    print(f"TEST: {description}")
    print(f"COMMAND: {command}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        print("STDOUT:")
        print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        print(f"Return code: {result.returncode}")
        return result.returncode == 0
    except Exception as e:
        print(f"Error running command: {e}")
        return False

def main():
    """Run various test scenarios."""
    print("Fragment Clustering Tool - Test Examples")
    print("This script demonstrates different usage scenarios.")
    
    # Check if test data exists
    if not os.path.exists('test_fragments.csv'):
        print("Error: test_fragments.csv not found. Please create it first.")
        return
    
    tests = [
        # Basic usage
        ("python cluster_fragments.py -i test_fragments.csv -v", 
         "Basic usage with verbose output"),
        
        # High-throughput processing
        ("python cluster_fragments.py -i test_fragments.csv -b 512 -r 1 -c 0.3 --props MW,logP -p 5 -v", 
         "High-throughput processing with minimal properties"),
        
        # Detailed analysis
        ("python cluster_fragments.py -i test_fragments.csv -b 2048 -r 3 -c 0.5 --props all -p 2 -s 42 -v", 
         "Detailed analysis with all properties and reproducible results"),
        
        # Fragment diversity analysis
        ("python cluster_fragments.py -i test_fragments.csv -c 0.4 --props MW,logP,HBD,HBA,TPSA -p 3 -v", 
         "Fragment diversity analysis with key properties"),
        
        # No property calculation
        ("python cluster_fragments.py -i test_fragments.csv --props none -p 4 -v", 
         "Clustering only (no property calculation)"),
        
        # Custom output prefix
        ("python cluster_fragments.py -i test_fragments.csv -o test_results -v", 
         "Custom output prefix"),
    ]
    
    successful_tests = 0
    total_tests = len(tests)
    
    for command, description in tests:
        if run_test(command, description):
            successful_tests += 1
        print("\n" + "-"*60)
    
    print(f"\n{'='*60}")
    print(f"TEST SUMMARY: {successful_tests}/{total_tests} tests passed")
    print(f"{'='*60}")
    
    if successful_tests == total_tests:
        print("All tests passed! The clustering tool is working correctly.")
    else:
        print("Some tests failed. Please check the output above for errors.")
    
    # Show generated files
    print("\nGenerated files:")
    for file in os.listdir('.'):
        if file.endswith('.csv') and ('test' in file or 'fragment' in file):
            print(f"  - {file}")

if __name__ == "__main__":
    main() 