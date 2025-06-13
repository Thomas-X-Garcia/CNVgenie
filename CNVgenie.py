#!/usr/bin/env python3
"""
Integrated CNV Detection Pipeline for T2T-aligned ONT data
Complete workflow: BAM/CRAM → mosdepth → baseline → CNV detection

This comprehensive pipeline integrates:
1. Batch mosdepth processing for coverage analysis
2. Robust baseline generation from multiple samples
3. Classic and optimized CNV detection algorithms
4. Quality control and reporting

Author: Thomas X. Garcia, PhD, HCLD
Version: 1.0
"""

import numpy as np
import pandas as pd
import gzip
import json
import argparse
import logging
import sys
import os
import pickle
import time
import gc
import subprocess
import shutil
import signal
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path
from scipy import stats
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
AUTOSOME_CHROMS = [f'chr{i}' for i in range(1, 23)]
SEX_CHROMS = ['chrX', 'chrY']
MITO_CHROM = 'chrM'
MIN_MAPPABILITY = 0.5  # Minimum mappability score to include region
MIN_SAMPLES_FOR_BASELINE = 20  # Minimum samples required for robust baseline
DEFAULT_BIN_SIZE = 500  # Default bin size from mosdepth
CHUNK_SIZE = 100000  # Process regions in chunks for memory efficiency

# Global utility functions
def validate_file_exists(filepath: str, description: str = "File") -> None:
    """Validate that a file exists and is readable."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{description} not found: {filepath}")
    if not os.access(filepath, os.R_OK):
        raise PermissionError(f"{description} is not readable: {filepath}")

def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """Safely divide two numbers, returning default if denominator is zero."""
    if abs(denominator) < 1e-10:  # Avoid division by very small numbers
        return default
    return numerator / denominator

def setup_logging_with_file(output_dir: str, log_name: str = "pipeline") -> None:
    """Set up logging with both file and console output."""
    log_file = Path(output_dir) / f"{log_name}_{time.strftime('%Y%m%d_%H%M%S')}.log"
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    root_logger.handlers = []  # Clear existing handlers
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)
    
    logger.info(f"Logging to file: {log_file}")

def read_regions_with_line_tracking(regions_file: str, debug: bool = False) -> pd.DataFrame:
    """
    Read mosdepth regions.bed.gz file with comprehensive error handling and line tracking.
    
    This is a centralized utility function used by all classes to ensure consistent
    file reading behavior and error reporting.
    
    Args:
        regions_file: Path to the regions.bed.gz file
        debug: Enable debug logging
        
    Returns:
        DataFrame with columns: chrom, start, end, coverage
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    validate_file_exists(regions_file, "Regions file")
    
    regions = []
    line_num = 0
    
    try:
        with gzip.open(regions_file, 'rt') as f:
            for line in f:
                line_num += 1
                parts = line.strip().split('\t')
                
                if len(parts) < 4:
                    if debug:
                        logger.debug(f"Line {line_num}: Skipping line with {len(parts)} columns")
                    continue
                
                try:
                    regions.append({
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'coverage': float(parts[3])
                    })
                except ValueError as e:
                    logger.warning(f"Line {line_num} in {regions_file}: {e}")
                    if debug:
                        logger.debug(f"Problem line: {line.strip()}")
                    continue
                    
    except Exception as e:
        logger.error(f"Error reading {regions_file} at line {line_num}: {e}")
        raise
        
    if debug:
        logger.debug(f"Successfully read {len(regions)} regions from {regions_file}")
        
    if len(regions) == 0:
        raise ValueError(f"No valid regions found in {regions_file}")
        
    return pd.DataFrame(regions)

@dataclass
class QCMetrics:
    """Quality control metrics for a sample with enhanced documentation."""
    sample_name: str
    total_reads: int
    mean_coverage: float
    median_coverage: float
    coverage_std: float
    coverage_mad: float
    gc_bias_score: float  # Coefficient of variation (CV), not actual GC bias
    mappability_score: float
    sex: str  # 'XX', 'XY', or 'unknown'
    outlier_bins: int
    quality_score: float  # Overall quality score 0-1

class MosdepthProcessor:
    """
    Integrated mosdepth processing functionality for batch BAM/CRAM analysis.
    
    Provides comprehensive coverage analysis capabilities with support for
    both BAM and CRAM files, parallel processing, and robust error handling.
    """
    
    def __init__(self, output_dir: str):
        """
        Initialize MosdepthProcessor.
        
        Args:
            output_dir: Directory for mosdepth output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def validate_alignment_file(self, alignment_path: str) -> Tuple[bool, str]:
        """Validate that alignment file (BAM/CRAM) exists and has an index."""
        alignment_path = Path(alignment_path)
        
        if not alignment_path.exists():
            return False, f"Alignment file does not exist: {alignment_path}"
        
        if not alignment_path.is_file():
            return False, f"Path is not a file: {alignment_path}"
        
        # Check file extension
        suffix = alignment_path.suffix.lower()
        if suffix not in ['.bam', '.cram']:
            return False, f"File is not a BAM or CRAM file: {alignment_path}"
        
        # Check for index files
        if suffix == '.bam':
            # Check for BAM index (.bai or .bam.bai)
            bai_path1 = alignment_path.with_suffix('.bam.bai')
            bai_path2 = alignment_path.with_suffix('.bai')
            
            if not (bai_path1.exists() or bai_path2.exists()):
                return False, f"BAM index not found for: {alignment_path}"
        
        elif suffix == '.cram':
            # Check for CRAM index (.crai or .cram.crai)
            crai_path1 = alignment_path.with_suffix('.cram.crai')
            crai_path2 = alignment_path.with_suffix('.crai')
            
            if not (crai_path1.exists() or crai_path2.exists()):
                return False, f"CRAM index not found for: {alignment_path}"
        
        return True, "OK"
    
    def check_mosdepth_installation(self) -> Tuple[bool, str]:
        """Check if mosdepth is installed and accessible."""
        try:
            result = subprocess.run(
                ['mosdepth', '--version'],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode == 0:
                version = result.stdout.strip()
                return True, version
            else:
                return False, "mosdepth command failed"
        except FileNotFoundError:
            return False, "mosdepth not found in PATH"
        except Exception as e:
            return False, f"Error checking mosdepth: {str(e)}"
    
    def run_mosdepth(self, alignment_path: str, bin_size: int = 500, threads: int = 4, 
                     fast_mode: bool = True, reference: Optional[str] = None) -> Tuple[bool, str, str]:
        """
        Run mosdepth on a single BAM/CRAM file.
        
        Args:
            alignment_path: Path to BAM/CRAM file
            bin_size: Bin size for coverage calculation
            threads: Number of threads to use
            fast_mode: Enable fast mode
            reference: Reference genome for CRAM files
            
        Returns:
            Tuple of (success, message, output_prefix)
        """
        alignment_path = Path(alignment_path)
        
        # Create output prefix based on alignment filename
        output_prefix = self.output_dir / alignment_path.stem
        
        # Build mosdepth command
        cmd = [
            'mosdepth',
            '-n',  # don't output per-base depth
            '-t', str(threads),  # number of threads
            '--by', str(bin_size),  # bin size
        ]
        
        # Add fast mode flag if enabled
        if fast_mode:
            cmd.append('--fast-mode')
        
        # Add reference for CRAM files
        if alignment_path.suffix.lower() == '.cram' and reference:
            cmd.extend(['--fasta', str(reference)])
        
        cmd.extend([
            str(output_prefix),  # output prefix
            str(alignment_path)  # input BAM/CRAM
        ])
        
        try:
            # Run mosdepth
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Check if output files were created
            expected_files = [
                f"{output_prefix}.regions.bed.gz",
                f"{output_prefix}.regions.bed.gz.csi",
                f"{output_prefix}.mosdepth.global.dist.txt",
                f"{output_prefix}.mosdepth.region.dist.txt",
                f"{output_prefix}.mosdepth.summary.txt"
            ]
            
            missing_files = [f for f in expected_files if not Path(f).exists()]
            
            if missing_files:
                return False, f"Some output files were not created: {missing_files}", str(output_prefix)
            
            return True, "Success", str(output_prefix)
            
        except subprocess.CalledProcessError as e:
            error_msg = f"mosdepth failed with exit code {e.returncode}\n"
            error_msg += f"STDOUT: {e.stdout}\n"
            error_msg += f"STDERR: {e.stderr}"
            return False, error_msg, str(output_prefix)
        
        except Exception as e:
            return False, f"Unexpected error: {str(e)}", str(output_prefix)
    
    def is_alignment_file(self, filepath: str) -> bool:
        """Check if a file is a BAM or CRAM file."""
        try:
            path = Path(filepath)
            return path.is_file() and path.suffix.lower() in ['.bam', '.cram']
        except:
            return False
    
    def has_cram_files(self, alignment_files: List[str]) -> bool:
        """Check if any of the alignment files are CRAM files."""
        return any(Path(f).suffix.lower() == '.cram' for f in alignment_files)
    
    def process_alignment_files(self, alignment_files: List[str], bin_size: int = 500, 
                               threads: int = 4, max_workers: int = 1, fast_mode: bool = True, 
                               reference: Optional[str] = None) -> int:
        """
        Process all alignment files in the list.
        
        Args:
            alignment_files: List of alignment file paths
            bin_size: Bin size for coverage calculation
            threads: Threads per mosdepth job
            max_workers: Maximum parallel jobs
            fast_mode: Enable fast mode
            reference: Reference genome for CRAM files
            
        Returns:
            Exit code (0 for success, 1 for failure)
        """
        logger.info(f"Starting batch mosdepth processing")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Bin size: {bin_size} bp")
        logger.info(f"Threads per mosdepth job: {threads}")
        logger.info(f"Maximum parallel jobs: {max_workers}")
        logger.info(f"Fast mode: {'enabled' if fast_mode else 'disabled'}")
        if reference:
            logger.info(f"Reference genome: {reference}")
        logger.info(f"Found {len(alignment_files)} alignment files to process")
        
        # Check if CRAM files are present but no reference provided
        if self.has_cram_files(alignment_files) and not reference:
            cram_files = [f for f in alignment_files if Path(f).suffix.lower() == '.cram']
            logger.warning(f"Found {len(cram_files)} CRAM files but no reference genome provided.")
            logger.warning("CRAM files require a reference genome. Use --reference option or set REF_PATH environment variable.")
        
        # Validate all alignment files first
        valid_files = []
        for alignment_file in alignment_files:
            is_valid, message = self.validate_alignment_file(alignment_file)
            if is_valid:
                # Additional check for CRAM files without reference
                if Path(alignment_file).suffix.lower() == '.cram' and not reference and not os.environ.get('REF_PATH'):
                    logger.warning(f"Skipping CRAM file (no reference provided): {alignment_file}")
                else:
                    valid_files.append(alignment_file)
                    logger.info(f"Validated: {alignment_file}")
            else:
                logger.warning(f"Skipping invalid alignment file: {message}")
        
        if not valid_files:
            logger.error("No valid alignment files to process")
            return 1
        
        logger.info(f"Processing {len(valid_files)} valid alignment files")
        
        # Process alignment files
        successful = 0
        failed = 0
        
        if max_workers == 1:
            # Sequential processing
            for i, alignment_file in enumerate(valid_files, 1):
                logger.info(f"Processing [{i}/{len(valid_files)}]: {alignment_file}")
                success, message, output_prefix = self.run_mosdepth(
                    alignment_file, bin_size, threads, fast_mode, reference
                )
                
                if success:
                    successful += 1
                    logger.info(f"Successfully processed: {alignment_file}")
                    logger.info(f"Output prefix: {output_prefix}")
                else:
                    failed += 1
                    logger.error(f"Failed to process {alignment_file}: {message}")
        else:
            # Parallel processing
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit all jobs
                future_to_file = {
                    executor.submit(
                        self.run_mosdepth, alignment_file, bin_size, threads, fast_mode, reference
                    ): alignment_file
                    for alignment_file in valid_files
                }
                
                # Process completed jobs
                for i, future in enumerate(as_completed(future_to_file), 1):
                    alignment_file = future_to_file[future]
                    try:
                        success, message, output_prefix = future.result()
                        
                        if success:
                            successful += 1
                            logger.info(f"[{i}/{len(valid_files)}] Successfully processed: {alignment_file}")
                            logger.info(f"Output prefix: {output_prefix}")
                        else:
                            failed += 1
                            logger.error(f"[{i}/{len(valid_files)}] Failed to process {alignment_file}: {message}")
                    
                    except Exception as e:
                        failed += 1
                        logger.error(f"[{i}/{len(valid_files)}] Exception processing {alignment_file}: {str(e)}")
        
        # Summary
        logger.info("=" * 60)
        logger.info(f"Mosdepth processing complete:")
        logger.info(f"  Total alignment files: {len(alignment_files)}")
        logger.info(f"  Valid alignment files: {len(valid_files)}")
        logger.info(f"  Successfully processed: {successful}")
        logger.info(f"  Failed: {failed}")
        
        return 0 if failed == 0 else 1
    
    def process_single_or_batch(self, alignment_input: str, bin_size: int = 500, 
                               threads: int = 4, max_workers: int = 1, fast_mode: bool = True, 
                               reference: Optional[str] = None) -> int:
        """
        Process either a single alignment file or a batch from a list file.
        
        Args:
            alignment_input: Single BAM/CRAM file or text file with file list
            bin_size: Bin size for coverage calculation
            threads: Threads per mosdepth job
            max_workers: Maximum parallel jobs
            fast_mode: Enable fast mode
            reference: Reference genome for CRAM files
            
        Returns:
            Exit code (0 for success, 1 for failure)
        """
        # Determine if input is a single alignment file or a list file
        alignment_files = []
        
        if self.is_alignment_file(alignment_input):
            # Single alignment file
            alignment_files = [alignment_input]
            logger.info(f"Processing single alignment file: {alignment_input}")
            
            # Check if it's a CRAM file without reference
            if Path(alignment_input).suffix.lower() == '.cram' and not reference and not os.environ.get('REF_PATH'):
                logger.error("CRAM files require a reference genome. Use --reference option or set REF_PATH environment variable.")
                return 1
        else:
            # Assume it's a text file with list of alignment files
            if not os.path.isfile(alignment_input):
                logger.error(f"Input file does not exist: {alignment_input}")
                return 1
            
            try:
                with open(alignment_input, 'r') as f:
                    alignment_files = [line.strip() for line in f if line.strip()]
                logger.info(f"Reading alignment file list from: {alignment_input}")
            except Exception as e:
                logger.error(f"Failed to read alignment list file: {str(e)}")
                return 1
        
        if not alignment_files:
            logger.error("No alignment files found to process")
            return 1
        
        return self.process_alignment_files(
            alignment_files, bin_size, threads, max_workers, fast_mode, reference
        )

class FlexibleBaselineParser:
    """
    Mixin class for flexible baseline parsing with enhanced error handling.
    
    Provides robust parsing capabilities for different baseline file formats
    with comprehensive validation and error reporting.
    """
    
    def _read_baseline_flexible(self, baseline_file: str) -> Tuple[pd.DataFrame, Dict]:
        """
        Read baseline with flexible format detection and comprehensive validation.
        
        Args:
            baseline_file: Path to baseline TSV file
            
        Returns:
            Tuple of (baseline_df, format_info)
            
        Raises:
            FileNotFoundError: If baseline file doesn't exist
            ValueError: If required columns are missing or format is invalid
        """
        validate_file_exists(baseline_file, "Baseline file")
        
        format_info = {}
        
        # First pass: analyze header
        with open(baseline_file, 'r') as f:
            header_line = f.readline()
            if not header_line.strip():
                raise ValueError(f"Empty baseline file: {baseline_file}")
                
            header_parts = header_line.strip().split('\t')
            
            logger.info(f"Baseline file header ({len(header_parts)} columns): {header_parts}")
            
            # Detect format
            format_info['n_columns'] = len(header_parts)
            format_info['has_region_column'] = 'region' in header_parts
            format_info['header'] = header_parts
            
            # Sample first few lines for format detection
            sample_lines = []
            for i in range(5):
                line = f.readline()
                if not line:
                    break
                sample_lines.append(line.strip().split('\t'))
            
            # Detect column positions
            format_info['column_map'] = self._detect_column_positions(header_parts, sample_lines)
        
        # Second pass: read with pandas
        try:
            baseline_df = pd.read_csv(baseline_file, sep='\t')
        except Exception as e:
            raise ValueError(f"Failed to read baseline file {baseline_file}: {e}")
        
        if baseline_df.empty:
            raise ValueError(f"Baseline file is empty: {baseline_file}")
        
        # Validate and remap if needed
        required = ['chrom', 'start', 'end', 'median_normalized', 'mad']
        if not all(col in baseline_df.columns for col in required):
            baseline_df = self._remap_columns(baseline_df, format_info)
        
        # Handle n_samples column variations
        n_samples_col = self._find_n_samples_column(baseline_df, format_info)
        if n_samples_col and n_samples_col != 'n_samples':
            baseline_df['n_samples'] = baseline_df[n_samples_col]
            
        # Final validation
        required_columns = ['chrom', 'start', 'end', 'median_normalized', 'mad']
        missing_columns = [col for col in required_columns if col not in baseline_df.columns]

        if missing_columns:
            available_columns = list(baseline_df.columns)
            error_msg = (f"Required columns missing from baseline file: {missing_columns}. "
                        f"Available columns: {available_columns}. "
                        f"Format detection found: {format_info.get('column_map', {})}")
            logger.error(error_msg)
            raise ValueError(error_msg)

        # Validate n_samples column
        if 'n_samples' not in baseline_df.columns:
            logger.warning("n_samples column not found, defaulting to 1 for all regions")
            baseline_df['n_samples'] = 1
        
        # Validate data types and ranges
        self._validate_baseline_data(baseline_df)
        
        return baseline_df, format_info
    
    def _validate_baseline_data(self, baseline_df: pd.DataFrame) -> None:
        """Validate baseline data for consistency and valid ranges."""
        # Check for negative values where they shouldn't exist
        if (baseline_df['median_normalized'] < 0).any():
            logger.warning("Found negative median_normalized values in baseline")
        
        if (baseline_df['mad'] < 0).any():
            logger.warning("Found negative MAD values in baseline")
            
        # Check for extreme outliers
        median_range = baseline_df['median_normalized'].quantile([0.01, 0.99])
        if median_range[0.99] / median_range[0.01] > 100:
            logger.warning("Extreme range in median_normalized values detected")
    
    def _detect_column_positions(self, header: List[str], sample_lines: List[List[str]]) -> Dict:
        """Intelligently detect column positions with enhanced pattern matching."""
        column_map = {}
        
        # Enhanced variation patterns
        variations = {
            'chrom': ['chrom', 'chr', 'chromosome', '#chrom', 'seqname'],
            'start': ['start', 'pos', 'begin', 'chromStart', 'position'],
            'end': ['end', 'stop', 'chromEnd', 'finish'],
            'median_normalized': ['median_normalized', 'median_norm', 'median', 'med_norm', 'normalized'],
            'mad': ['mad', 'MAD', 'std', 'stdev', 'deviation'],
            'n_samples': ['n_samples', 'num_samples', 'samples', 'n', 'count', 'sample_count']
        }
        
        # Find by name with case-insensitive matching
        for i, col in enumerate(header):
            col_lower = col.lower().strip()
            for expected, variants in variations.items():
                if any(v.lower() in col_lower for v in variants):
                    column_map[expected] = i
                    break
        
        # Find n_samples by data type if not found by name
        if 'n_samples' not in column_map and sample_lines:
            for i in range(len(header)):
                if i < min(len(line) for line in sample_lines if line):
                    try:
                        values = [int(line[i]) for line in sample_lines if i < len(line) and line[i].strip()]
                        if values and all(v > 0 for v in values):
                            column_map['n_samples'] = i
                            logger.info(f"Detected n_samples in column {i} by data type")
                            break
                    except (ValueError, IndexError):
                        continue
        
        return column_map
    
    def _find_n_samples_column(self, df: pd.DataFrame, format_info: Dict) -> Optional[str]:
        """Find the n_samples column in the dataframe with enhanced search."""
        # Check common names with priority order
        for col_name in ['n_samples', 'n_samples_clean', 'num_samples', 'samples', 'count']:
            if col_name in df.columns:
                return col_name
        
        # Check column map from detection
        if 'column_map' in format_info and 'n_samples' in format_info['column_map']:
            col_idx = format_info['column_map']['n_samples']
            if col_idx < len(df.columns):
                return df.columns[col_idx]
        
        return None
    
    def _remap_columns(self, df: pd.DataFrame, format_info: Dict) -> pd.DataFrame:
        """Remap columns to expected names with validation."""
        column_map = format_info.get('column_map', {})
        rename_dict = {}
        
        for expected_name, col_idx in column_map.items():
            if col_idx < len(df.columns):
                actual_name = df.columns[col_idx]
                if actual_name != expected_name:
                    rename_dict[actual_name] = expected_name
        
        if rename_dict:
            logger.info(f"Remapping columns: {rename_dict}")
            df = df.rename(columns=rename_dict)
        
        return df

class GCBiasCorrector:
    """
    Handles GC bias correction with robust statistical methods.
    
    Implements LOESS-like local regression for GC bias correction
    with enhanced robustness for ONT data characteristics.
    """
    
    def __init__(self, gc_bins: int = 101):
        """
        Initialize GC bias corrector.
        
        Args:
            gc_bins: Number of GC content bins for correction
        """
        self.gc_bins = gc_bins
        self.gc_model = None
        
    def fit(self, gc_content: np.ndarray, coverage: np.ndarray, weights: Optional[np.ndarray] = None) -> None:
        """
        Fit LOESS-like model for GC bias correction.
        
        Args:
            gc_content: Array of GC content values (0-1)
            coverage: Array of coverage values
            weights: Optional weights for fitting
        """
        if len(gc_content) != len(coverage):
            raise ValueError("GC content and coverage arrays must have same length")
        
        # Bin GC content
        gc_binned = np.digitize(gc_content, np.linspace(0, 1, self.gc_bins))
        
        # Calculate robust statistics per GC bin
        self.gc_model = {}
        for i in range(1, self.gc_bins + 1):
            mask = gc_binned == i
            if np.sum(mask) > 10:  # Need sufficient data points
                bin_coverage = coverage[mask]
                if weights is not None:
                    bin_weights = weights[mask]
                    # Weighted median calculation
                    sorted_idx = np.argsort(bin_coverage)
                    sorted_coverage = bin_coverage[sorted_idx]
                    sorted_weights = bin_weights[sorted_idx]
                    cumsum = np.cumsum(sorted_weights)
                    median_idx = np.searchsorted(cumsum, cumsum[-1] / 2)
                    median_cov = sorted_coverage[median_idx]
                else:
                    median_cov = np.median(bin_coverage)
                
                mad_cov = np.median(np.abs(bin_coverage - median_cov))
                self.gc_model[i] = {
                    'median': median_cov,
                    'mad': mad_cov,
                    'n': len(bin_coverage)
                }
        
        # Smooth the model using local regression
        self._smooth_gc_model()
        
    def _smooth_gc_model(self) -> None:
        """Apply smoothing to GC correction factors."""
        if not self.gc_model:
            return
            
        bins = sorted(self.gc_model.keys())
        medians = [self.gc_model[b]['median'] for b in bins]
        
        # Simple moving average smoothing
        window = 5
        if len(medians) >= window:
            smoothed = np.convolve(medians, np.ones(window)/window, mode='same')
            
            for i, b in enumerate(bins):
                self.gc_model[b]['smoothed_median'] = smoothed[i]
    
    def correct(self, gc_content: np.ndarray, coverage: np.ndarray) -> np.ndarray:
        """
        Apply GC bias correction.
        
        Args:
            gc_content: Array of GC content values
            coverage: Array of coverage values to correct
            
        Returns:
            Corrected coverage array
        """
        if self.gc_model is None:
            raise ValueError("GC bias model not fitted")
            
        gc_binned = np.digitize(gc_content, np.linspace(0, 1, self.gc_bins))
        corrected = np.zeros_like(coverage, dtype=float)
        
        # Global median for normalization
        valid_medians = [self.gc_model[b].get('smoothed_median', self.gc_model[b]['median']) 
                        for b in self.gc_model.keys()]
        global_median = np.median(valid_medians) if valid_medians else 1.0
        
        for i in range(len(coverage)):
            bin_idx = gc_binned[i]
            if bin_idx in self.gc_model:
                bin_median = self.gc_model[bin_idx].get('smoothed_median', self.gc_model[bin_idx]['median'])
                correction_factor = safe_divide(global_median, bin_median, 1.0)
                corrected[i] = coverage[i] * correction_factor
            else:
                corrected[i] = coverage[i]
                
        return corrected

class BaselineGenerator:
    """
    Creates robust baseline from multiple samples with advanced QC and memory management.
    
    Supports checkpoint/resume functionality and chunked processing for large datasets.
    """
    
    def __init__(self, min_samples: int = MIN_SAMPLES_FOR_BASELINE, checkpoint_dir: Optional[str] = None):
        """
        Initialize baseline generator.
        
        Args:
            min_samples: Minimum number of samples required for robust baseline
            checkpoint_dir: Directory for checkpoint files (enables resume capability)
        """
        self.min_samples = min_samples
        self.baseline = None
        self.sample_metrics = {}
        self.gc_corrector = GCBiasCorrector()
        self.excluded_regions = set()
        self.checkpoint_dir = checkpoint_dir
        self.processed_samples = set()
        
        # Create checkpoint directory if specified
        if self.checkpoint_dir:
            os.makedirs(self.checkpoint_dir, exist_ok=True)
            self._load_checkpoint()
        
    def _load_checkpoint(self) -> None:
        """Load previous progress from checkpoint files."""
        if not self.checkpoint_dir:
            return
            
        checkpoint_file = os.path.join(self.checkpoint_dir, 'baseline_checkpoint.pkl')
        metrics_file = os.path.join(self.checkpoint_dir, 'sample_metrics.pkl')
        processed_file = os.path.join(self.checkpoint_dir, 'processed_samples.txt')
        
        try:
            if os.path.exists(checkpoint_file):
                with open(checkpoint_file, 'rb') as f:
                    self.baseline = pickle.load(f)
                logger.info(f"Loaded baseline checkpoint with {len(self.baseline) if self.baseline else 0} regions")
                
            if os.path.exists(metrics_file):
                with open(metrics_file, 'rb') as f:
                    self.sample_metrics = pickle.load(f)
                logger.info(f"Loaded metrics for {len(self.sample_metrics)} samples")
                
            if os.path.exists(processed_file):
                with open(processed_file, 'r') as f:
                    self.processed_samples = set(line.strip() for line in f if line.strip())
                logger.info(f"Found {len(self.processed_samples)} previously processed samples")
                    
        except Exception as e:
            logger.warning(f"Failed to load checkpoint: {e}")
            # Reset to clean state on checkpoint load failure
            self.baseline = None
            self.sample_metrics = {}
            self.processed_samples = set()
    
    def _save_checkpoint(self) -> None:
        """Save current progress to checkpoint files."""
        if not self.checkpoint_dir:
            return
            
        try:
            # Save baseline data
            checkpoint_file = os.path.join(self.checkpoint_dir, 'baseline_checkpoint.pkl')
            with open(checkpoint_file, 'wb') as f:
                pickle.dump(self.baseline, f)
                
            # Save sample metrics
            metrics_file = os.path.join(self.checkpoint_dir, 'sample_metrics.pkl')
            with open(metrics_file, 'wb') as f:
                pickle.dump(self.sample_metrics, f)
                
            # Save processed sample list
            processed_file = os.path.join(self.checkpoint_dir, 'processed_samples.txt')
            with open(processed_file, 'w') as f:
                for sample in self.processed_samples:
                    f.write(f"{sample}\n")
                    
        except Exception as e:
            logger.warning(f"Failed to save checkpoint: {e}")
    
    def _clear_checkpoint(self) -> None:
        """Clear checkpoint files after successful completion."""
        if not self.checkpoint_dir:
            return
            
        try:
            for filename in ['baseline_checkpoint.pkl', 'sample_metrics.pkl', 'processed_samples.txt']:
                filepath = os.path.join(self.checkpoint_dir, filename)
                if os.path.exists(filepath):
                    os.remove(filepath)
            logger.info("Checkpoint files cleared")
        except Exception as e:
            logger.warning(f"Failed to clear checkpoint files: {e}")
        
    def add_sample(self, sample_name: str, regions_file: str, summary_file: str, 
                   gc_content: Optional[Dict[str, float]] = None,
                   mappability: Optional[Dict[str, float]] = None) -> bool:
        """
        Add a sample to the baseline with comprehensive QC checks.
        
        Args:
            sample_name: Name of the sample
            regions_file: Path to regions.bed.gz file
            summary_file: Path to summary.txt file
            gc_content: Optional GC content annotation
            mappability: Optional mappability annotation
            
        Returns:
            True if sample was successfully added, False otherwise
        """
        # Skip if already processed (for checkpoint resume)
        if sample_name in self.processed_samples:
            logger.info(f"Skipping already processed sample: {sample_name}")
            return True
            
        try:
            # Validate input files
            validate_file_exists(regions_file, f"Regions file for {sample_name}")
            validate_file_exists(summary_file, f"Summary file for {sample_name}")
            
            # Read genome-wide statistics
            genome_stats = self._read_summary_stats(summary_file)
            
            # Read per-region coverage
            regions_df = read_regions_with_line_tracking(regions_file, debug=logger.isEnabledFor(logging.DEBUG))
            
            # Calculate QC metrics
            qc_metrics = self._calculate_qc_metrics(sample_name, regions_df, genome_stats)
            
            # Check if sample passes QC
            if not self._passes_qc(qc_metrics):
                logger.warning(f"Sample {sample_name} failed QC: quality_score={qc_metrics.quality_score:.3f}")
                self.processed_samples.add(sample_name)  # Mark as processed even if failed
                return False
                
            # Store normalized coverage
            normalized_coverage = self._normalize_coverage(regions_df, genome_stats, gc_content, mappability)
            
            if self.baseline is None:
                self.baseline = defaultdict(list)
                
            for _, row in normalized_coverage.iterrows():
                key = f"{row['chrom']}:{row['start']}-{row['end']}"
                self.baseline[key].append(row['normalized_coverage'])
                
            self.sample_metrics[sample_name] = qc_metrics
            self.processed_samples.add(sample_name)
            
            # Save checkpoint every 10 samples
            if len(self.processed_samples) % 10 == 0:
                self._save_checkpoint()
                logger.info(f"Checkpoint saved after processing {len(self.processed_samples)} samples")
            
            return True
            
        except Exception as e:
            logger.error(f"Error processing sample {sample_name}: {e}")
            return False
    
    def add_samples_chunked(self, sample_files: List[str], chunk_size: int = 50,
                           gc_content: Optional[Dict[str, float]] = None,
                           mappability: Optional[Dict[str, float]] = None) -> None:
        """
        Process samples in chunks for memory efficiency.
        
        Args:
            sample_files: List of summary file paths
            chunk_size: Number of samples to process per chunk
            gc_content: Optional GC content annotation
            mappability: Optional mappability annotation
        """
        total_samples = len(sample_files)
        
        for chunk_start in range(0, total_samples, chunk_size):
            chunk_end = min(chunk_start + chunk_size, total_samples)
            chunk = sample_files[chunk_start:chunk_end]
            
            logger.info(f"Processing chunk {chunk_start//chunk_size + 1} "
                       f"(samples {chunk_start + 1}-{chunk_end} of {total_samples})")
            
            # Process chunk
            chunk_data = defaultdict(list)
            for summary_file in chunk:
                sample_name = os.path.basename(summary_file).replace('.mosdepth.summary.txt', '')
                regions_file = summary_file.replace('.mosdepth.summary.txt', '.regions.bed.gz')
                
                if not os.path.exists(regions_file):
                    logger.warning(f"Regions file not found for {summary_file}, skipping")
                    continue
                
                # Skip if already processed
                if sample_name in self.processed_samples:
                    logger.info(f"Skipping already processed sample: {sample_name}")
                    continue
                
                try:
                    # Read genome-wide statistics
                    genome_stats = self._read_summary_stats(summary_file)
                    
                    # Read per-region coverage
                    regions_df = read_regions_with_line_tracking(regions_file, debug=logger.isEnabledFor(logging.DEBUG))
                    
                    # Calculate QC metrics
                    qc_metrics = self._calculate_qc_metrics(sample_name, regions_df, genome_stats)
                    
                    # Check if sample passes QC
                    if not self._passes_qc(qc_metrics):
                        logger.warning(f"Sample {sample_name} failed QC: quality_score={qc_metrics.quality_score:.3f}")
                        self.processed_samples.add(sample_name)
                        continue
                    
                    # Store normalized coverage
                    normalized_coverage = self._normalize_coverage(regions_df, genome_stats, gc_content, mappability)
                    
                    for _, row in normalized_coverage.iterrows():
                        key = f"{row['chrom']}:{row['start']}-{row['end']}"
                        chunk_data[key].append(row['normalized_coverage'])
                    
                    self.sample_metrics[sample_name] = qc_metrics
                    self.processed_samples.add(sample_name)
                    
                except Exception as e:
                    logger.error(f"Error processing sample {sample_name}: {e}")
                    continue
            
            # Merge chunk into baseline
            if self.baseline is None:
                self.baseline = defaultdict(list)
            
            for region, values in chunk_data.items():
                self.baseline[region].extend(values)
            
            # Clear chunk data to free memory
            chunk_data.clear()
            
            # Force garbage collection for large chunks
            if chunk_end - chunk_start > 100:
                gc.collect()
            
            # Save checkpoint after each chunk
            self._save_checkpoint()
            logger.info(f"Checkpoint saved after processing chunk {chunk_start//chunk_size + 1}")
    
    def _read_summary_stats(self, summary_file: str) -> Dict:
        """Read mosdepth summary statistics with validation."""
        validate_file_exists(summary_file, "Summary file")
        
        stats = {}
        try:
            with open(summary_file, 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        if parts[0] == 'total':
                            stats['mean_coverage'] = float(parts[3])
                            stats['bases'] = int(parts[2])
                            stats['length'] = int(parts[1])
                        elif parts[0] in ['chrX', 'chrY']:
                            stats[f'{parts[0]}_mean'] = float(parts[3])
        except Exception as e:
            logger.error(f"Error reading summary file {summary_file}: {e}")
            raise
            
        if 'mean_coverage' not in stats:
            raise ValueError(f"Invalid summary file format: {summary_file}")
            
        return stats
    
    def _calculate_qc_metrics(self, sample_name: str, regions_df: pd.DataFrame, 
                             genome_stats: Dict) -> QCMetrics:
        """Calculate comprehensive QC metrics with enhanced validation."""
        if regions_df.empty:
            raise ValueError(f"No regions data for sample {sample_name}")
        
        coverage_values = regions_df['coverage'].values
        
        # Basic statistics with safety checks
        mean_cov = np.mean(coverage_values)
        median_cov = np.median(coverage_values)
        std_cov = np.std(coverage_values)
        mad_cov = np.median(np.abs(coverage_values - median_cov))
        
        # Detect sex from X/Y coverage ratio
        sex = 'unknown'
        if 'chrX_mean' in genome_stats and 'chrY_mean' in genome_stats:
            autosome_coverage = genome_stats.get('mean_coverage', 1.0)
            x_ratio = safe_divide(genome_stats['chrX_mean'], autosome_coverage, 0.0)
            y_ratio = safe_divide(genome_stats['chrY_mean'], autosome_coverage, 0.0)
            
            if x_ratio > 0.8 and y_ratio < 0.1:
                sex = 'XX'
            elif x_ratio < 0.6 and y_ratio > 0.3:
                sex = 'XY'
        
        # Coverage variability (coefficient of variation)
        cv_score = safe_divide(std_cov, mean_cov, 1.0)
        
        # Count outlier bins (beyond 3 MAD)
        outlier_threshold = 3 * mad_cov
        outliers = np.sum(np.abs(coverage_values - median_cov) > outlier_threshold)
        
        # Calculate overall quality score
        quality_score = self._calculate_quality_score(
            mean_cov, std_cov, mad_cov, cv_score, safe_divide(outliers, len(coverage_values), 0.0)
        )
        
        return QCMetrics(
            sample_name=sample_name,
            total_reads=int(genome_stats.get('bases', 0) / 150),  # Estimate
            mean_coverage=mean_cov,
            median_coverage=median_cov,
            coverage_std=std_cov,
            coverage_mad=mad_cov,
            gc_bias_score=cv_score,  # Now represents CV, not actual GC bias
            mappability_score=0.95,  # Placeholder
            sex=sex,
            outlier_bins=outliers,
            quality_score=quality_score
        )
    
    def _calculate_quality_score(self, mean_cov: float, std_cov: float, mad_cov: float,
                                cv_score: float, outlier_frac: float) -> float:
        """Calculate overall sample quality score (0-1) optimized for ONT data."""
        scores = []
        
        # Coverage depth score - more lenient for ONT
        if mean_cov >= 25:
            scores.append(1.0)
        elif mean_cov >= 15:
            scores.append(0.9)
        elif mean_cov >= 10:
            scores.append(0.7)
        elif mean_cov >= 7:
            scores.append(0.5)
        else:
            scores.append(0.2)
        
        # Coverage uniformity score - more lenient for ONT's natural variance
        if cv_score <= 0.5:
            uniformity_score = 1.0
        elif cv_score <= 0.8:
            uniformity_score = 0.8
        elif cv_score <= 1.2:
            uniformity_score = 0.6
        elif cv_score <= 1.5:
            uniformity_score = 0.4
        else:
            uniformity_score = 0.2
        scores.append(uniformity_score)
        
        # MAD-based consistency score
        mad_ratio = safe_divide(mad_cov, mean_cov, 1.0)
        if mad_ratio <= 0.3:
            consistency_score = 1.0
        elif mad_ratio <= 0.5:
            consistency_score = 0.8
        elif mad_ratio <= 0.8:
            consistency_score = 0.6
        elif mad_ratio <= 1.0:
            consistency_score = 0.4
        else:
            consistency_score = 0.2
        scores.append(consistency_score)
        
        # Outlier score - slightly more lenient
        if outlier_frac <= 0.05:  # 5% outliers
            outlier_score = 1.0
        elif outlier_frac <= 0.10:  # 10% outliers
            outlier_score = 0.8
        elif outlier_frac <= 0.15:  # 15% outliers
            outlier_score = 0.6
        elif outlier_frac <= 0.20:  # 20% outliers
            outlier_score = 0.4
        else:
            outlier_score = 0.2
        scores.append(outlier_score)
        
        return np.mean(scores)
    
    def _passes_qc(self, metrics: QCMetrics) -> bool:
        """Determine if sample passes QC thresholds - adjusted for ONT data."""
        if metrics.quality_score < 0.5:  # More lenient threshold
            return False
        if metrics.mean_coverage < 7:  # Lower minimum coverage for ONT
            return False
        if metrics.gc_bias_score > 1.5:  # More lenient CV threshold
            return False
        return True
    
    def _normalize_coverage(self, regions_df: pd.DataFrame, genome_stats: Dict,
                          gc_content: Optional[Dict[str, float]] = None,
                          mappability: Optional[Dict[str, float]] = None) -> pd.DataFrame:
        """Normalize coverage with GC and mappability correction."""
        regions_df = regions_df.copy()
        
        # Basic normalization by genome-wide mean
        genome_mean = genome_stats.get('mean_coverage', 1.0)
        if genome_mean <= 0:
            logger.warning(f"Invalid genome mean coverage: {genome_mean}, using 1.0")
            genome_mean = 1.0
            
        regions_df['normalized_coverage'] = regions_df['coverage'] / genome_mean
        
        # GC bias correction if available
        if gc_content is not None:
            gc_values = []
            for _, row in regions_df.iterrows():
                key = f"{row['chrom']}:{row['start']}-{row['end']}"
                gc_values.append(gc_content.get(key, 0.5))
            
            if hasattr(self.gc_corrector, 'gc_model') and self.gc_corrector.gc_model:
                corrected = self.gc_corrector.correct(
                    np.array(gc_values),
                    regions_df['normalized_coverage'].values
                )
                regions_df['normalized_coverage'] = corrected
        
        return regions_df
    
    def finalize_baseline(self, output_file: str) -> pd.DataFrame:
        """Calculate final baseline statistics with outlier removal and enhanced validation."""
        if not self.baseline:
            raise ValueError("No samples in baseline")
            
        n_samples = len(self.sample_metrics)
        if n_samples < self.min_samples:
            raise ValueError(f"Insufficient samples for baseline: {n_samples} < {self.min_samples}")
            
        logger.info(f"Finalizing baseline with {n_samples} samples")
        
        baseline_data = []
        
        for region, values in self.baseline.items():
            if len(values) < 3:
                continue
                
            values_arr = np.array(values)
            
            # Remove outliers using IQR method
            q1, q3 = np.percentile(values_arr, [25, 75])
            iqr = q3 - q1
            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr
            mask = (values_arr >= lower_bound) & (values_arr <= upper_bound)
            clean_values = values_arr[mask]
            
            if len(clean_values) < 3:
                continue
            
            # Calculate robust statistics
            median_cov = np.median(clean_values)
            mad = np.median(np.abs(clean_values - median_cov)) * 1.4826
            
            # Skip regions with very low coverage or no variation
            if median_cov < 0.1 or mad < 0.01:
                self.excluded_regions.add(region)
                continue
            
            # Parse region
            try:
                chrom, pos = region.split(':')
                start, end = pos.split('-')
                start, end = int(start), int(end)
            except ValueError as e:
                logger.warning(f"Invalid region format: {region}, skipping")
                continue
            
            baseline_data.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'median_normalized': median_cov,
                'mad': mad,
                'iqr': iqr,
                'n_samples': len(values),
                'n_samples_clean': len(clean_values),
                'excluded': False
            })
        
        if not baseline_data:
            raise ValueError("No valid baseline regions after filtering")
        
        baseline_df = pd.DataFrame(baseline_data)
        
        # Additional filtering for problematic regions
        self._filter_problematic_regions(baseline_df)
        
        # Save baseline
        baseline_df.to_csv(output_file, sep='\t', index=False)
        
        # Save metadata with proper type conversion
        metadata = self._create_metadata(n_samples)
        metadata_file = output_file.replace('.tsv', '_metadata.json')
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
            
        logger.info(f"Baseline saved to {output_file}")
        logger.info(f"Metadata saved to {metadata_file}")
        
        # Clear checkpoint files on successful completion
        self._clear_checkpoint()
        
        return baseline_df
    
    def _create_metadata(self, n_samples: int) -> Dict:
        """Create metadata dictionary with proper type conversion."""
        def convert_numpy_types(obj):
            """Convert numpy types to standard Python types for JSON serialization."""
            if isinstance(obj, dict):
                return {k: convert_numpy_types(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy_types(v) for v in obj]
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            else:
                return obj
        
        sample_metrics_serializable = {}
        for k, v in self.sample_metrics.items():
            sample_metrics_serializable[k] = convert_numpy_types(v.__dict__)
        
        return {
            'n_samples': int(n_samples),
            'sample_metrics': sample_metrics_serializable,
            'excluded_regions': len(self.excluded_regions),
            'total_regions': len(self.baseline) if self.baseline else 0,
            'sex_distribution': self._get_sex_distribution()
        }
    
    def _filter_problematic_regions(self, baseline_df: pd.DataFrame) -> None:
        """Filter regions with high variance or systematic biases."""
        # Mark regions with very high variance relative to median
        for idx, row in baseline_df.iterrows():
            cv = safe_divide(row['mad'], row['median_normalized'], 0.0)
            if cv > 0.5:  # High coefficient of variation
                baseline_df.at[idx, 'excluded'] = True
    
    def _get_sex_distribution(self) -> Dict[str, int]:
        """Get distribution of sexes in baseline samples."""
        sex_counts = defaultdict(int)
        for metrics in self.sample_metrics.values():
            sex_counts[metrics.sex] += 1
        return dict(sex_counts)

class CNVDetector(FlexibleBaselineParser):
    """
    Classic CNV detector with recursive binary segmentation (original algorithm).
    
    Implements the proven changepoint detection approach with enhanced
    error handling and validation.
    """
    
    def __init__(self, baseline_file: str, metadata_file: Optional[str] = None):
        """
        Initialize CNV detector with baseline data.
        
        Args:
            baseline_file: Path to baseline TSV file
            metadata_file: Optional metadata file path
        """
        self.baseline_df, self.format_info = self._read_baseline_flexible(baseline_file)
        self.baseline_dict = self._create_baseline_dict()
        
        # Load metadata if available
        self.metadata = {}
        if metadata_file and os.path.exists(metadata_file):
            with open(metadata_file, 'r') as f:
                self.metadata = json.load(f)
        elif baseline_file.endswith('.tsv'):
            default_metadata = baseline_file.replace('.tsv', '_metadata.json')
            if os.path.exists(default_metadata):
                with open(default_metadata, 'r') as f:
                    self.metadata = json.load(f)

    def _create_baseline_dict(self) -> Dict:
        """Create dictionary for fast baseline lookup."""
        baseline_dict = {}
        for _, row in self.baseline_df.iterrows():
            if not row.get('excluded', False):
                key = f"{row['chrom']}:{row['start']}-{row['end']}"
                baseline_dict[key] = row
        logger.info(f"Created baseline dictionary with {len(baseline_dict)} regions")
        return baseline_dict
    
    def detect_cnvs(self, summary_file: str, regions_file: str,
                   sample_sex: Optional[str] = None,
                   gc_content: Optional[Dict[str, float]] = None,
                   mappability: Optional[Dict[str, float]] = None) -> pd.DataFrame:
        """
        Detect CNVs in a test sample using classic algorithm.
        
        Args:
            summary_file: Path to sample summary file
            regions_file: Path to sample regions file
            sample_sex: Sample sex ('XX', 'XY', or None for auto-detection)
            gc_content: Optional GC content annotation
            mappability: Optional mappability annotation
            
        Returns:
            DataFrame of CNV calls
        """
        sample_name = os.path.basename(summary_file).replace('.mosdepth.summary.txt', '')
        logger.info(f"Detecting CNVs in {sample_name} using classic algorithm")
        
        # Validate input files
        validate_file_exists(summary_file, f"Summary file for {sample_name}")
        validate_file_exists(regions_file, f"Regions file for {sample_name}")
        
        # Read sample data
        genome_stats = self._read_summary_stats(summary_file)
        regions_df = read_regions_with_line_tracking(regions_file, debug=logger.isEnabledFor(logging.DEBUG))
        
        # Detect sex if not provided
        if sample_sex is None:
            sample_sex = self._detect_sex(genome_stats)
            logger.info(f"Detected sex: {sample_sex}")
        
        # Calculate log2 ratios
        log2_ratios = self._calculate_log2_ratios(regions_df, genome_stats, gc_content, mappability)
        
        # Perform segmentation
        segments = self._segment_genome(log2_ratios)
        
        # Call CNVs with proper ploidy handling
        cnv_calls = self._call_cnvs(segments, sample_sex)
        
        return cnv_calls
    
    def _read_summary_stats(self, summary_file: str) -> Dict:
        """Read mosdepth summary statistics with validation."""
        validate_file_exists(summary_file, "Summary file")
        
        stats = {}
        try:
            with open(summary_file, 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        if parts[0] == 'total':
                            stats['mean_coverage'] = float(parts[3])
                        elif parts[0] in AUTOSOME_CHROMS + SEX_CHROMS + [MITO_CHROM]:
                            stats[f'{parts[0]}_mean'] = float(parts[3])
        except Exception as e:
            raise ValueError(f"Error reading summary file {summary_file}: {e}")
            
        if 'mean_coverage' not in stats:
            raise ValueError(f"Invalid summary file format: {summary_file}")
            
        return stats
    
    def _detect_sex(self, genome_stats: Dict) -> str:
        """Detect sample sex from coverage ratios."""
        if 'chrX_mean' not in genome_stats or 'chrY_mean' not in genome_stats:
            return 'unknown'
            
        # Calculate coverage ratios relative to autosomes
        autosome_means = [genome_stats.get(f'{chrom}_mean', 0) for chrom in AUTOSOME_CHROMS]
        autosome_means = [x for x in autosome_means if x > 0]
        
        if not autosome_means:
            return 'unknown'
            
        autosome_mean = np.mean(autosome_means)
        
        x_ratio = safe_divide(genome_stats['chrX_mean'], autosome_mean, 0.0)
        y_ratio = safe_divide(genome_stats['chrY_mean'], autosome_mean, 0.0)
        
        if x_ratio > 0.8 and y_ratio < 0.1:
            return 'XX'
        elif x_ratio < 0.6 and y_ratio > 0.3:
            return 'XY'
        else:
            return 'unknown'
    
    def _calculate_log2_ratios(self, regions_df: pd.DataFrame, genome_stats: Dict,
                              gc_content: Optional[Dict[str, float]] = None,
                              mappability: Optional[Dict[str, float]] = None) -> pd.DataFrame:
        """Calculate log2 ratios with proper normalization and validation."""
        regions_df = regions_df.copy()
        genome_mean = genome_stats.get('mean_coverage', 1.0)
        
        if genome_mean <= 0:
            raise ValueError(f"Invalid genome mean coverage: {genome_mean}")
        
        # Initialize log2 ratios
        regions_df['log2_ratio'] = np.nan
        regions_df['z_score'] = np.nan
        
        matched_regions = 0
        for idx, row in regions_df.iterrows():
            key = f"{row['chrom']}:{row['start']}-{row['end']}"
            
            if key in self.baseline_dict:
                baseline = self.baseline_dict[key]
                
                # Normalize coverage
                normalized_cov = row['coverage'] / genome_mean
                
                # Calculate log2 ratio
                if baseline['median_normalized'] > 0:
                    ratio = normalized_cov / baseline['median_normalized']
                    # Add small pseudocount to avoid log(0)
                    log2_ratio = np.log2(max(ratio, 0.001))
                    regions_df.at[idx, 'log2_ratio'] = log2_ratio
                    
                    # Calculate Z-score
                    if baseline['mad'] > 0:
                        z_score = (normalized_cov - baseline['median_normalized']) / baseline['mad']
                        regions_df.at[idx, 'z_score'] = z_score
                    
                    matched_regions += 1
        
        logger.info(f"Matched {matched_regions} regions with baseline out of {len(regions_df)}")
        
        if matched_regions == 0:
            raise ValueError("No regions matched with baseline")
        
        # Apply median centering for diploid regions
        self._apply_median_centering(regions_df)
        
        return regions_df
    
    def _apply_median_centering(self, regions_df: pd.DataFrame) -> None:
        """Center log2 ratios around diploid baseline."""
        # Use only autosomal regions for centering
        autosomal_mask = regions_df['chrom'].isin(AUTOSOME_CHROMS)
        autosomal_log2 = regions_df.loc[autosomal_mask & ~regions_df['log2_ratio'].isna(), 'log2_ratio']
        
        if len(autosomal_log2) > 0:
            # Use robust estimator (median) for centering
            median_log2 = np.median(autosomal_log2)
            regions_df['log2_ratio'] -= median_log2
            logger.debug(f"Applied median centering: {median_log2:.3f}")
    
    def _segment_genome(self, regions_df: pd.DataFrame) -> List[Dict]:
        """Perform genome segmentation using CBS-like algorithm."""
        segments = []
        
        # Process each chromosome separately
        for chrom in regions_df['chrom'].unique():
            chrom_data = regions_df[regions_df['chrom'] == chrom].copy()
            chrom_data = chrom_data.sort_values('start')
            
            # Skip if too few data points
            if len(chrom_data) < 5:
                continue
            
            # Simple segmentation based on changepoint detection
            chrom_segments = self._segment_chromosome(chrom_data)
            segments.extend(chrom_segments)
        
        logger.info(f"Created {len(segments)} segments across {len(regions_df['chrom'].unique())} chromosomes")
        return segments
    
    def _segment_chromosome(self, chrom_data: pd.DataFrame) -> List[Dict]:
        """Segment a single chromosome using changepoint detection."""
        segments = []
        log2_values = chrom_data['log2_ratio'].values
        positions = chrom_data['start'].values
        
        # Skip if too many missing values
        valid_mask = ~np.isnan(log2_values)
        if np.sum(valid_mask) < 5:
            return segments
        
        # Simple iterative segmentation
        changepoints = self._find_changepoints(log2_values[valid_mask])
        changepoints = [0] + changepoints + [len(log2_values)]
        
        for i in range(len(changepoints) - 1):
            start_idx = changepoints[i]
            end_idx = changepoints[i + 1]
            
            if end_idx - start_idx < 3:  # Skip very small segments
                continue
            
            segment_log2 = log2_values[start_idx:end_idx]
            valid_segment = segment_log2[~np.isnan(segment_log2)]
            
            if len(valid_segment) == 0:
                continue
            
            segment = {
                'chrom': chrom_data.iloc[0]['chrom'],
                'start': int(positions[start_idx]),
                'end': int(positions[end_idx - 1]) + DEFAULT_BIN_SIZE,
                'n_bins': end_idx - start_idx,
                'n_valid_bins': len(valid_segment),
                'mean_log2': np.mean(valid_segment),
                'median_log2': np.median(valid_segment),
                'mad_log2': np.median(np.abs(valid_segment - np.median(valid_segment)))
            }
            segments.append(segment)
        
        return segments
    
    def _find_changepoints(self, values: np.ndarray, min_segment_size: int = 3) -> List[int]:
        """Find changepoints using recursive binary segmentation."""
        changepoints = []
        n = len(values)
        
        if n < 2 * min_segment_size:
            return changepoints
        
        # Use recursive binary segmentation
        def binary_segment(start: int, end: int, depth: int = 0):
            if depth > 10 or end - start < 2 * min_segment_size:
                return
            
            max_stat = 0
            best_pos = -1
            
            # Test each possible split point
            for pos in range(start + min_segment_size, end - min_segment_size):
                left = values[start:pos]
                right = values[pos:end]
                
                # T-statistic for difference in means
                if len(left) > 1 and len(right) > 1:
                    try:
                        t_stat, p_value = stats.ttest_ind(left, right, equal_var=False)
                        if abs(t_stat) > max_stat:
                            max_stat = abs(t_stat)
                            best_pos = pos
                    except:
                        continue
            
            # If significant changepoint found, recursively segment
            if max_stat > 3.0 and best_pos > 0:  # Threshold for significance
                changepoints.append(best_pos)
                binary_segment(start, best_pos, depth + 1)
                binary_segment(best_pos, end, depth + 1)
        
        binary_segment(0, n)
        return sorted(changepoints)
    
    def _call_cnvs(self, segments: List[Dict], sample_sex: str) -> pd.DataFrame:
        """Call CNVs from segments with proper ploidy handling."""
        cnv_calls = []
        
        for segment in segments:
            # Determine expected ploidy
            expected_ploidy = self._get_expected_ploidy(segment['chrom'], sample_sex)
            
            # Calculate copy number from log2 ratio
            # log2(observed/expected) = log2(CN/ploidy)
            # CN = ploidy * 2^(log2_ratio)
            copy_number = expected_ploidy * (2 ** segment['mean_log2'])
            
            # Determine CNV type with ploidy-aware thresholds
            cnv_type = self._classify_cnv(copy_number, expected_ploidy, segment['mean_log2'])
            
            if cnv_type != 'neutral':
                cnv_calls.append({
                    'chrom': segment['chrom'],
                    'start': segment['start'],
                    'end': segment['end'],
                    'size': segment['end'] - segment['start'],
                    'type': cnv_type,
                    'copy_number': round(copy_number, 2),
                    'expected_copies': expected_ploidy,
                    'log2_ratio': round(segment['mean_log2'], 3),
                    'n_bins': segment['n_bins'],
                    'quality': self._calculate_segment_quality(segment)
                })
        
        # Merge adjacent CNVs of same type
        cnv_df = pd.DataFrame(cnv_calls)
        if not cnv_df.empty:
            cnv_df = self._merge_adjacent_cnvs(cnv_df)
        
        logger.info(f"Called {len(cnv_df)} CNVs from {len(segments)} segments")
        return cnv_df
    
    def _get_expected_ploidy(self, chrom: str, sample_sex: str) -> int:
        """Get expected copy number for chromosome given sex."""
        if chrom in AUTOSOME_CHROMS:
            return 2
        elif chrom == 'chrX':
            return 2 if sample_sex == 'XX' else 1
        elif chrom == 'chrY':
            return 0 if sample_sex == 'XX' else 1
        elif chrom == MITO_CHROM:
            return 100  # Approximate, highly variable
        else:
            return 2  # Default
    
    def _classify_cnv(self, copy_number: float, expected_ploidy: int, log2_ratio: float) -> str:
        """Classify CNV type based on copy number and expected ploidy."""
        # Use both absolute copy number and log2 ratio for classification
        
        # For log2 ratio thresholds
        del_threshold = -0.3  # log2(0.8) ≈ -0.32, allowing 20% deviation
        dup_threshold = 0.26  # log2(1.2) ≈ 0.26, allowing 20% deviation
        
        # For copy number thresholds (ploidy-aware)
        if expected_ploidy == 2:  # Diploid
            if copy_number < 1.5 and log2_ratio < del_threshold:
                return 'deletion'
            elif copy_number > 2.5 and log2_ratio > dup_threshold:
                return 'duplication'
        elif expected_ploidy == 1:  # Haploid (male X/Y)
            if copy_number < 0.5 and log2_ratio < del_threshold:
                return 'deletion'
            elif copy_number > 1.5 and log2_ratio > dup_threshold:
                return 'duplication'
        elif expected_ploidy == 0:  # Expected absence (female Y)
            if copy_number > 0.5:
                return 'duplication'
        
        return 'neutral'
    
    def _calculate_segment_quality(self, segment: Dict) -> float:
        """Calculate quality score for segment."""
        # Based on segment size and consistency
        size_score = min(1.0, segment['n_valid_bins'] / 10)
        
        # Based on deviation magnitude
        magnitude_score = min(1.0, abs(segment['mean_log2']) / 0.5)
        
        # Based on consistency (low MAD relative to mean)
        consistency_score = 1.0
        if segment['mad_log2'] > 0 and abs(segment['mean_log2']) > 0:
            consistency_score = max(0, 1 - segment['mad_log2'] / abs(segment['mean_log2']))
        
        return round(np.mean([size_score, magnitude_score, consistency_score]), 3)
    
    def _merge_adjacent_cnvs(self, cnv_df: pd.DataFrame, max_gap: int = 10000) -> pd.DataFrame:
        """Merge adjacent CNVs of the same type."""
        if cnv_df.empty:
            return cnv_df
        
        cnv_df = cnv_df.sort_values(['chrom', 'start'])
        merged = []
        
        for chrom in cnv_df['chrom'].unique():
            chrom_cnvs = cnv_df[cnv_df['chrom'] == chrom]
            
            current = None
            for _, cnv in chrom_cnvs.iterrows():
                if current is None:
                    current = cnv.to_dict()
                elif (cnv['type'] == current['type'] and 
                      cnv['start'] - current['end'] <= max_gap):
                    # Merge
                    current['end'] = cnv['end']
                    current['size'] = current['end'] - current['start']
                    current['n_bins'] += cnv['n_bins']
                    # Weighted average for numeric fields
                    w1 = current['n_bins'] - cnv['n_bins']
                    w2 = cnv['n_bins']
                    total_weight = current['n_bins']
                    current['copy_number'] = (current['copy_number'] * w1 + cnv['copy_number'] * w2) / total_weight
                    current['log2_ratio'] = (current['log2_ratio'] * w1 + cnv['log2_ratio'] * w2) / total_weight
                    current['quality'] = (current['quality'] * w1 + cnv['quality'] * w2) / total_weight
                else:
                    merged.append(current)
                    current = cnv.to_dict()
            
            if current is not None:
                merged.append(current)
        
        return pd.DataFrame(merged)

class OptimizedCNVDetector(FlexibleBaselineParser):
    """
    Optimized CNV detector using vectorized operations and efficient data structures.
    
    Provides fast and full modes for different performance/accuracy trade-offs.
    """
    
    def __init__(self, baseline_file: str, metadata_file: Optional[str] = None):
        """
        Initialize optimized CNV detector.
        
        Args:
            baseline_file: Path to baseline TSV file
            metadata_file: Optional metadata file path
        """
        logger.info("Loading baseline data...")
        self.baseline_df, self.format_info = self._read_baseline_flexible(baseline_file)
        
        # Create region key for efficient merging
        self.baseline_df['region_key'] = (self.baseline_df['chrom'] + ':' + 
                                          self.baseline_df['start'].astype(str) + '-' + 
                                          self.baseline_df['end'].astype(str))
        
        # Filter out excluded regions
        self.baseline_df = self.baseline_df[~self.baseline_df.get('excluded', False)]
        logger.info(f"Loaded {len(self.baseline_df)} baseline regions")
        
        # Load metadata if available
        self.metadata = {}
        if metadata_file and os.path.exists(metadata_file):
            with open(metadata_file, 'r') as f:
                self.metadata = json.load(f)
        elif baseline_file.endswith('.tsv'):
            default_metadata = baseline_file.replace('.tsv', '_metadata.json')
            if os.path.exists(default_metadata):
                with open(default_metadata, 'r') as f:
                    self.metadata = json.load(f)
    
    def detect_cnvs(self, summary_file: str, regions_file: str,
                   sample_sex: Optional[str] = None,
                   gc_content: Optional[Dict[str, float]] = None,
                   mappability: Optional[Dict[str, float]] = None,
                   fast_mode: bool = True) -> pd.DataFrame:
        """
        Detect CNVs using optimized algorithm with choice of speed vs accuracy.
        
        Args:
            summary_file: Path to sample summary file
            regions_file: Path to sample regions file
            sample_sex: Sample sex ('XX', 'XY', or None for auto-detection)
            gc_content: Optional GC content annotation
            mappability: Optional mappability annotation
            fast_mode: If True, use fast sliding window approach; if False, use full segmentation
            
        Returns:
            DataFrame of CNV calls
        """
        sample_name = os.path.basename(summary_file).replace('.mosdepth.summary.txt', '')
        logger.info(f"Detecting CNVs in {sample_name} using optimized algorithm")
        
        # Validate input files
        validate_file_exists(summary_file, f"Summary file for {sample_name}")
        validate_file_exists(regions_file, f"Regions file for {sample_name}")
        
        # Read sample data
        genome_stats = self._read_summary_stats(summary_file)
        
        # Read regions with progress tracking
        logger.info("Reading sample regions...")
        regions_df = self._read_regions_optimized(regions_file)
        logger.info(f"Loaded {len(regions_df)} sample regions")
        
        # Detect sex if not provided
        if sample_sex is None:
            sample_sex = self._detect_sex(genome_stats)
            logger.info(f"Detected sex: {sample_sex}")
        
        # Calculate log2 ratios using vectorized operations
        logger.info("Calculating log2 ratios...")
        log2_ratios_df = self._calculate_log2_ratios_vectorized(regions_df, genome_stats, gc_content, mappability)
        
        if fast_mode:
            # Use simplified bin-based CNV calling for speed
            logger.info("Using fast mode CNV detection...")
            cnv_calls = self._fast_cnv_detection(log2_ratios_df, sample_sex)
        else:
            # Perform full segmentation
            logger.info("Performing segmentation...")
            segments = self._segment_genome_optimized(log2_ratios_df)
            
            # Call CNVs with proper ploidy handling
            logger.info("Calling CNVs...")
            cnv_calls = self._call_cnvs(segments, sample_sex)
        
        return cnv_calls
    
    def _read_summary_stats(self, summary_file: str) -> Dict:
        """Read mosdepth summary statistics with validation."""
        validate_file_exists(summary_file, "Summary file")
        
        stats = {}
        try:
            with open(summary_file, 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        if parts[0] == 'total':
                            stats['mean_coverage'] = float(parts[3])
                        elif parts[0] in AUTOSOME_CHROMS + SEX_CHROMS + [MITO_CHROM]:
                            stats[f'{parts[0]}_mean'] = float(parts[3])
        except Exception as e:
            raise ValueError(f"Error reading summary file {summary_file}: {e}")
            
        if 'mean_coverage' not in stats:
            raise ValueError(f"Invalid summary file format: {summary_file}")
            
        return stats
    
    def _read_regions_optimized(self, regions_file: str) -> pd.DataFrame:
        """Read mosdepth regions file with optimized memory usage."""
        validate_file_exists(regions_file, "Regions file")
        
        # Check file size to decide strategy
        file_size = os.path.getsize(regions_file)
        
        # If file is small enough (< 100MB compressed), use simple reading
        if file_size < 100 * 1024 * 1024:
            regions_df = read_regions_with_line_tracking(regions_file, debug=logger.isEnabledFor(logging.DEBUG))
        else:
            # For large files, read in chunks
            chunks = []
            chunk_data = []
            line_num = 0
            
            try:
                with gzip.open(regions_file, 'rt') as f:
                    for i, line in enumerate(f):
                        line_num = i + 1
                        parts = line.strip().split('\t')
                        
                        if len(parts) < 4:
                            logger.debug(f"Line {line_num}: Skipping line with {len(parts)} columns")
                            continue
                        
                        try:
                            chunk_data.append({
                                'chrom': parts[0],
                                'start': int(parts[1]),
                                'end': int(parts[2]),
                                'coverage': float(parts[3])
                            })
                        except ValueError as e:
                            logger.warning(f"Line {line_num} in {regions_file}: {e}")
                            continue
                        
                        # Process in chunks to manage memory
                        if (i + 1) % CHUNK_SIZE == 0:
                            chunks.append(pd.DataFrame(chunk_data))
                            chunk_data = []
                    
                    # Add remaining data
                    if chunk_data:
                        chunks.append(pd.DataFrame(chunk_data))
                        
            except Exception as e:
                logger.error(f"Error reading {regions_file} at line {line_num}: {e}")
                raise
            
            # Concatenate all chunks
            regions_df = pd.concat(chunks, ignore_index=True)
            logger.debug(f"Successfully read {len(regions_df)} regions from {regions_file}")
        
        # Create region key for merging
        regions_df['region_key'] = (regions_df['chrom'] + ':' + 
                                regions_df['start'].astype(str) + '-' + 
                                regions_df['end'].astype(str))
        
        return regions_df
        
    def _detect_sex(self, genome_stats: Dict) -> str:
        """Detect sample sex from coverage ratios."""
        if 'chrX_mean' not in genome_stats or 'chrY_mean' not in genome_stats:
            return 'unknown'
            
        # Calculate coverage ratios relative to autosomes
        autosome_means = [genome_stats.get(f'{chrom}_mean', 0) for chrom in AUTOSOME_CHROMS]
        autosome_means = [x for x in autosome_means if x > 0]
        
        if not autosome_means:
            return 'unknown'
            
        autosome_mean = np.mean(autosome_means)
        
        x_ratio = safe_divide(genome_stats['chrX_mean'], autosome_mean, 0.0)
        y_ratio = safe_divide(genome_stats['chrY_mean'], autosome_mean, 0.0)
        
        if x_ratio > 0.8 and y_ratio < 0.1:
            return 'XX'
        elif x_ratio < 0.6 and y_ratio > 0.3:
            return 'XY'
        else:
            return 'unknown'
    
    def _calculate_log2_ratios_vectorized(self, regions_df: pd.DataFrame, genome_stats: Dict,
                                         gc_content: Optional[Dict[str, float]] = None,
                                         mappability: Optional[Dict[str, float]] = None) -> pd.DataFrame:
        """Calculate log2 ratios using vectorized operations."""
        
        # Normalize coverage
        genome_mean = genome_stats.get('mean_coverage', 1.0)
        if genome_mean <= 0:
            raise ValueError(f"Invalid genome mean coverage: {genome_mean}")
            
        regions_df['normalized_coverage'] = regions_df['coverage'] / genome_mean
        
        # Merge with baseline using vectorized operation
        logger.info("Merging with baseline...")
        merged_df = pd.merge(
            regions_df[['chrom', 'start', 'end', 'region_key', 'normalized_coverage']],
            self.baseline_df[['region_key', 'median_normalized', 'mad']],
            on='region_key',
            how='inner'
        )
        
        logger.info(f"Merged {len(merged_df)} regions with baseline")
        
        if merged_df.empty:
            raise ValueError("No regions matched with baseline")
        
        # Vectorized log2 ratio calculation
        merged_df['ratio'] = merged_df['normalized_coverage'] / merged_df['median_normalized']
        merged_df['log2_ratio'] = np.log2(np.maximum(merged_df['ratio'], 0.001))
        
        # Vectorized Z-score calculation
        merged_df['z_score'] = ((merged_df['normalized_coverage'] - merged_df['median_normalized']) / 
                                np.maximum(merged_df['mad'], 0.001))
        
        # Apply median centering for diploid regions
        self._apply_median_centering_vectorized(merged_df)
        
        return merged_df
    
    def _apply_median_centering_vectorized(self, merged_df: pd.DataFrame) -> None:
        """Center log2 ratios around diploid baseline using vectorized operations."""
        # Use only autosomal regions for centering
        autosomal_mask = merged_df['chrom'].isin(AUTOSOME_CHROMS)
        autosomal_log2 = merged_df.loc[autosomal_mask, 'log2_ratio']
        
        if len(autosomal_log2) > 0:
            # Use robust estimator (median) for centering
            median_log2 = np.median(autosomal_log2)
            merged_df['log2_ratio'] -= median_log2
            logger.debug(f"Applied median centering: {median_log2:.3f}")
    
    def _segment_genome_optimized(self, log2_ratios_df: pd.DataFrame) -> List[Dict]:
        """Perform genome segmentation with progress tracking."""
        segments = []
        
        # Group by chromosome for efficient processing
        grouped = log2_ratios_df.groupby('chrom')
        
        # Process each chromosome with progress bar
        for chrom, chrom_data in tqdm(grouped, desc="Segmenting chromosomes"):
            chrom_data = chrom_data.sort_values('start').reset_index(drop=True)
            
            # Skip if too few data points
            if len(chrom_data) < 5:
                continue
            
            # Vectorized segmentation
            chrom_segments = self._segment_chromosome_vectorized(chrom_data)
            segments.extend(chrom_segments)
        
        return segments
    
    def _segment_chromosome_vectorized(self, chrom_data: pd.DataFrame) -> List[Dict]:
        """Segment a single chromosome using vectorized changepoint detection."""
        segments = []
        log2_values = chrom_data['log2_ratio'].values
        
        # Skip if too many missing values
        valid_mask = ~np.isnan(log2_values)
        if np.sum(valid_mask) < 5:
            return segments
        
        # Use vectorized changepoint detection
        changepoints = self._find_changepoints_vectorized(log2_values[valid_mask])
        changepoints = [0] + changepoints + [len(log2_values)]
        
        # Create segments
        for i in range(len(changepoints) - 1):
            start_idx = changepoints[i]
            end_idx = changepoints[i + 1]
            
            if end_idx - start_idx < 3:  # Skip very small segments
                continue
            
            segment_data = chrom_data.iloc[start_idx:end_idx]
            segment_log2 = segment_data['log2_ratio'].values
            valid_segment = segment_log2[~np.isnan(segment_log2)]
            
            if len(valid_segment) == 0:
                continue
            
            segment = {
                'chrom': chrom_data.iloc[0]['chrom'],
                'start': int(segment_data.iloc[0]['start']),
                'end': int(segment_data.iloc[-1]['end']),
                'n_bins': end_idx - start_idx,
                'n_valid_bins': len(valid_segment),
                'mean_log2': np.mean(valid_segment),
                'median_log2': np.median(valid_segment),
                'mad_log2': np.median(np.abs(valid_segment - np.median(valid_segment)))
            }
            segments.append(segment)
        
        return segments
    
    def _find_changepoints_vectorized(self, values: np.ndarray, min_segment_size: int = 10) -> List[int]:
        """Find changepoints using a fast sliding window approach."""
        changepoints = []
        n = len(values)
        
        if n < 2 * min_segment_size:
            return changepoints
        
        # Use a sliding window approach for speed
        window_size = max(min_segment_size * 2, 50)
        step_size = max(min_segment_size, 25)
        
        # Calculate rolling statistics
        for i in range(window_size, n - window_size, step_size):
            left_window = values[max(0, i - window_size):i]
            right_window = values[i:min(n, i + window_size)]
            
            # Quick statistical test
            if len(left_window) > 5 and len(right_window) > 5:
                left_mean = np.mean(left_window)
                right_mean = np.mean(right_window)
                pooled_var = (np.var(left_window) + np.var(right_window)) / 2
                
                # T-statistic approximation
                if pooled_var > 0:
                    pooled_std = np.sqrt(pooled_var)
                    t_stat = abs(left_mean - right_mean) / (pooled_std * np.sqrt(2 / window_size))
                    if t_stat > 3.0:  # Significant difference
                        changepoints.append(i)
        
        # Merge nearby changepoints
        if changepoints:
            merged = [changepoints[0]]
            for cp in changepoints[1:]:
                if cp - merged[-1] > min_segment_size * 2:
                    merged.append(cp)
            changepoints = merged
        
        return changepoints
    
    def _fast_cnv_detection(self, log2_ratios_df: pd.DataFrame, sample_sex: str) -> pd.DataFrame:
        """Fast CNV detection using sliding windows without full segmentation."""
        cnv_calls = []
        
        # Parameters for sliding window
        window_size = 10  # 10 bins = 5kb
        threshold_del = -0.3
        threshold_dup = 0.26
        min_cnv_size = 5000  # 5kb minimum
        
        # Group by chromosome
        for chrom, chrom_data in log2_ratios_df.groupby('chrom'):
            chrom_data = chrom_data.sort_values('start').reset_index(drop=True)
            
            if len(chrom_data) < window_size:
                continue
                
            # Calculate rolling statistics
            chrom_data['mean_log2_roll'] = chrom_data['log2_ratio'].rolling(
                window=window_size, center=True, min_periods=window_size//2
            ).mean()
            
            # Identify CNV regions
            chrom_data['cnv_type'] = 'neutral'
            chrom_data.loc[chrom_data['mean_log2_roll'] <= threshold_del, 'cnv_type'] = 'deletion'
            chrom_data.loc[chrom_data['mean_log2_roll'] >= threshold_dup, 'cnv_type'] = 'duplication'
            
            # Group consecutive CNV bins
            chrom_data['cnv_group'] = (chrom_data['cnv_type'] != chrom_data['cnv_type'].shift()).cumsum()
            
            # Create CNV calls from groups
            for (cnv_type, cnv_group), group_data in chrom_data.groupby(['cnv_type', 'cnv_group']):
                if cnv_type == 'neutral':
                    continue
                    
                start = group_data.iloc[0]['start']
                end = group_data.iloc[-1]['end']
                size = end - start
                
                if size >= min_cnv_size:
                    expected_ploidy = self._get_expected_ploidy(chrom, sample_sex)
                    mean_log2 = group_data['log2_ratio'].mean()
                    copy_number = expected_ploidy * (2 ** mean_log2)
                    
                    cnv_calls.append({
                        'chrom': chrom,
                        'start': int(start),
                        'end': int(end),
                        'size': int(size),
                        'type': cnv_type,
                        'copy_number': round(copy_number, 2),
                        'expected_copies': expected_ploidy,
                        'log2_ratio': round(mean_log2, 3),
                        'n_bins': len(group_data),
                        'quality': round(min(1.0, len(group_data) / 10), 3)
                    })
        
        cnv_df = pd.DataFrame(cnv_calls)
        if not cnv_df.empty:
            cnv_df = self._merge_adjacent_cnvs_vectorized(cnv_df)
        
        return cnv_df
    
    def _call_cnvs(self, segments: List[Dict], sample_sex: str) -> pd.DataFrame:
        """Call CNVs from segments with proper ploidy handling."""
        cnv_calls = []
        
        for segment in segments:
            # Determine expected ploidy
            expected_ploidy = self._get_expected_ploidy(segment['chrom'], sample_sex)
            
            # Calculate copy number from log2 ratio
            copy_number = expected_ploidy * (2 ** segment['mean_log2'])
            
            # Determine CNV type with ploidy-aware thresholds
            cnv_type = self._classify_cnv(copy_number, expected_ploidy, segment['mean_log2'])
            
            if cnv_type != 'neutral':
                cnv_calls.append({
                    'chrom': segment['chrom'],
                    'start': segment['start'],
                    'end': segment['end'],
                    'size': segment['end'] - segment['start'],
                    'type': cnv_type,
                    'copy_number': round(copy_number, 2),
                    'expected_copies': expected_ploidy,
                    'log2_ratio': round(segment['mean_log2'], 3),
                    'n_bins': segment['n_bins'],
                    'quality': self._calculate_segment_quality(segment)
                })
        
        # Merge adjacent CNVs of same type
        cnv_df = pd.DataFrame(cnv_calls)
        if not cnv_df.empty:
            cnv_df = self._merge_adjacent_cnvs_vectorized(cnv_df)
        
        return cnv_df
    
    def _get_expected_ploidy(self, chrom: str, sample_sex: str) -> int:
        """Get expected copy number for chromosome given sex."""
        if chrom in AUTOSOME_CHROMS:
            return 2
        elif chrom == 'chrX':
            return 2 if sample_sex == 'XX' else 1
        elif chrom == 'chrY':
            return 0 if sample_sex == 'XX' else 1
        elif chrom == MITO_CHROM:
            return 100  # Approximate, highly variable
        else:
            return 2  # Default
    
    def _classify_cnv(self, copy_number: float, expected_ploidy: int, log2_ratio: float) -> str:
        """Classify CNV type based on copy number and expected ploidy."""
        # Use both absolute copy number and log2 ratio for classification
        
        # For log2 ratio thresholds
        del_threshold = -0.3  # log2(0.8) ≈ -0.32, allowing 20% deviation
        dup_threshold = 0.26  # log2(1.2) ≈ 0.26, allowing 20% deviation
        
        # For copy number thresholds (ploidy-aware)
        if expected_ploidy == 2:  # Diploid
            if copy_number < 1.5 and log2_ratio < del_threshold:
                return 'deletion'
            elif copy_number > 2.5 and log2_ratio > dup_threshold:
                return 'duplication'
        elif expected_ploidy == 1:  # Haploid (male X/Y)
            if copy_number < 0.5 and log2_ratio < del_threshold:
                return 'deletion'
            elif copy_number > 1.5 and log2_ratio > dup_threshold:
                return 'duplication'
        elif expected_ploidy == 0:  # Expected absence (female Y)
            if copy_number > 0.5:
                return 'duplication'
        
        return 'neutral'
    
    def _calculate_segment_quality(self, segment: Dict) -> float:
        """Calculate quality score for segment."""
        # Based on segment size and consistency
        size_score = min(1.0, segment['n_valid_bins'] / 10)
        
        # Based on deviation magnitude
        magnitude_score = min(1.0, abs(segment['mean_log2']) / 0.5)
        
        # Based on consistency (low MAD relative to mean)
        consistency_score = 1.0
        if segment['mad_log2'] > 0 and abs(segment['mean_log2']) > 0:
            consistency_score = max(0, 1 - segment['mad_log2'] / abs(segment['mean_log2']))
        
        return round(np.mean([size_score, magnitude_score, consistency_score]), 3)
    
    def _merge_adjacent_cnvs_vectorized(self, cnv_df: pd.DataFrame, max_gap: int = 10000) -> pd.DataFrame:
        """Merge adjacent CNVs using vectorized operations."""
        if cnv_df.empty:
            return cnv_df
        
        # Sort by chromosome and start position
        cnv_df = cnv_df.sort_values(['chrom', 'start']).reset_index(drop=True)
        
        # Mark groups of adjacent CNVs
        cnv_df['gap_to_prev'] = cnv_df['start'] - cnv_df['end'].shift(1)
        cnv_df['same_type'] = cnv_df['type'] == cnv_df['type'].shift(1)
        cnv_df['same_chrom'] = cnv_df['chrom'] == cnv_df['chrom'].shift(1)
        cnv_df['should_merge'] = (cnv_df['gap_to_prev'] <= max_gap) & cnv_df['same_type'] & cnv_df['same_chrom']
        
        # Create merge groups
        cnv_df['merge_group'] = (~cnv_df['should_merge']).cumsum()
        
        # Aggregate by merge group
        merged = cnv_df.groupby('merge_group').agg({
            'chrom': 'first',
            'start': 'min',
            'end': 'max',
            'type': 'first',
            'expected_copies': 'first',
            'copy_number': lambda x: np.average(x, weights=cnv_df.loc[x.index, 'n_bins']),
            'log2_ratio': lambda x: np.average(x, weights=cnv_df.loc[x.index, 'n_bins']),
            'n_bins': 'sum',
            'quality': lambda x: np.average(x, weights=cnv_df.loc[x.index, 'n_bins'])
        }).reset_index(drop=True)
        
        # Recalculate size
        merged['size'] = merged['end'] - merged['start']
        
        # Round numeric values
        merged['copy_number'] = merged['copy_number'].round(2)
        merged['log2_ratio'] = merged['log2_ratio'].round(3)
        merged['quality'] = merged['quality'].round(3)
        
        # Reorder columns
        return merged[['chrom', 'start', 'end', 'size', 'type', 'copy_number', 
                      'expected_copies', 'log2_ratio', 'n_bins', 'quality']]

def print_enhanced_cnv_summary(cnv_df: pd.DataFrame, sample_name: str) -> None:
    """Print detailed CNV summary statistics with enhanced formatting."""
    print(f"\n{'='*60}")
    print(f"CNV Detection Summary for {sample_name}")
    print(f"{'='*60}")
    
    if cnv_df.empty:
        print("\nNo CNVs detected")
        return
    
    print(f"\nTotal CNVs detected: {len(cnv_df)}")
    print(f"  Deletions: {len(cnv_df[cnv_df['type'] == 'deletion'])}")
    print(f"  Duplications: {len(cnv_df[cnv_df['type'] == 'duplication'])}")
    
    # Per-chromosome breakdown
    print("\nPer-chromosome breakdown:")
    chrom_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    for chrom in chrom_order:
        if chrom in cnv_df['chrom'].values:
            chrom_cnvs = cnv_df[cnv_df['chrom'] == chrom]
            dels = len(chrom_cnvs[chrom_cnvs['type'] == 'deletion'])
            dups = len(chrom_cnvs[chrom_cnvs['type'] == 'duplication'])
            total_size = chrom_cnvs['size'].sum()
            print(f"  {chrom}: {len(chrom_cnvs)} CNVs ({dels} del, {dups} dup) - {total_size:,} bp total")
    
    # Size distribution
    print("\nSize distribution:")
    size_ranges = [(0, 1000), (1000, 10000), (10000, 100000), 
                   (100000, 1000000), (1000000, float('inf'))]
    labels = ['<1kb', '1-10kb', '10-100kb', '100kb-1Mb', '>1Mb']
    
    for (min_size, max_size), label in zip(size_ranges, labels):
        count = len(cnv_df[(cnv_df['size'] >= min_size) & (cnv_df['size'] < max_size)])
        if count > 0:
            print(f"  {label}: {count} CNVs")
    
    # Quality distribution
    print("\nQuality distribution:")
    if 'quality' in cnv_df.columns:
        q_bins = [0, 0.5, 0.7, 0.9, 1.0]
        q_labels = ['Low (0-0.5)', 'Medium (0.5-0.7)', 'High (0.7-0.9)', 'Very High (0.9-1.0)']
        for i in range(len(q_bins)-1):
            count = len(cnv_df[(cnv_df['quality'] >= q_bins[i]) & 
                              (cnv_df['quality'] < q_bins[i+1])])
            if count > 0:
                print(f"  {q_labels[i]}: {count} CNVs")
    
    # Top CNVs
    print(f"\nLargest CNVs:")
    for _, cnv in cnv_df.nlargest(5, 'size').iterrows():
        print(f"  {cnv['chrom']}:{cnv['start']:,}-{cnv['end']:,} "
              f"({cnv['size']:,} bp) {cnv['type']} "
              f"CN={cnv.get('copy_number', 'N/A'):.2f} "
              f"Quality={cnv.get('quality', 'N/A')}")

def main():
    """Main entry point for the integrated CNV pipeline with comprehensive error handling."""
    parser = argparse.ArgumentParser(
        description='Integrated CNV Detection Pipeline: BAM/CRAM → mosdepth → baseline → CNV detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Complete Workflow Examples:

# 1. Run mosdepth on BAM/CRAM files
python cnv_pipeline_integrated.py mosdepth sample.bam ./mosdepth_output/
python cnv_pipeline_integrated.py mosdepth bam_files.txt ./mosdepth_output/ --parallel 4

# 2. Create baseline from mosdepth output
python cnv_pipeline_integrated.py baseline baseline.tsv ./mosdepth_output/*.summary.txt

# 3. Detect CNVs using different modes
python cnv_pipeline_integrated.py detect sample.summary.txt baseline.tsv --mode fast
python cnv_pipeline_integrated.py detect sample.summary.txt baseline.tsv --mode classic --show-summary

# 4. Complete workflow
python cnv_pipeline_integrated.py mosdepth samples.txt ./mosdepth_out/ --parallel 2
python cnv_pipeline_integrated.py baseline baseline.tsv ./mosdepth_out/*.summary.txt --min-samples 20
python cnv_pipeline_integrated.py detect sample.summary.txt baseline.tsv --mode fast --show-summary

Key Features:
- Integrated mosdepth processing with parallel support
- Robust baseline generation with checkpoint/resume
- Classic and optimized CNV detection algorithms  
- Comprehensive QC and error handling
- Memory-efficient processing for large datasets
        """
    )
    
    # Add global debug flag
    parser.add_argument('--debug', action='store_true', 
                       help='Enable debug mode with verbose output')
    
    subparsers = parser.add_subparsers(dest='command', help='Pipeline commands')
    
    # Mosdepth processing command
    mosdepth_parser = subparsers.add_parser(
        'mosdepth', 
        help='Run mosdepth on BAM/CRAM files with specified bin size',
        description='Run mosdepth on BAM/CRAM files with specified bin size',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    # Process a single BAM file
    python cnv_pipeline_integrated.py mosdepth sample.bam output_dir/
    
    # Process a single CRAM file (requires reference)
    python cnv_pipeline_integrated.py mosdepth sample.cram output_dir/ --reference genome.fa
    
    # Process multiple files from a list
    python cnv_pipeline_integrated.py mosdepth alignment_files.txt output_dir/
    
    # With custom parameters
    python cnv_pipeline_integrated.py mosdepth --bin-size 1000 --threads 8 alignment_files.txt output_dir/
    python cnv_pipeline_integrated.py mosdepth --parallel 4 sample.bam output_dir/
    python cnv_pipeline_integrated.py mosdepth --no-fast-mode sample.cram output_dir/ --reference genome.fa

Note: CRAM files require a reference genome. Either use --reference or set the REF_PATH 
environment variable to the reference genome path.
        """
    )
    mosdepth_parser.add_argument('alignment_input', 
                                help='Single BAM/CRAM file or text file containing alignment file paths (one per line)')
    mosdepth_parser.add_argument('output_dir', 
                                help='Directory where mosdepth output files will be saved')
    mosdepth_parser.add_argument('--reference', '--fasta', '-f',
                                help='Reference genome FASTA file (required for CRAM files)')
    mosdepth_parser.add_argument('--bin-size', type=int, default=500,
                                help='Bin size for coverage calculation (default: 500)')
    mosdepth_parser.add_argument('--threads', type=int, default=4,
                                help='Number of threads per mosdepth job (default: 4)')
    mosdepth_parser.add_argument('--parallel', type=int, default=1,
                                help='Number of parallel mosdepth jobs (default: 1, sequential)')
    mosdepth_parser.add_argument('--no-fast-mode', action='store_true',
                                help='Disable fast mode for more accurate coverage calculation (slower)')
    
    # Baseline creation command
    baseline_parser = subparsers.add_parser('baseline', help='Create baseline from multiple samples')
    baseline_parser.add_argument('output', help='Output baseline file (TSV format)')
    baseline_parser.add_argument('samples', nargs='+', help='Sample summary files (.mosdepth.summary.txt)')
    baseline_parser.add_argument('--min-samples', type=int, default=20,
                                help='Minimum samples for baseline (default: 20)')
    baseline_parser.add_argument('--gc-content', help='GC content file (optional)')
    baseline_parser.add_argument('--mappability', help='Mappability file (optional)')
    baseline_parser.add_argument('--checkpoint-dir', help='Directory for saving checkpoints (enables resume capability)')
    baseline_parser.add_argument('--resume', action='store_true', help='Resume from checkpoint (requires --checkpoint-dir)')
    baseline_parser.add_argument('--chunk-size', type=int, default=50, 
                                help='Number of samples to process per chunk (default: 50)')
    
    # CNV detection command
    detect_parser = subparsers.add_parser('detect', help='Detect CNVs in a sample')
    detect_parser.add_argument('summary', help='Sample .mosdepth.summary.txt file')
    detect_parser.add_argument('baseline', help='Baseline TSV file')
    detect_parser.add_argument('--output', help='Output CNV file (default: sample.cnvs.tsv)')
    detect_parser.add_argument('--sex', choices=['XX', 'XY'], help='Sample sex (auto-detected if not provided)')
    detect_parser.add_argument('--gc-content', help='GC content file (optional)')
    detect_parser.add_argument('--mappability', help='Mappability file (optional)')
    detect_parser.add_argument('--mode', choices=['classic', 'fast', 'full'], default='fast',
                              help='Detection mode: classic (v2/v4 algorithm), fast (optimized sliding window), full (optimized with segmentation)')
    detect_parser.add_argument('--show-summary', action='store_true',
                              help='Show detailed CNV summary statistics')
    
    # QC command
    qc_parser = subparsers.add_parser('qc', help='Run QC on samples')
    qc_parser.add_argument('samples', nargs='+', help='Sample summary files')
    qc_parser.add_argument('--output', help='QC report output file')
    
    args = parser.parse_args()
    
    # Set debug level if requested
    if args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug mode enabled")
    
    try:
        if args.command == 'mosdepth':
            # Run mosdepth processing
            
            # Validate reference file if provided
            if args.reference:
                ref_path = Path(args.reference)
                if not ref_path.exists():
                    logger.error(f"Reference file does not exist: {args.reference}")
                    sys.exit(1)
                if not ref_path.is_file():
                    logger.error(f"Reference path is not a file: {args.reference}")
                    sys.exit(1)
            
            # Create output directory
            output_dir = Path(args.output_dir)
            try:
                output_dir.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                logger.error(f"Failed to create output directory: {str(e)}")
                sys.exit(1)
            
            # Set up logging with file output
            setup_logging_with_file(args.output_dir, "mosdepth_batch")
            
            # Initialize processor
            processor = MosdepthProcessor(args.output_dir)
            
            # Check mosdepth installation
            is_installed, version_or_error = processor.check_mosdepth_installation()
            if not is_installed:
                logger.error(f"Mosdepth installation check failed: {version_or_error}")
                logger.error("Please install mosdepth before running this script")
                sys.exit(1)
            
            logger.info(f"Found mosdepth: {version_or_error}")
            
            # Run processing
            exit_code = processor.process_single_or_batch(
                args.alignment_input,
                bin_size=args.bin_size,
                threads=args.threads,
                max_workers=args.parallel,
                fast_mode=not args.no_fast_mode,
                reference=args.reference
            )
            
            sys.exit(exit_code)
            
        elif args.command == 'baseline':
            # Create baseline
            checkpoint_dir = getattr(args, 'checkpoint_dir', None)
            generator = BaselineGenerator(min_samples=args.min_samples, checkpoint_dir=checkpoint_dir)
            
            # Load optional annotation files
            gc_content = None
            mappability = None
            if args.gc_content:
                # TODO: Implement GC content loading based on file format
                logger.info(f"GC content file specified but not yet implemented: {args.gc_content}")
            if args.mappability:
                # TODO: Implement mappability loading based on file format
                logger.info(f"Mappability file specified but not yet implemented: {args.mappability}")
            
            # Use chunked processing
            chunk_size = getattr(args, 'chunk_size', 50)
            logger.info(f"Processing samples in chunks of {chunk_size}")
            generator.add_samples_chunked(args.samples, chunk_size=chunk_size, gc_content=gc_content, mappability=mappability)
            
            # Finalize baseline
            generator.finalize_baseline(args.output)
            logger.info("Baseline creation completed successfully")
                
        elif args.command == 'detect':
            # Detect CNVs
            validate_file_exists(args.baseline, "Baseline file")
                
            regions_file = args.summary.replace('.mosdepth.summary.txt', '.regions.bed.gz')
            validate_file_exists(regions_file, "Regions file")
            
            # Load optional annotation files
            gc_content = None
            mappability = None
            if args.gc_content:
                # TODO: Implement GC content loading
                logger.info(f"GC content file specified but not yet implemented: {args.gc_content}")
            if args.mappability:
                # TODO: Implement mappability loading
                logger.info(f"Mappability file specified but not yet implemented: {args.mappability}")
            
            # Choose detector based on mode
            start_time = time.time()
            mode = getattr(args, 'mode', 'fast')
            
            if mode == 'classic':
                # Use classic detector with recursive binary segmentation
                detector = CNVDetector(args.baseline)
                cnvs = detector.detect_cnvs(
                    args.summary,
                    regions_file,
                    sample_sex=args.sex,
                    gc_content=gc_content,
                    mappability=mappability
                )
            else:
                # Use optimized detector
                detector = OptimizedCNVDetector(args.baseline)
                fast_mode = (mode == 'fast')
                cnvs = detector.detect_cnvs(
                    args.summary,
                    regions_file,
                    sample_sex=args.sex,
                    gc_content=gc_content,
                    mappability=mappability,
                    fast_mode=fast_mode
                )
            
            # Save results
            output_file = args.output or args.summary.replace('.mosdepth.summary.txt', '.cnvs.tsv')
            cnvs.to_csv(output_file, sep='\t', index=False)
            
            elapsed_time = time.time() - start_time
            logger.info(f"Detected {len(cnvs)} CNVs in {elapsed_time:.1f} seconds using {mode} mode")
            logger.info(f"Results saved to {output_file}")
            
            # Get sample name for summary
            sample_name = os.path.basename(args.summary).replace('.mosdepth.summary.txt', '')
            
            # Print enhanced summary if requested or if CNVs found
            if getattr(args, 'show_summary', False) or not cnvs.empty:
                print_enhanced_cnv_summary(cnvs, sample_name)
                
        elif args.command == 'qc':
            # Run QC analysis
            logger.info("QC analysis functionality not yet implemented")
            # TODO: Implement comprehensive QC analysis
                
        else:
            parser.print_help()
            
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    # Handle Ctrl+C gracefully
    signal.signal(signal.SIGINT, lambda x, y: sys.exit(130))
    main()
