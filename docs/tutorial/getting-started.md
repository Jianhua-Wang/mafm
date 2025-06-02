# Getting Started

## What is MAFM?

MAFM (Multi-Ancestry Fine-Mapping) is a comprehensive pipeline for performing statistical fine-mapping analysis across multiple ancestries and cohorts. It provides a unified framework for:

- **Quality Control**: Assess the reliability of summary statistics and LD matrices
- **Meta-Analysis**: Combine data across populations and cohorts  
- **Fine-Mapping**: Identify causal variants using multiple statistical methods
- **Post-Processing**: Combine and interpret results across studies

## The MAFM Framework

MAFM's workflow can be visualized as follows:

```mermaid
graph TD
    A[Multiple GWAS<br/>Summary Statistics] --> B[Quality Control]
    C[Multiple LD<br/>Matrices] --> B
    B --> D{Meta-Analysis Strategy}
    D -->|Cross-Ancestry| E[Cross-Ancestry Meta]
    D -->|Within-Ancestry| F[Within-Ancestry Meta]
    D -->|No Meta| G[Keep Separate]
    E --> H[Fine-Mapping]
    F --> H
    G --> H
    H --> I{Fine-Mapping Strategy}
    I -->|Single Input| J[Run on Each<br/>Study Separately]
    I -->|Multi Input| K[Run on Combined<br/>Studies Together]
    I -->|Post-hoc Combine| L[Run Separately<br/>Then Combine Results]
    J --> M[Results Integration]
    K --> M
    L --> M
    M --> N[Credible Sets &<br/>Posterior Probabilities]
```

## Core Concepts

### 1. Input Data Types

MAFM works with two main types of input data:

- **Summary Statistics**: GWAS results containing effect sizes, standard errors, and p-values
- **LD Matrices**: Linkage disequilibrium correlation matrices from reference panels

### 2. Meta-Analysis Strategies

!!! info "Meta-Analysis Methods"
    MAFM supports three meta-analysis approaches:
    
    - **`meta_all`**: Combine all studies regardless of ancestry
    - **`meta_by_population`**: Combine studies within each ancestry separately  
    - **`no_meta`**: Keep all studies separate

### 3. Fine-Mapping Strategies

MAFM offers three complementary fine-mapping strategies:

!!! note "Strategy Comparison"
    
    **Single Input** (`single_input`)
    : Use traditional fine-mapping tools that analyze one study at a time
    : Best for: Well-powered single studies, ancestry-specific analysis
    
    **Multi Input** (`multi_input`)  
    : Use tools designed to analyze multiple studies simultaneously
    : Best for: Leveraging shared causal architecture across populations
    
    **Post-hoc Combine** (`post_hoc_combine`)
    : Run single-input tools on each study, then combine results
    : Best for: Maximum flexibility and interpretability

### 4. Supported Fine-Mapping Tools

| Tool | Strategy | Best Use Case |
|------|----------|---------------|
| **SuSiE** | Single/Post-hoc | General purpose, robust |
| **FINEMAP** | Single/Post-hoc | Bayesian model averaging |
| **ABF** | Single/Post-hoc | Simple, fast Bayes factors |
| **CARMA** | Single/Post-hoc | Model uncertainty quantification |
| **RSparsePro** | Single/Post-hoc | Sparse regression approach |
| **MultiSuSiE** | Multi-input | Multi-population analysis |
| **SuSiEx** | Multi-input | Cross-ancestry fine-mapping |

## Key Features

### Quality Control Metrics

MAFM provides comprehensive QC including:

- **Consistency checks**: Kriging RSS to detect allele switches
- **LD structure**: Eigenvalue decomposition and 4th moment analysis  
- **Cross-study comparisons**: Cochran's Q test for heterogeneity
- **Frequency comparisons**: MAF consistency across studies

### Flexible Output

MAFM generates:

- **Credible sets**: Sets of variants likely to contain causal variants
- **Posterior inclusion probabilities (PIPs)**: Individual variant probabilities
- **QC reports**: Detailed quality control metrics
- **Meta-analysis results**: Combined summary statistics and LD matrices

## When to Use MAFM

!!! tip "Use Cases"
    
    **Multi-ancestry studies**
    : Leverage power across populations while accounting for LD differences
    
    **Multiple cohorts per ancestry**  
    : Combine studies within ancestry groups for increased power
    
    **Heterogeneous effect sizes**
    : Use post-hoc combination to preserve ancestry-specific signals
    
    **Quality control focus**
    : Extensive QC metrics help identify problematic data

!!! warning "Considerations"
    
    - Requires matched summary statistics and LD matrices
    - Computational requirements scale with number of studies
    - Some tools require specific data formats or parameter tuning

## Next Steps

Now that you understand the MAFM framework, let's jump into a practical example:

👉 **[Quick Start Guide](quick-start.md)** - Run your first MAFM analysis

Or explore specific scenarios:

- [Single-Input Fine-Mapping](single-input.md) for single-study analysis
- [Multi-Input Fine-Mapping](multi-input.md) for multi-study analysis 