# kinshipLR

## Description

An R package for evaluating the efficiency of autosomal markers in complex kinship testing between two individuals by likelihood ratio method. This package provides functions for calculating likelihood ratios for various relationship types including:

- Parent-Child relationship
- Full Sibling relationship
- Half-Sibling relationship
- Grandparent-Grandchild relationship
- Uncle-Nephew relationship
- Unrelated individuals

## Installation

You can install the development version of kinshipLR from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("OllieChueng/kinshipLR")
```

## Usage

```R
library(kinshipLR)

# Example frequency data
freq_data <- data.frame(
  Locus = c("A1", "A2", "A3"),
  D1S1656 = c(0.1, 0.2, 0.3),
  D2S441 = c(0.15, 0.25, 0.35)
)

# Calculate LR for Parent-Child relationship,
# n is simulation times
result <- LR_incl_PC(freq_data, n = 1000) 
```

## Functions

- `LR_incl_PC()`: Calculate likelihood ratio for Parent-Child relationship
- `LR_incl_FS()`: Calculate likelihood ratio for Full Sibling relationship
- `LR_incl_HS()`: Calculate likelihood ratio for Half-Sibling relationship
- `LR_incl_GP()`: Calculate likelihood ratio for Grandparent-Grandchild relationship
- `LR_incl_UA()`: Calculate likelihood ratio for Uncle-Nephew relationship
- `LR_excl_FK()`: Calculate likelihood ratio for unrelated individuals (first degree kinship)
- `LR_excl_SK()`: Calculate likelihood ratio for unrelated individuals (second degree kinship)
- `LR_excl_FS()`: Calculate likelihood ratio for unrelated individuals using FS method

## License

MIT
