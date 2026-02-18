# ElPGx Documentation Index

## Documentation Overview

| Document | Description | Audience |
|----------|-------------|----------|
| [COMPLETE_DOCUMENTATION.md](COMPLETE_DOCUMENTATION.md) | Full technical documentation | Developers, Bioinformaticians |
| [Technical_Rationale.md](Technical_Rationale.md) | Design decisions and architecture | Developers |
| [PGx_Pipeline_Presentation.md](PGx_Pipeline_Presentation.md) | Overview presentation | All users |
| [FUTURE_DEVELOPMENT.md](FUTURE_DEVELOPMENT.md) | Roadmap and planned features | Contributors, Developers |
| [VALIDATION.md](VALIDATION.md) | Validation framework and metrics | QA, Clinical |
| [SAFEGUARDS_STATUS.md](SAFEGUARDS_STATUS.md) | Clinical safeguards implementation status | Clinical, Regulatory |

---

## Quick Links

### For Users
- [Installation Guide](COMPLETE_DOCUMENTATION.md#installation)
- [Running the Pipeline](COMPLETE_DOCUMENTATION.md#usage)
- [Interpreting Results](COMPLETE_DOCUMENTATION.md#output-format)

### For Developers
- [Architecture Overview](Technical_Rationale.md)
- [Contributing Guide](FUTURE_DEVELOPMENT.md#contributing)
- [Testing Framework](FUTURE_DEVELOPMENT.md#validation-and-quality-assurance)

### For Contributors
- [Priority Matrix](FUTURE_DEVELOPMENT.md#implementation-priority-matrix)
- [Confidence Score Enhancements](FUTURE_DEVELOPMENT.md#confidence-score-enhancements)
- [CNV Detection Plans](FUTURE_DEVELOPMENT.md#copy-number-variation-detection)

---

## Current Version: 1.1.0

### What's Included
- 18 pharmacogenes supported
- CPIC v1.24 official data
- SNV/Indel variant calling
- Star allele diplotype assignment
- Phenotype translation
- Clinical prescribing implications

### Known Limitations
- No CNV detection (affects CYP2D6 accuracy)
- Limited haplotype phasing
- No structural variant support

See [FUTURE_DEVELOPMENT.md](FUTURE_DEVELOPMENT.md) for planned enhancements.

---

## Data Sources

| Source | Version | Usage |
|--------|---------|-------|
| CPIC | v1.24 | Allele definitions, phenotypes |
| PharmVar | Current | Star allele nomenclature |
| GRCh38 | hg38 | Reference genome coordinates |

Data provenance: `data/PROVENANCE.json`
