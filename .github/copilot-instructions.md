# Pyteomics Copilot Instructions

## Project Overview
Pyteomics is a proteomics data analysis library focused on parsing mass spectrometry file formats and peptide/protein analysis. The architecture centers around format-specific parsers, including those built on a common XML parsing framework, as well as other specialized parsers for non-XML formats.

## Project Structure
The project follows a standard Python package layout:

### Main Code (`pyteomics/`)
- **Format parsers**: Individual modules for each supported format (e.g., `mzml.py`, `pepxml.py`, `mgf.py`, `fasta.py`)
- **Core modules**: `xml.py` (XML parsing framework), `parser.py` (peptide sequences), `proforma.py` (modern sequence notation)
- **Specialized subpackages**:
  - `auxiliary/`: Core utilities, base classes, and helper functions (`file_helpers.py`, `structures.py`, `utils.py`, `target_decoy.py`)
  - `mass/`: Mass calculations and UniMod integration (`mass.py`, `unimod.py`)
  - `openms/`: OpenMS-specific formats (`featurexml.py`, `idxml.py`, `trafoxml.py`)
- **Analysis modules**: `achrom.py` (chromatography), `electrochem.py` (electrochemistry), `pylab_aux.py` (plotting utilities)

### Tests (`tests/`)
- Test files: `test_*.py` for each module with corresponding test data files (`test.mzML`, `test.pepxml`, etc.)
- Controlled vocabularies: OBO files for CV validation (`psi-ms.obo`, `PSI-MOD.obo`)
- Test databases: `unimod.db` and `unimod.xml.gz` for UniMod testing

### Documentation (`doc/`)
- **Source**: RST files organized by topic (`data.rst`, `parser.rst`, `mass.rst`, etc.)
- **API reference**: Comprehensive module documentation in `source/api/`
- **Examples**: Practical usage examples in `source/examples/`
- **Build system**: Sphinx configuration and Makefile for HTML generation

## Core Architecture

### Unified Parser Hierarchy
All parsers inherit from common base classes in `pyteomics/auxiliary/file_helpers.py`:
- `FileReader`: Base class for all format parsers (XML and non-XML)
- `IndexedReaderMixin`: Adds indexing capabilities for random access
- `TaskMappingMixin`: Enables parallel processing via `map()` method

### XML Parser Hierarchy
XML-based format parsers additionally inherit from `pyteomics/xml.py`:
- `XML`: Base iterative parser with schema support
- `IndexedXML`: Adds byte-offset indexing for random access
- `MultiProcessingXML`: Combines indexing with parallel processing
- `IndexSavingXML`: Persistent index caching to disk

**XML Pattern**: Each XML format (mzML, pepXML, mzIdentML, etc.) combines some of these mixins:
```python
class MzML(BinaryArrayConversionMixin, CVParamParser, TimeOrderedIndexedReaderMixin, MultiProcessingXML, IndexSavingXML):
```

### Non-XML Parser Hierarchy
Non-XML formats (MGF, FASTA, MS1/MS2) use the same base infrastructure:
```python
class MGF(MGFBase, FileReader)                    # Sequential parser
class IndexedMGF(MGFBase, TaskMappingMixin, TimeOrderedIndexedReaderMixin, IndexSavingTextReader)  # Indexed parser
```

### Common Parser Interface
Every format parser provides:
- **Class-based**: `FormatClass(file, **kwargs)` - full-featured parser (iterative or indexed)
- **Functional**: `read(file, **kwargs)` - old functional interface, returns an instance of the parser class
- **Chain**: `chain(files, **kwargs)` - multi-file processing
- **Indexed access**: `parser[element_id]` when `use_index=True` or when directly instantiating an indexed parser

**Critical Parameters**:
- `iterative=True`: Memory-efficient streaming (default)
- `use_index=False`: Enable random access indexing
- `read_schema=False`: Parse XSD schemas for type conversion
- `retrieve_refs=True`: Auto-resolve ID references (mzIdentML)

## Key Modules

### `pyteomics/auxiliary/`
Core utilities and base classes:
- `file_helpers.py`: `FileReader`, indexing mixins, multiprocessing
- `structures.py`: `BasicComposition`, charge handling, unit types
- `utils.py`: Binary array conversion, base64 decoding
- `target_decoy.py`: FDR calculation algorithms

### `pyteomics/mass/`
- `mass.py`: Amino acid masses, isotope calculations
- `unimod.py`: UniMod database integration (requires SQLAlchemy)

### `pyteomics/parser.py`
Peptide sequence parsing using "modX" notation:
- `parse()`: "H-PEPTIDE-OH" → structured format
- `cleave()`: Enzymatic digestion simulation
- `amino_acid_composition()`: AA counting

### `pyteomics/proforma.py`
ProForma peptide sequence notation support:
- Modern standard for modified peptide sequences
- Controlled vocabulary integration
- Cross-conversion with modX format via `to_proforma()`

**Two Peptide Sequence Formats**:
- **modX**: Legacy format (e.g., "H-oxMPEPTIDE-OH") - simple, readable
- **ProForma**: Modern standard (e.g., "[Oxidation]M-PEPTIDE") - CV-based, precise

## Testing Strategy

### Test Structure
- Individual test files: `tests/test_*.py` for each module
- Test data: `tests/test.*` files (mzML, pepXML, etc.)
- Run individual tests: `python test_mzml.py` (from tests/ directory)
- CI runs: `find . -name 'test_*.py' -print0 | xargs -0 -n1 python`

### Test Data Dependencies
Tests expect specific test files in `tests/` directory. Use existing test data as examples for new format tests.

## Development Patterns

### Adding New XML Format Support
1. Create parser class inheriting appropriate XML mixins
2. Define `_default_iter_tag` and `_indexed_tags`
3. Implement `_get_info_smart()` for format-specific element processing
4. Add `read()`, `chain`, and `version_info` functions
5. Create comprehensive test with existing test data pattern

### Adding New Non-XML Format Support
1. Create base parser class inheriting from appropriate auxiliary mixins
2. Create sequential parser: `Format(FormatBase, FileReader)`
3. Create indexed parser: `IndexedFormat(FormatBase, TaskMappingMixin, IndexedTextReader)`
4. Implement `read()` function that dispatches based on `use_index` parameter
5. Follow same `read()`, `chain()`, test pattern as XML formats

### Binary Data Handling
Use `BinaryArrayConversionMixin` for formats with binary arrays (mzML, mzXML):
- Handles base64 decoding, compression (zlib, numpress)
- Automatic numpy array conversion
- `decode_binary=True` parameter controls this behavior

### Controlled Vocabularies
XML formats use psims library for CV validation:
- Enable OBO cache: `obo_cache.enabled = True` in tests
- Pass pre-created CV objects to avoid repeated downloads
- `CVParamParser` mixin handles CV parameter processing

## Build & Dependencies

### Core Dependencies
- `lxml`: XML parsing (required for most functionality)
- `numpy`: Numerical operations, binary arrays
- `psims`: Controlled vocabulary handling (XML formats)

### Optional Feature Dependencies
- `pandas`: DataFrame output (`DF` extra)
- `matplotlib`: Plotting utilities (`graphics` extra)
- `sqlalchemy>=1.4`: UniMod database (`Unimod` extra)
- `h5py`: mzMLb format support (`mzMLb` extra)

### Installation
```bash
pip install -e .  # Development install
pip install -e .[all]  # All optional dependencies
```

## Common Pitfalls

1. **Memory Usage**: Always use `iterative=True` for large files unless you need tree operations
2. **Index Performance**: Enable `use_index=True` for random access, but builds index on first access
3. **CV Downloads**: Enable psims caching to avoid repeated CV downloads in tests
4. **Binary Decoding**: Set `decode_binary=False` if you only need metadata, not spectra
5. **Schema Parsing**: `read_schema=True` enables automatic type conversion but can be slow

## File Format Conventions
- Module naming: `format.py` (e.g., `mzml.py`, `pepxml.py`)
- Class naming: `FormatXML` (e.g., `MzML`, `PepXML`)
- Test naming: `test_format.py` with `FormatTest` class
- Test data: `test.format` files in `tests/` directory