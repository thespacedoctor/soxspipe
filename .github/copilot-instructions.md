# SOXS Data Reduction Pipeline (soxspipe)

soxspipe is a Python data-reduction pipeline for the SOXS astronomical instrument, providing command-line tools for processing spectroscopic data through various calibration and science recipes.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

### Bootstrap and Environment Setup
- **ALWAYS use conda for installation** (recommended approach):
  ```bash
  source /usr/share/miniconda/etc/profile.d/conda.sh
  conda create -n soxspipe python=3.12 soxspipe -c conda-forge -y
  conda activate soxspipe
  ```
  - Takes 94 seconds to complete. NEVER CANCEL. Set timeout to 120+ minutes.
  - Conda activation requires sourcing the profile script first

- **Development installation** (from source):
  ```bash
  source /usr/share/miniconda/etc/profile.d/conda.sh
  conda activate soxspipe  # after creating conda environment above
  pip install -e .
  ```
  - Takes 5 seconds
  - Use this for making code changes

### Testing
- **Install test dependencies first**:
  ```bash
  pip install pytest nose2 pytest-profiling gprof2dot pyprof2calltree
  mkdir -p prof
  ```

- **Run quick individual tests**:
  ```bash
  pytest soxspipe/commonutils/tests/test_detector_lookup.py::test_detector_lookup::test_soxs_detector_lookup_function -v
  ```
  - Individual tests take ~2 seconds

- **Run test categories**:
  ```bash
  pytest -m "not full" -v    # Run non-full tests (lite tests)
  pytest -m "not slow" -v    # Skip slow tests  
  ```

- **WARNING: Full test suite limitations**:
  - Many tests require external test data from `/home/runner/xshooter-pipeline-data/` not available in repository
  - Some test files have import issues (missing `import pytest`)
  - Use individual test files for validation: `pytest soxspipe/commonutils/tests/test_detector_lookup.py -v`

### CLI Tool Usage and Validation
- **Verify installation**: `soxspipe -v` (should show version 0.13.4)
- **Basic help**: `soxspipe --help` 
- **Test workspace creation**:
  ```bash
  mkdir -p /tmp/test_workspace
  soxspipe prep /tmp/test_workspace
  # Expected: "There are no FITS files in this directory. Please add your data before running \`soxspipe prep\`"
  ```

### Documentation
- **Build documentation** (has known issues):
  ```bash
  cd docs
  pip install sphinx
  make html  # Will show extension errors but should not fail completely
  ```

## Validation Scenarios

**ALWAYS test these scenarios after making changes**:

1. **Environment activation and CLI availability**:
   ```bash
   source /usr/share/miniconda/etc/profile.d/conda.sh
   conda activate soxspipe
   soxspipe -v  # Should display version
   ```

2. **Basic CLI functionality**:
   ```bash
   soxspipe --help | head -30  # Should show usage information
   soxspipe mbias --help | head -10  # Should show recipe help
   ```

3. **Development installation test**:
   ```bash
   pip install -e .  # Should complete successfully
   soxspipe -v  # Should still work
   ```

4. **Basic test execution**:
   ```bash
   pytest soxspipe/commonutils/tests/test_detector_lookup.py::test_detector_lookup::test_soxs_detector_lookup_function -v
   # Should pass in ~2 seconds
   ```

## Common Tasks

### Install and Activate Environment
```bash
# Always source conda first
source /usr/share/miniconda/etc/profile.d/conda.sh

# Create environment (first time only) - NEVER CANCEL, 120+ second timeout
conda create -n soxspipe python=3.12 soxspipe -c conda-forge -y

# Activate environment  
conda activate soxspipe

# Verify installation
soxspipe -v
```

### Development Workflow
```bash
# After environment setup above
cd /path/to/soxspipe/repository
pip install -e .  # Install in development mode

# Install test dependencies
pip install pytest nose2 pytest-profiling gprof2dot pyprof2calltree
mkdir -p prof

# Run targeted tests
pytest soxspipe/commonutils/tests/test_detector_lookup.py -v

# Test CLI after changes
soxspipe -v
soxspipe --help
```

### Key CLI Commands and Recipes
- `soxspipe prep <workspaceDirectory>` - Prepare workspace for data reduction
- `soxspipe session ls` - List data reduction sessions
- `soxspipe mbias <inputFrames>` - Master bias recipe
- `soxspipe mdark <inputFrames>` - Master dark recipe  
- `soxspipe mflat <inputFrames>` - Master flat recipe
- `soxspipe disp_sol <inputFrames>` - Dispersion solution recipe
- `soxspipe order_centres <inputFrames>` - Order centers recipe
- `soxspipe stare <inputFrames>` - Process stare mode science frames
- `soxspipe nod <inputFrames>` - Process nodding mode science frames

## Repository Structure Reference

### Root Directory
```
.
├── README.md               # Installation and basic usage
├── setup.py               # Python package configuration
├── environment.yml        # Conda environment specification
├── Makefile              # Test targets (litetest, fulltest)
├── pytest.ini           # Test configuration with markers
├── nose2.cfg            # Alternative test runner config
├── Jenkinsfile          # CI configuration (Jenkins, not GitHub Actions)
├── docs/                # Sphinx documentation
├── soxspipe/           # Main package directory
│   ├── cl_utils.py     # Command-line interface
│   ├── recipes/        # Data reduction recipes
│   ├── commonutils/    # Shared utilities
│   └── tests/          # Unit tests
└── prof/               # Profiling output directory (create manually)
```

### Key Files
- `soxspipe/__version__.py` - Version information
- `soxspipe/cl_utils.py` - Main CLI entry point
- `setup.py` - Dependencies and package metadata
- `pytest.ini` - Test markers: 'slow', 'full', 'serial'

## Common Issues and Workarounds

1. **Conda activation fails**: Always run `source /usr/share/miniconda/etc/profile.d/conda.sh` first

2. **pkg_resources warnings**: Expected deprecation warnings from fundamentals package, can be ignored

3. **Test data missing**: Full test suite needs external data, use individual test files for validation

4. **Documentation build errors**: Known extension issues, focus on code functionality

5. **Import errors in tests**: Some test files have missing imports, skip problematic test files

6. **Makefile test targets fail**: Need pytest-profiling and prof/ directory setup

## Expected Timing (with 50% safety buffer)
- **Conda environment creation**: 94 seconds (use 120+ second timeout)
- **Pip development install**: 5 seconds (use 30+ second timeout)  
- **Individual tests**: 2 seconds (use 30+ second timeout)
- **CLI operations**: Immediate (use 10+ second timeout)
- **Documentation build attempt**: 2 seconds (use 30+ second timeout)

## CRITICAL Reminders
- **NEVER CANCEL** conda environment creation - it takes time but will complete
- **ALWAYS** source conda profile script before activation
- **DO NOT** expect full test suite to pass without external test data
- **ALWAYS** test CLI functionality after making changes to validate the pipeline works
- Use conda-forge channel for all conda operations
- Development changes require `pip install -e .` to take effect