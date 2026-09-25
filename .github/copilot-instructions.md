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
- **Install and run the required offline suite**:
  ```bash
  python -m pip install -e ".[tests]"
  python -m pytest tests/unit tests/integration -m "not slow"
  ```
- The required suite is configured in `pyproject.toml`, uses only synthetic local data, blocks network access, and must not write to package directories.
- Tests in `tests/real_data` are opt-in. They require a verified, prepared workspace supplied through `SOXSPIPE_REAL_DATA_DIR`; the `Real-data tests` GitHub Actions workflow creates it from the repository-owned manifest.

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
  python -m pip install -e ".[tests]"
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
├── pyproject.toml       # Package, pytest, and coverage configuration
├── .github/workflows/   # Required and opt-in GitHub Actions workflows
├── docs/                # Sphinx documentation
├── soxspipe/           # Main package directory
│   ├── cl_utils.py     # Command-line interface
│   ├── recipes/        # Data reduction recipes
│   ├── commonutils/    # Shared utilities
└── tests/               # Unit, integration, and opt-in real-data tests
```

### Key Files
- `soxspipe/__version__.py` - Version information
- `soxspipe/cl_utils.py` - Main CLI entry point
- `setup.py` - Dependencies and package metadata
- `pyproject.toml` - Test markers: `unit`, `integration`, `real_data`, `slow`, and `serial`

## Common Issues and Workarounds

1. **Conda activation fails**: Always run `source /usr/share/miniconda/etc/profile.d/conda.sh` first

2. **pkg_resources warnings**: Expected deprecation warnings from fundamentals package, can be ignored

3. **Real-data test data missing**: Run the required synthetic suite, or use the opt-in real-data workflow to create a verified workspace before running `tests/real_data`.

4. **Documentation build errors**: Known extension issues, focus on code functionality

## Expected Timing (with 50% safety buffer)
- **Conda environment creation**: 94 seconds (use 120+ second timeout)
- **Pip development install**: 5 seconds (use 30+ second timeout)  
- **Individual tests**: 2 seconds (use 30+ second timeout)
- **CLI operations**: Immediate (use 10+ second timeout)
- **Documentation build attempt**: 2 seconds (use 30+ second timeout)

## CRITICAL Reminders
- **NEVER CANCEL** conda environment creation - it takes time but will complete
- **ALWAYS** source conda profile script before activation
- **DO NOT** run `tests/real_data` without a verified workspace supplied through `SOXSPIPE_REAL_DATA_DIR`
- **ALWAYS** test CLI functionality after making changes to validate the pipeline works
- Use conda-forge channel for all conda operations
- Development changes require `pip install -e .` to take effect
