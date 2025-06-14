# CATChem

**CATChem** (Configurable ATmospheric Chemistry) is a modelling component that includes all chemical and aerosol processes needed to perform atmospheric chemistry and composition simulations within a model through a flexible, easy to modify, and well-documented infrastructure.

## Features

- **NUOPC-compliant interface** for integration with Earth system models
- **ESMF I/O integration** for efficient parallel file operations
- **Flexible configuration system** with YAML-based setup
- **Automatic regridding** with weight file management
- **CF-compliant input/output** supporting Climate and Forecast conventions
- **Comprehensive utility tools** for configuration and weight generation

## Documentation

Comprehensive documentation is available in the `docs/` directory:

- **[NUOPC Integration Overview](docs/nuopc-overview.rst)** - Architecture and capabilities
- **[Examples and Tutorials](docs/nuopc-examples.rst)** - Step-by-step guides and common configurations
- **[Utility Scripts](docs/utilities.rst)** - Weight generation and configuration tools
- **[Configuration Reference](docs/nuopc-configuration.rst)** - Complete YAML configuration guide

## Quick Start

### NUOPC Interface

```bash
# Build the NUOPC interface
cd drivers/nuopc
make

# Set up utilities and weight files
../../util/generate_esmf_weights.sh setup

# Validate configuration
python ../../util/validate_weight_config.py catchem_input_config.yml

# Run standalone test
./catchem_nuopc_driver
```

### Utility Scripts

The `util/` directory contains helpful scripts:

- **`generate_esmf_weights.sh`** - ESMF weight file generation and management
- **`manage_weights.sh`** - Weight file optimization and maintenance
- **`validate_weight_config.py`** - Configuration validation and checking

## Development

### pre-commit

If you don't have the command-line tool `pre-commit` available, [install it](https://pre-commit.com/#install).

Install the pre-commit hooks with

```
pre-commit install --install-hooks
```

Now, some checks and auto-formatting will run automatically when you commit.

> [!NOTE]
> For source files with their own formatting that we don't intend to modify (or only modify slightly),
> e.g. "vendored" modules, we generally don't want `findent` to be applied.
> For such files, update the `exclude` section in `.pre-commit-config.yaml` accordingly.

> [!TIP]
> ```
> pre-commit run --all-files
> ```
> can be used to check the pre-commit config
> (e.g. to make sure that the `findent` `exclude` section is working as you intended)
> or to catch up if commits were made without pre-commit or the pre-commit config was updated.

### Build

To test the build locally, first configure, using `FC` to specify the compiler you want to use, then build. For example:

```
FC=gfortran-12 cmake -B build
```
```
cmake --build build -j
```

To clean, you can use

```
cmake --build build --target clean
```

or remove the build directory (`./build`).

### Test

To run the tests, after building, use

```
ctest --test-dir build/tests
```

There [are options](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html#testing-using-ctest) for selecting specific tests.

Edit `tests/CMakelists.txt` to add new tests.
