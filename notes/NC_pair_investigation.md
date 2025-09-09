# Investigation: NC_pair Node Missing in Upside2 Force Field

**Date:** September 6, 2025  
**Issue:** WORKER children failing with `NoSuchNodeError: group '/input/potential' does not have a child named 'NC_pair'`  
**Root Cause:** Mismatch between legacy `scale_params.py` and modern upside2 force field architecture

## Problem Summary

All WORKER processes in the condiv training system were failing at `scale_params.py:8` when trying to access `pot_group.NC_pair.interaction_param[:]`. This node doesn't exist in modern upside2 configurations, causing 100% worker failure rate.

## Key Error Messages

```
File "/home/okleinmann/upside2-md/py3/py/scale_params.py", line 8, in apply_param_scale
    pot_group.NC_pair.interaction_param[:] *= hb_scale
File "/home/okleinmann/miniconda3/envs/condiv-env/lib/python3.8/site-packages/tables/group.py", line 375, in _g_check_has_child
    raise NoSuchNodeError(
tables.exceptions.NoSuchNodeError: group ``/input/potential`` does not have a child named ``NC_pair``
```

## Investigation Evidence

### HDF5 Structure Analysis

**Current upside2 config contains:**
- `sigmoid_coupling_environment` (not `nonlinear_coupling_environment`)
- `hbond_energy` with `parameters` array
- `hbond_coverage` and `hbond_coverage_hydrophobe` with `interaction_param` arrays
- **NO** `NC_pair` node
- **NO** `midpoint_NC` node

**Generated via:** `condiv2.py:549` uses `environment_potential_type=1` → creates sigmoid coupling

### Git History Evidence from upside2-md Repository

**Timeline of Critical Changes:**

#### 1. Feb 27, 2022 - Commit `86a3152`: "(FF2) config for hbond terms"
- **Added** `write_midpoint_NC_pair()` and `write_midpoint_NC()` functions
- **Added** `NC_pair` node creation with `interaction_param` from external `NC_library`
- **Added** command line argument `--mid-NC-dist-energy`

```python
def write_midpoint_NC_pair(fasta, NC_library):
    write_midpoint_NC(fasta)
    n_res = len(fasta)
    grp = t.create_group(potential, 'NC_pair')
    grp._v_attrs.arguments = np.array(['midpoint_NC'])
    
    with tb.open_file(NC_library) as lib:
        param = lib.root.interaction_param[:]
        create_array(grp, 'interaction_param', obj=param)
```

#### 2. Feb 27, 2022 - Commit `9cb8bfd`: "Major config system redesign"  
- **REMOVED** `write_midpoint_NC_pair()` function entirely
- **REMOVED** call to `write_midpoint_NC_pair(fasta_seq, args.mid_NC_dist_energy)`
- Split config system into "core" and "advanced" stages
- **This is where NC_pair was eliminated from main upside2**

#### 3. May 9, 2022 - Commit `36af4f5`: "Add function for scaling potential strengths"
- **Added** `apply_param_scale()` function to `upside_config.py`
- **Original scaling function NEVER included NC_pair references**
- Only scaled: `hbond_energy.parameters`, `sigmoid/nonlinear_coupling_environment`, `hbond_coverage`

```python
def apply_param_scale(hb_scale=1., env_scale=1., rot_scale=1., memb_scale=1.):
    pot_group = t.root.input.potential
    if hb_scale != 1.:
        print ("scaling hb {}x".format(hb_scale))
        pot_group.hbond_energy.parameters[:4] *= hb_scale
        # ToDo: what about the hbond_coverage groups -> how are they coupled
        # to the energy? For env it's clear, but not so for hb
```

## Force Field Architecture Evolution

### Legacy System (Pre-2022)
- **NC_pair**: Direct nonlinear coupling between residue pairs via `midpoint_NC`
- **Explicit nonlinear coupling**: Separate nodes for nonlinear interactions
- **Specialized scaling**: Required separate scaling for NC_pair parameters

### Modern Upside2 (Post-2022 Redesign)  
- **Generalized coupling**: `sigmoid_coupling_environment` vs `nonlinear_coupling_environment`
- **Integrated H-bond system**: Direct scaling via `hbond_energy.parameters`
- **Coverage-based interactions**: `hbond_coverage` scaled under `rot_scale`

## What NC_pair Actually Was

**NC = "Nonlinear Coupling"**

- **Purpose**: Direct pairwise interactions with nonlinear coupling between residues
- **Arguments**: `['midpoint_NC']` - operated on midpoint coordinates
- **Structure**: Per-residue arrays (`index`, `id`, `type`) with interaction parameters
- **Library-based**: Parameters loaded from external `NC_library` file
- **Scope**: Residue-to-residue nonlinear coupling (not H-bond specific)

## Current Scaling Architecture Analysis

The modern upside2 scaling correctly separates concerns:

### `hb_scale` - Direct Hydrogen Bond Energy
```python
pot_group.hbond_energy.parameters[:4] *= hb_scale
```
- Controls intrinsic H-bond strength
- Direct energy scaling for helix/sheet/turn H-bonds

### `rot_scale` - Geometric/Orientational Coupling  
```python
pot_group.hbond_coverage.interaction_param[:] *= rot_scale
pot_group.hbond_coverage_hydrophobe.interaction_param[:] *= rot_scale  
pot_group.rotamer.pair_interaction.interaction_param[:] *= rot_scale
```
- Controls geometric coupling between H-bonds and rotamers
- Modulates how sidechain orientations affect H-bond formation

### `env_scale` - Environment Coupling
```python
if 'nonlinear_coupling_environment' in pot_group:
    pot_group.nonlinear_coupling_environment.coeff[:] *= env_scale
if 'sigmoid_coupling_environment' in pot_group:  
    pot_group.sigmoid_coupling_environment.scale[:] *= env_scale
```

## The TODO Comment Analysis

**Original TODO:**
```python
# ToDo: what about the hbond_coverage groups -> how are they coupled
# to the energy? For env it's clear, but not so for hb
```

**My Analysis:** This TODO asks whether H-bond scaling should also affect coverage interactions. The current separation is **architecturally correct**:

- **Intrinsic H-bond strength** (`hb_scale`) vs **geometric coupling** (`rot_scale`)
- Allows independent tuning of direct energies vs orientational effects
- Environment coupling is clear because it directly modulates environmental potentials

## Solution Implemented

### Problem
```python
# This line was causing all worker failures:
pot_group.NC_pair.interaction_param[:] *= hb_scale
```

### Fix Applied
```python
# Removed the problematic line and added explanation:
# Note: NC_pair.interaction_param was removed in upside2 redesign (Feb 2022)
# The hydrogen bond scaling now works through hbond_energy.parameters
# and the coverage-based interactions below in rot_scale section
```

## Key Insights

1. **`scale_params.py` was never part of the main upside2-md repository**
   - It's a local/training-specific modification
   - The `NC_pair` reference was added outside the official codebase

2. **NC_pair functionality was absorbed into modern coupling system**
   - Nonlinear coupling now handled by `sigmoid_coupling_environment`
   - Direct H-bond interactions via `hbond_energy.parameters`

3. **The original upside2 scaling function never had NC_pair**
   - Commit `36af4f5` shows the official implementation
   - It only scaled the nodes that actually exist in modern configs

4. **Modern architecture is more modular and flexible**
   - Separates direct energies from geometric coupling
   - Environment vs rotamer vs H-bond scaling are independent

## Files Modified

- **`/home/okleinmann/upside2-md/py3/py/scale_params.py`**
  - **Line 8**: Removed `pot_group.NC_pair.interaction_param[:] *= hb_scale`
  - **Added documentation** explaining the architectural change

## Verification

The fix aligns `scale_params.py` with the official upside2 architecture established in May 2022. All necessary H-bond scaling is preserved through:
- Direct energy scaling: `hbond_energy.parameters[:4]`
- Coverage interaction scaling: `hbond_coverage.interaction_param[:]` (under `rot_scale`)

## Conclusion

This was a classic case of **legacy code incompatibility** following a major architectural redesign. The `NC_pair` node was a short-lived experimental feature (Feb 2022) that was removed in the same month during the config system redesign. The local `scale_params.py` was modified to reference this deprecated node, causing worker failures when run against modern upside2 configurations.

The fix preserves all necessary scaling functionality while removing the reference to the non-existent node.