# Strategy for Eliminating Complex Nested Data Structures in CoPro

## Current Problems

The CoPro package currently uses deeply nested list structures that are:
- **Error-prone**: Easy to make indexing mistakes
- **Hard to debug**: Difficult to trace issues in nested access
- **Memory inefficient**: Redundant storage and poor cache locality
- **Maintenance burden**: Complex logic for single vs multi-slide cases

### Examples of Current Complexity:
```r
# 4-level nesting!
object@cellScores[["sigma_0.1"]][["slide1"]][["TypeA"]][, 1]

# Complex kernel access with symmetric fallback
K_list[[ct1]][[ct2]] or t(K_list[[ct2]][[ct1]])

# Nested distance structures  
object@distances[["slide1"]][["TypeA"]][["TypeB"]]
```

## Recommended Solutions (Ranked by Impact)

### 🥇 **Solution 1: Accessor Functions Strategy** (Immediate Impact)

**Why Best**: Can be implemented incrementally without breaking existing code.

```r
# BEFORE: object@cellScores[["sigma_0.1"]][["slide1"]][["TypeA"]]
# AFTER:  getCellScores(object, sigma = 0.1, cellType = "TypeA", slide = "slide1")

# BEFORE: K_list[[ct1]][[ct2]] with manual symmetric fallback  
# AFTER:  getKernelMatrix(object, sigma = 0.1, cellType1 = ct1, cellType2 = ct2)
```

**Benefits**:
- ✅ Hide complexity behind clean API
- ✅ Centralized validation and error handling
- ✅ Automatic symmetric matrix handling
- ✅ Can be added alongside existing code

**Implementation Steps**:
1. Create generic accessor functions
2. Add method dispatch for Single vs Multi objects
3. Gradually replace direct slot access in functions
4. Eventually deprecate direct access

### 🥈 **Solution 2: Flat Structure with Composite Keys** (Medium-term)

**Why Good**: Eliminates nesting while maintaining simple R data structures.

```r
# BEFORE: cellScores[["sigma_0.1"]][["slide1"]][["TypeA"]]
# AFTER:  cellScores_flat[["sigma_0.1_slide_slide1_celltype_TypeA"]]

# Helper function for clean access
getCellScoresFlat(cellScores_flat, sigma = 0.1, cellType = "TypeA", slide = "slide1")
```

**Benefits**:
- ✅ Single level of indexing
- ✅ Easier to debug and iterate over
- ✅ More memory efficient
- ✅ Works with existing R infrastructure

### 🥉 **Solution 3: Data.frame Registry + Storage** (Long-term)

**Why Powerful**: Most structured approach with metadata management.

```r
# Registry tracks what's stored where
registry <- data.frame(
  sigma = c(0.1, 0.1, 0.2),
  slide = c("slide1", "slide2", "slide1"), 
  cellType = c("TypeA", "TypeA", "TypeB"),
  matrixID = c("scores_1", "scores_2", "scores_3")
)

# Flat storage for actual matrices
storage@matrices[["scores_1"]] <- actual_matrix
```

**Benefits**:
- ✅ Most structured and queryable
- ✅ Built-in metadata tracking
- ✅ Easy to implement caching/lazy loading
- ✅ Natural fit for database-like operations

## Implementation Roadmap

### Phase 1: Immediate (Accessor Functions)
1. ✅ Create `getCellScores()` and `getKernelMatrix()` functions
2. Add method dispatch for CoProSingle vs CoProMulti
3. Replace direct access in plotting functions first
4. Add comprehensive input validation

### Phase 2: Medium-term (Hybrid Approach)
1. Implement flat storage as internal optimization
2. Keep accessor functions as external API
3. Gradually migrate internal functions to flat storage
4. Add performance benchmarks

### Phase 3: Long-term (Full Restructure)
1. Implement data.frame registry system
2. Add lazy loading for large matrices
3. Consider database backend for very large datasets
4. Migrate all slot structures to new system

## Specific Code Changes Needed

### 1. Update Class Definitions
```r
# Add accessor-friendly slots
setClass("CoPro", 
  slots = list(
    # ... existing slots ...
    .cellScoresFlat = "list",        # Internal flat storage
    .kernelMatricesFlat = "list",    # Internal flat storage  
    .scoreRegistry = "data.frame"    # Optional: registry system
  )
)
```

### 2. Replace Direct Access Throughout Codebase
```r
# Files to update (priority order):
# 1. R/80_get_cs_in_situ.r - main user-facing function
# 2. R/70_get_norm_corr_plot.R - plotting functions
# 3. R/16_gene_and_cell_score.R - core computation
# 4. examples/plotting_functions.R - user examples
```

### 3. Add Validation and Error Handling
```r
.validateInputs <- function(object, sigma, cellType, slide = NULL) {
  # Centralized validation logic
  # Better error messages
  # Suggest corrections when possible
}
```

## Benefits Summary

| Approach | Implementation | Maintenance | Performance | User Experience |
|----------|---------------|-------------|-------------|-----------------|
| Accessor Functions | Easy | Easy | Same | Much Better |
| Flat Keys | Medium | Easy | Better | Better |
| Registry System | Hard | Medium | Best | Best |

## Recommended Next Steps

1. **Start with Accessor Functions** - Implement `getCellScores()` and `getKernelMatrix()`
2. **Update User-Facing Functions** - Replace direct access in plotting and analysis functions  
3. **Add Comprehensive Tests** - Ensure accessor functions work for all edge cases
4. **Gradual Migration** - Replace direct slot access throughout codebase
5. **Performance Testing** - Benchmark to ensure no performance regression

This strategy provides a clear path from the current complex nested structures to a cleaner, more maintainable system while allowing incremental implementation. 