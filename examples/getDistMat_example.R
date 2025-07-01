# Example demonstrating the getDistMat function
# This example shows how to use getDistMat to access distance matrices
# from CoPro objects using the new flat structure

library(CoPro)

# Load or create a CoPro object (example assumes object already exists)
# object <- newCoProSingle(...)

# Example 1: Single slide object with one cell type (within-cell-type distances)
# After running computeDistance() on a CoProSingle object
# dist_matrix <- getDistMat(object, cellType1 = "TypeA", cellType2 = "TypeA")

# Example 2: Single slide object with two cell types (pairwise distances)
# dist_matrix <- getDistMat(object, cellType1 = "TypeA", cellType2 = "TypeB")

# Example 3: Multi-slide object requiring slide parameter
# dist_matrix <- getDistMat(object, cellType1 = "TypeA", cellType2 = "TypeB", slide = "slide1")

# Example 4: Force return of transpose
# dist_matrix_T <- getDistMat(object, cellType1 = "TypeA", cellType2 = "TypeB", 
#                            returnTranspose = TRUE)

# Example 5: Inspect available distance matrices
# Check the names of stored flat distance matrices
# names(object@distances)

# Example names in flat structure:
# Single slide: "dist|TypeA|TypeB", "dist|TypeA|TypeA"
# Multi slide: "dist|slide1|TypeA|TypeB", "dist|slide2|TypeA|TypeA"

# Example 6: Error handling - missing distance matrix
# This will give a helpful error message if the distance matrix doesn't exist:
# try({
#   dist_matrix <- getDistMat(object, cellType1 = "NonExistent", cellType2 = "TypeA")
# })

# Example 7: Symmetric access
# If you request TypeB -> TypeA but only TypeA -> TypeB exists,
# getDistMat will automatically return the transpose with a message
# dist_matrix <- getDistMat(object, cellType1 = "TypeB", cellType2 = "TypeA", verbose = TRUE)

print("getDistMat examples loaded successfully!")
print("Make sure to run computeDistance() before using getDistMat()") 