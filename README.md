# Trees

## Raw data structure

In the `testdata` directory you’ll find the provided input files.  
Each ENSG number appears to correspond to a different tree — so far I’ve focused on `ENSG00000013016`.

The `.table.dat` file has two columns: the first lists the **parent node**, and the second lists the **child node**.  
The `.branchlength.dat` file contains the **branch length** for each edge.  
If a node is a terminal node (a tip of the tree), it also appears in the `.msa.dat` file with its **species name** and **sequence**.

## The Node and Tree objects

The `tree.py` script parses these files and defines two classes:  

- **`Node`** – represents a single node with its name, parent, children, branch length, and (for terminal nodes) the sequence.  
- **`Tree`** – initializes all nodes, stores them in a dictionary, and provides helper functions to find the root and the tips.

## Tree visualizer

The `visualize_tree.py` script contains a simple **ASCII tree visualizer**.  
It takes a `Tree` object and prints the structure in a readable tree format.  
I don’t expect to use it frequently, but it’s handy for quick inspection.

## Testing

In `tree.py`, the helper function **`parse_test_data()`** loads the example dataset from `testdata` and returns a ready-to-use `Tree` object.
