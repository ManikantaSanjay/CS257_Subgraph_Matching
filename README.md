
# Subgraph Matching
## Introduction
This project analyzes eight leading in-memory subgraph matching algorithms, including QuickSI, GraphQL, CFL, CECI, DP-iso, RI, and VF2++. We evaluate them based on:
1. Vertex filtering methods in data graphs.
2. Techniques for ordering query vertices.
3. Strategies for enumerating partial results.
4. Other optimization approaches.

We also compare these algorithms against Glasgow, a constraint programming-based algorithm. Our findings show the strengths of different methods in various aspects, leading to recommendations for specific scenarios based on graph density and query size.

## Compilation
To compile the source code, navigate to the project's root directory and execute:

```zsh
mkdir build
cd build
cmake ..
make
```

## Testing
To verify the binary's correctness, run the following:

```zsh
cd test
python test.py ../build/matching/SubgraphMatching.out
```

## Execution
Find the compiled 'SubgraphMatching.out' in 'build/matching'. Execute it using:

```zsh
./SubgraphMatching.out -d [data_graphs] -q [query_graphs] -filter [filter_method] -order [order_method] -engine [enum_method] -num [embeddings_number]
```

Set `-num` to 'MAX' for all results. The program stops upon reaching the limit or finding all results.

Example:

```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter GQL -order GQL -engine LFTJ -num MAX
```

## Input Format
The input graphs should be vertex-labeled, starting with 't N M' (N = vertices, M = edges). Format vertices as 'v VertexID LabelId Degree' and edges as 'e VertexId VertexId'. Vertex ids start from 0.

Example:

```zsh
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

## Supported Techniques
### Filtering Methods
- LDF: Label and degree filter (default).
- NLF: Neighborhood label frequency filter.
- GQL: GraphQL's method.
- TSO: TurboIso's method.
- CFL: CFL's method.
- DPiso: DP-iso's method.
- CECI: CECI's method.

### Ordering Methods
- QSI: QuickSI's method.
- GQL: GraphQL's method.
- TSO: TurboIso's method.
- CFL: CFL's method.
- DPiso: DP-iso's method.
- CECI: CECI's method.
- RI: RI's method.
- VF2++: VF2++'s method.

### Enumeration Methods
- QSI: QuickSI and RI's method.
- GQL: GraphQL's method.
- EXPLORE: CFL's method.
- DPiso: DP-iso's method.
- CECI: CECI's method.
- LFTJ: Local candidates computation.

Examples of command-line execution are provided for various configurations.

## Recommendations
We suggest using GQL, CFL, or DPiso for filtering; GQL or RI for ordering; and always LFTJ for enumeration. Enable failing set pruning for large queries and QFilter for dense graphs.

## Command-Line Execution Examples

### Basic Execution
Run a basic query using GraphQL's filtering and ordering methods with local candidate computation for enumeration:

```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter GQL -order GQL -engine LFTJ -num MAX
```

### Custom Execution
Execute with custom settings, such as using label and degree filter, CFL's ordering method, and set-intersection based local candidates computation:

```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter LDF -order CFL -engine LFTJ -num MAX
```

### Advanced Configuration
For complex queries, configure the set intersection algorithms and optimization techniques in 'configuration/config.h'.

### Spectrum Analysis
Perform spectrum analysis on a query with time limits:

```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter LDF -order Spectrum -engine Spectrum -order_num 100 -time_limit 60 -num 100000
```
