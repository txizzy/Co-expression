# Description
Fast Pearson correlation calculation with multi-threading suitable for co-expression analysis.

# Install
1. git clone https://github.com/txizzy/Co-expression.git
2. cd Co-expression
3. cargo build -r

# Usage
1. Prepare your gene expression matrix. Use a table to split the columns. The first column should contain the gene names and a header is required.
2. target/release/Pearson -i infile1.xls -I infile2.xls -o output.xls -c 0.95 -p 0.05 -t 4

# Output
The output file contains four columns: gene name from infile1, gene name from infile2, correlation coefficient, and p-value, using table to split.


