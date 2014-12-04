BW_Mapper
=========

Our BW implementation for a Read mapper. Bio465

To Use:
python bbw_mapper.py input (prefix for .ref and .bwref files) read (full name of read file) -o output (optional output) 

Examples:
1) For 1data
python bbw_mapper.py 1data/chr1 1data/chr1.reads -o output1

2) For 2data
python bbw_mapper.py 2data/Chr19Partial 1data/Chr19Partial.bwref.reads -o output2
