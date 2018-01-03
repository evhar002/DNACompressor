README
=============
## DNACompressor: Fast and Deep Compression for Raw and Compressed Sequencing Data

### Download and Install: 
**Requirements:** <br>
`Python 2.7.*` <br>
`bz2file` <br>


```Bash
## Linux or Macox command line
(sudo) pip install bz2file
```
### Main arguments:
usage: DNACompressor <-f file> <-d> [-o output prefix] <br>

Example: <br>
`DNACompressor -f ./test/sample.fasta -o sample_out `<br>
`DNACompressor -f ./test/sample_out.bc.bz2 -d -o sample_2.fasta `<br>

### optional arguments: 
|  Parameter   |  Introduction |
| :---------- | :-------- |
|  -f INPUT_FILE, --file INPUT_FILE |   input filename
|  -o OUTPUT_FILE, --output OUTPUT_FILE |  Output filename prefix(Default:noname)
|  -d, --decode   |       extract from the compressed file
|  -l LEVEL, --level LEVEL |  Compression level. Default:6
