# kdatapack

Kirk's set of python programs and libraries for data analysis and sequence tools

## dbminers

_routines to extract protein or DNA sequences, focused on glycoside hydrolase applications_

* dbminers/pullycazy.py - scrape CAZY website for a given GH family, extracting metadata then calling Entrez e-utilities (via Biopython) to obtain fasta sequence and outputting results to a subfolder
  * for example, `python pullcazy.py 11 <youremailaddress>` will download of GH family 11, type `python pullycazy.py --help` to see options
* 