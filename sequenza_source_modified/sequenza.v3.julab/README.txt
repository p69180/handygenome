2020-03-17
Replace package "copynumber" -> "copynumbersypark"

2021-09-01
- init by pjh
- Wrote "sequenza.extract.modified" function (revised version of "sequenza.extract").
  - Entire *.seqz file is loaded onto a single tibble named "seqz.data.all".
  - Downstream functions using seqz file data (gc.sample.stats, read.seqz, etc.) utilizes the tibble.
  - Wrote "gc.sample.stats.modified" (which uses "seqz.data.all") function which is used instead of "gc.sample.stats".
  - Those measures result in less memory occupancy.
- Added 'sequenza.result2' and 'sequenza.result3' functions from "/home/users/sypark/01_Python_files/sequenza/sequenza_report_210127.R" into "R/results.R" script.

2021-09-02
- Wrote "windowValues.mod" and "windowBf.mod" functions. These are used in "sequenza.extract.modified".
  - Faster and uses less cpu than "windowValues" and "windowBf" respectively.
- Modified "sequenza.extract.modified" function.
  - "seqz.data" tibble is removed and garbage collected after every chromsome looping.
  - "seqz.data.all" tibble is removed and garbage collected after all chromosome looping has finished.