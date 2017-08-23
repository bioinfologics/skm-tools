# Implementation notes

skipmers are (as of now) just stored into a uint64_t in 2-bit-per-base encoding.

All tools basically work over sorted collections of skip-mers, with different extra data attached.

These collections are dumped in gz'd files. The readers automatically detect if the input is compressed and decompress
on-the-fly.

# Todo

* Review aligner
* Read multiple files
* Multithread count
* Disk batches (implement merge)
* Test taxonomic classification