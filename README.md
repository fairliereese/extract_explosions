# extract_explosions_extensible
Extensible and less complicated version of extract explosions

Can handle a variety of data inputs; any permutation of the following:
- have detection file
- need to create detection file from Tethys database
- xwavs on one drive
- xwavs on two drives
- xwavs on three drives
- phantom detections coming in from Tethys will be disregarded

Other:
- includes a test case to see how code is called
- splitALL, splitALL_3 responsible for splitting detection file up along temporal boundaries generated from xwav files
- expCharacteristics called for every xwav drive there is

