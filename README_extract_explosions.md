# extract_explosions
Acoustic characteristic extractor for explosions collected by HARPS (High-frequency acoustic recording packages)

Edit 2:
- dcOffset implemented for the entire execution of the program/analysis
- *testing* SNR calculated using the average RL of the signal
- Noise recieved level calculation added, along with a histogram. Includes TF 
Edit 1:
- Correct number of samples loaded and analyzed
- Correct smoothing 
- Duration
- Histogram creation for each of the acoustic characteristics at the end of execution
- Identification of start and end times as pulled from the detector
- Acoustic characteristics determined:
    - peak frequency
    - center frequency
    - 10 and 3 dB bandwidth
    - SEL
    - SPLrms
    - inter explosion interval
    - received level
    - duration
    - signal to noise ratio
