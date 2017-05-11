# auditoryPsychExp
Auditory Similarity and MDS Project

Psychophysics project analysing participants perceptions of the similarity/differences between sound tones

Still in the piloting stages with 3 different versions created thus far:
1) 12 frequencies only (no changes in the sound intensity)(250 Hz to 8kHz)
2) 12 frequencies and 6 intensities (250 Hz to 8 kHz, 11 to 61 dB)
3) 9 frequencies (250 to 500 Hz, linearly spaced) and 4 types of sound waves (sine, square, sawtooth, triangle)

Scripts here use tones created offline (may need to change this to create them in Matlab instead)

All responses by participants are presented on a 8 AFC wheel from -4 to 4 (this uses the response_screen.m script)

Auxilary scripts are provided to:
a) Patch sessions that are incomplete or where an error occurred (ensures results are concatenated over multiple sessions and that no same tone pairings are played)
b) Graphing scripts are provided to visualise the data grouped by either Hz, dB or wavetype
c) Other code is provided that helps with other ways of visualising/rearranging the data (randomCode_AuxExp.m)
