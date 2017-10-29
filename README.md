# Acoustic Similarity

We created custom software in MATLAB using many features created for Sound Analysis Pro (Tchernichovski, O. 2000 _Anim Behav._ 59:1167-1176) in order to 

* assign renditions of immature vocal utterances in juvenile songbirds to different syllable types, and 
* measure the acoustic similarity of these immature syllable renditions to mature syllables learned from an adult tutor.  

To assign syllable renditions to different types, we employed a combination of two measures of the acoustic distance between syllables that was then used to cluster syllables.  

1. The first distance measure was based on summary statistics of 10 different acoustic features.
2. The second distance measure was based on time-varying changes in 5 of those features.  

Similarity to tutor was calculated as the acoustic distance between each syllable rendition and its closest tutor syllable.  An alternate way of defining tutor similarity is as the distance between the center of the syllable cluster to which each rendition belonged and the closest tutor syllable (this is commented out in the code); for our data these two methods yielded highly similar results.

## Prerequisites

1. MATLAB R2013a or higher.

## Installation

1. Download the package by cloning this repository.
2. Open MATLAB, navigate to the download location, and add the root folder of the repository + all subfolders to the MATLAB path.

## Example Usage

You can segment syllables, perform clustering on syllables, or calculate tutor similarity using the provided functions. An example workflow is given below utilizing a sample data file that has been given.

1. Open MATLAB, and run segmentAndCluster.m, located in Acoustic_Similarity/workflow.
2. When prompted, select the example data file - Acoustic_Similarity/data/Gy242/Gy242_08_31_039.mat
3. As the script runs, check its progress and follow the onscreen prompts as they appear. 
   **Important:** Variables that the script prompts you to save must be saved in the same directory as the original data file (e.g. Acoustic_Similarity/data/Gy242).
4. (Instructions to run compareTutorJuvenile.m will go here. Need to check with John first to ensure all the inputs to the script have been generated.)

## Authors

Can fill out this section with John's info

## License

This project is licensed under the GNU GPLv3 License - see the [LICENSE.md](LICENSE) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

