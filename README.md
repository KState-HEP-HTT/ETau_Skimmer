# Information about the skimmer

## Choosing which samples to skim
The skimming executable requires a nickname to be passed as a command-line argument. This nickname maps to the path to the directory containing the root files you wish to skim. All *.root files in the given path will be skimmed. Users may add new samples to util.h to process samples not already listed. The skimmer must be recompiled after the addition.

    {"sample_nickname", "path/to/directory"},

## Compiling the skimmer
The skimmer can be compiled using the build script. The build script gives the user the ability to name the output executable so that multiple versions of the skimmer can be executed. This name must be included as a command-line option.

    ./build <exe_name>

## Running the skimmer 
The skimmer is executed just like any other executable. The only difference is that the nickname must be provided to choose which samples to skim. For example,

    ./Skim ZZ

will take an executable named Skim and process all root files in the directory mapped by ZZ.

# To do
1. Make sure the skim selection is correct (using etau baseline for now)
2. Make sure all necessary branches are being kept
3. Remove all unnecessary branches
4. If others are to use this skimmer, the path to store output files must be made more customizable
