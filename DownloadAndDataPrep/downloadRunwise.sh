# Script that downloads files runwise from the GRID and merges them
# In principle, you only need to set arguments in the first NOTE section below, but also carefully check the initial terminal output before you continue running the script

# NOTE: Set input arguments
    downloadFolder="$1" # e.g. "vAN-20180727-1"
    trainNr="$2"        # e.g. "2445"
    periodName="$3"     # e.g. "LHC18d"
    # * DPG runlists *
        softwareDir=$(pwd)
        runlists=(
            DPGTracks
            DPGTracksAndCalo
            DPGTracksIncAcc
        )
    # * Other arguments
        startsWith="HeavyNeutralMesonToGG_"
        passNr="1"
        searchForFile="root_archive.zip"
        URLs="/alice/data/*/${periodName}/*/pass${passNr}/PWGGA/GA_pp/${trainNr}_201?????-*_child_*/${searchForFile}"
        checkFinalFormat="^HeavyNeutralMesonToGG_LHC1.*-DPG.*-pass${passNr}-.*.root\$" # REGEX, used to check if upward copy file name is correct
#

# * Set script parameters *
    sampleHit=$(alien_ls -b ${URLs} | head -n 1)
    year=$(echo ${sampleHit} | cut -d "/" -f 4)
    child=$(echo ${sampleHit} | cut -d "/" -f 10)
    child=${child/${child%child_*}}
    subfolder="${trainNr}-${child}"
#

# * Check input *
    echo -e "\nSearching for files of format:"
    echo -e "\"${URLs}\"\n"
    echo "Results in following parameters:"
    echo " - Train number: ${trainNr}"
    echo " - Year:         ${year}"
    echo " - Period:       ${periodName} (${child})"
    echo " - Pass number:  ${passNr}"
    echo " - Starts with:  \"${startsWith}\""
#

# * Make download directory and navigate to it *
    # Remove possible final slash from input argument
    folder="$(dirname $1)/$(basename $1)"
    # Make directory and navigate
    if [ ! -d "${folder}/${subfolder}" ]; then
        mkdir "${folder}/${subfolder}"
    fi
    cd "${folder}/${subfolder}"
#

# * Construct list of run numbers *
    echo -e "\nFound $(alien_ls ${URLs} | wc -l) \"${searchForFile}\" files:"
    echo "" > temp_GridFiles.txt # empty file
    for i in $(alien_ls -b ${URLs}); do
        if [[ ${i} =~ ^/.*/${searchForFile}+$ ]]; then
            echo "${i}"
            echo "${i}" >> temp_GridFiles.txt
        fi
    done
    echo ""
    read -p "Press ENTER to continue..."
    echo ""
#

# * Download files from GRID and extract *
    echo "" > temp_Modes.txt
    for i in $(cat temp_GridFiles.txt); do
        # * Get run number and make subdirectory in non-existent
        runNr=$(echo ${i} | cut -d "/" -f 6)
        if [ ! -d "${runNr}" ]; then
            mkdir "${runNr}"
        fi
        # * Copy file from GRID
        if [ -s "${runNr}/$(basename ${i})" ]; then
            echo "File \"${runNr}/$(basename ${i})\" already exists -- no need to copy from GRID"
        else
            echo -e "\nCopying run number \"${runNr}\" from GRID..."
            alien_cp alien:${i} file:${runNr}/
        fi
        # * Extract files
        unzip -uoqq "${runNr}/$(basename ${i})" -d "${runNr}"
        # * Make list of ROOT filenames (representing modes) to be merged
        for filename in $(ls ${runNr}/${startsWith}*); do
            filenameToAdd=$(basename ${filename})
            # * Check if name already exists in temp_Modes.txt
            for i in $(cat temp_Modes.txt); do
                filenameToAdd=${filenameToAdd/${i}}
            done
            # * Add mode (filename) if it does not exist yet
            if [ "${filenameToAdd}" != "" ]; then
                echo "$(basename ${filename})" >> temp_Modes.txt
            fi
        done
    done
#

# * Download files from GRID and extract *
    # * LOOP over DPG runlists
    for runlist in ${runlists[@]}; do
        # * Get list of DPG checked run numbers
        runlistFile="${softwareDir}/runlists/runNumbers${periodName}_${runlist}.txt"
        if [ -s "${runlistFile}" ]; then
            echo -e "\n\n---=== Merging according to runlist file \"$(basename ${runlistFile})\" ===---"
            # * LOOP over file names (modes)
            for mode in $(cat temp_Modes.txt); do
                # * Create empty file
                echo "" > temp_FilesToMerge.txt
                # * LOOP over run numbers
                for run in $(cat "${runlistFile}"); do
                    # * Add leading zeros to run number
                    runNr=$(printf "%09d\n" ${run})
                    # * Add file name to respective list files
                    if [ -s "${runNr}/${mode}" ]; then
                        echo "${runNr}/${mode}" >> temp_FilesToMerge.txt
                    # else
                    #     echo "\"${runNr}/${mode}\" does not exist!"
                    fi
                done # end of ${run} LOOP
                # * Generate output merge filename
                outputFile=${mode/${startsWith}/${startsWith}${periodName}-pass${passNr}-${runlist}_}
                # * Merge files
                echo ""
                hadd -f ${outputFile} $(cat temp_FilesToMerge.txt)
            done # end of ${mode} LOOP
        else
            echo "Runlist file \"$(basename ${runlistFile})\" does not exist!"
        fi
    done # end of ${runlist} LOOP
#

# * Remove temporary files *
    rm -f temp_FilesToMerge.txt
    rm -f temp_GridFiles.txt
    rm -f temp_Modes.txt
#

# * Navigate back to original working directory *
    echo -e "\nNavigating back to original working directory:"
    cd -
    echo ""
#