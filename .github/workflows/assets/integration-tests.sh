#!/bin/bash

set -e

ACCEPT_DEFAULTS=false
for arg in "$@"; do
    case $arg in
        -y|--yes) ACCEPT_DEFAULTS=true ;;
    esac
done

echo "## CHECKLIST BEFORE YOU START"
echo "1. Have you checked in all the code into the soxspipe repo?"
echo "2. Have you updated the default soxspipe.yaml file?"
echo "3. Have you updated the master soxspipe.db database?"

# ASK USER TO CHOOSE TO CONTINUE OR NOT
if $ACCEPT_DEFAULTS; then REPLY=Y; else read -p "Do you want to continue? ([Y]/n) " -r; fi
REPLY=${REPLY:-Y}
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    echo "Aborting."
    exit 1
fi
sleep 1

MINIMUM_TESTING="https://www.dropbox.com/scl/fo/6io4fkryrx13la2z9p33q/ALN6daExPQ8E0zWbBEgaM0U?rlkey=qw4i648kkw494poyo41ewkyrz&dl=1"
QUICKSTART_DATA_TESTING="https://www.dropbox.com/scl/fo/bpmh99yy0iyg0bzo05dqh/AD6eclegt7qpgLibqDRMmRQ?rlkey=w37onmwxkberl52k241iokw1k&st=5a77gb7l&dl=1"

# ASK USER TO CHOOSE WHICH DATASET TO USE
echo "Which dataset do you want to use for testing?"
echo "1. Minimum Testing Dataset"
echo "2. Quickstart Testing Dataset"
if $ACCEPT_DEFAULTS; then REPLY=1; else read -p "Enter the number of the dataset you want to use: [1]" -r; fi
REPLY=${REPLY:-1}
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[1]$ ]]
then
    DSNUM=1
    DROPBOX_LINK="$MINIMUM_TESTING"
elif [[ $REPLY =~ ^[2]$ ]]
then
    DSNUM=2
    DROPBOX_LINK="$QUICKSTART_DATA_TESTING"
else
    echo "Invalid input. Aborting."
    exit 1
fi




# Download test data
DOWNLOAD_DIR="/tmp/soxs_data_$DSNUM"
ZIP_FILE="/tmp/soxs_data_$DSNUM.zip"



download_and_unzip() {
    echo "Downloading $DOWNLOAD_DIR..."
    rm -rf "$ZIP_FILE"
    curl -L "$DROPBOX_LINK" -o "$ZIP_FILE"
    sleep 1
    unzip "$ZIP_FILE" -d "$DOWNLOAD_DIR" -x /
    sleep 1
    rm -rf "$ZIP_FILE"
    echo "Download complete."
}

reset_workspace() {
    echo "" > /tmp/soxspipe_prep.log
    sleep 3
    if [ -f "$HOME/anaconda/etc/profile.d/conda.sh" ]; then
        open /tmp/soxspipe_prep.log
        . "$HOME/anaconda/etc/profile.d/conda.sh"
        conda activate soxspipe
    fi
    cd "$DOWNLOAD_DIR"
    if [ -d "raw" ]; then
        mv raw /tmp
        rm -rf *
        mv /tmp/raw .
    fi
    soxspipe prep . >> /tmp/soxspipe_prep.log
    sleep 1
}


if [ -d "$DOWNLOAD_DIR" ]; then
    if $ACCEPT_DEFAULTS; then REPLY=N; else read -p "'$DOWNLOAD_DIR' already exists. Delete and redownload? (y/[N]) " -r; fi
    REPLY=${REPLY:-N}
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$DOWNLOAD_DIR"
        download_and_unzip
    else
        echo "Proceeding with existing '$DOWNLOAD_DIR'."
    fi
else
    echo "Downloading $DOWNLOAD_DIR..."
    download_and_unzip  
fi


echo "## RUNNING HELPER COMMANDS"
reset_workspace
soxspipe list sof  >> /tmp/soxspipe_list_sof.log
soxspipe list ob  >> /tmp/soxspipe_list_sof.log
soxspipe raw sof 20250514T040311_NIR_3_OFFSET_OBJ_SLIT1_0_300_0S_SOXS_CD-4412736.sof >> /tmp/soxspipe_list_sof.log

echo "## TESTING SINGLE SOF REDUCTIONS"
soxspipe -q reduce sof 20250514T040311_NIR_3_OFFSET_OBJ_SLIT1_0_300_0S_SOXS_CD-4412736.sof  >> /tmp/soxspipe_prep.log
echo "## TESTING BATCH REDUCTION"
soxspipe -q reduce all . -b 25 >> /tmp/soxspipe_prep.log
echo "## TESTING SINGLE-THREAD REDUCTION"
soxspipe -q reduce all . >> /tmp/soxspipe_prep.log

echo "## TESTING MULTIPROCESSING REDUCTIONS"
reset_workspace
soxspipe reduce all . -m >> /tmp/soxspipe_prep.log

# ASK USER TO CHOOSE TO CONTINUE OR NOT
if $ACCEPT_DEFAULTS; then REPLY=N; else read -p "Test daemon mode reductions. (y/[N]) " -r; fi
REPLY=${REPLY:-N}
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "## TESTING DAEMON MODE REDUCTIONS"
    reset_workspace
    soxspipe watch start  >> /tmp/soxspipe_prep.log

    code /tmp/soxspipe_prep.log
    if [ -f "$HOME/anaconda/etc/profile.d/conda.sh" ]; then
        open ~/.config/soxspipe/daemon.log
    fi
    
    sleep 8
    soxspipe watch status >> /tmp/soxspipe_prep.log
    sleep 8
    soxspipe watch stop >> /tmp/soxspipe_prep.log
fi


echo "## READY TO RELEASE!"

