#!/bin/bash

# ==============================================================================
# Automated run for: Hartree, Hartree+Exchange, Hartree+Exchange+Correlation
# ==============================================================================

# Usage: ./run_dft_suit.sh [python_executable]
# If no argument is provided, it defaults to 'python'.
PYTHON_EXE="${1:-python}"

CONFIG_FILE="config.yaml"
BACKUP_FILE="config.yaml.bak"
PYTHON_DFT_PROGRAM="DFT_helium.py"
PYTHON_PLOT_PROGRAM="visualization/plot_results.py"

# Backup of the original config
cp "$CONFIG_FILE" "$BACKUP_FILE"

echo "" > run.log

# Function to update config and run simulation
run_simulation() {
    local folder_name=$1
    local use_exchange=$2
    local use_correlation=$3
    local img_spec_name=$4

    echo "----------------------------------------------------------------"
    echo "Running Simulation: $folder_name"
    echo " > Exchange:    $use_exchange"
    echo " > Correlation: $use_correlation"
    
    # Modify config.yaml    
    sed -i "s|out_dir: .*|out_dir: \"$folder_name\"|" "$CONFIG_FILE"
    sed -i "s|use_exchange: .*|use_exchange: $use_exchange|" "$CONFIG_FILE"
    sed -i "s|use_correlation: .*|use_correlation: $use_correlation|" "$CONFIG_FILE"

    # Run the Python program
    "$PYTHON_EXE" -X utf8 "$PYTHON_DFT_PROGRAM" >> run.log

    # Run also the plotting script to generate the plots
    "$PYTHON_EXE" -X utf8 "$PYTHON_PLOT_PROGRAM" "$folder_name" "$img_spec_name" "$folder_name" > /dev/null
    
    echo " > Done."
}

# Case 1: Hartree Only (No Exchange, No Correlation)
run_simulation "outputs/01_Hartree" "false" "false" "Hartree"

# Case 2: Hartree + Exchange (No Correlation)
run_simulation "outputs/02_Hartree_Exchange" "true" "false" "Hartree_Exchange"

# Case 3: Hartree + Exchange + Correlation (Full DFT-LDA)
run_simulation "outputs/03_Hartree_Exchange_Correlation" "true" "true" "Hartree_Exchange_Correlation"

# CLEANUP

# Restore the original configuration
mv "$BACKUP_FILE" "$CONFIG_FILE"
echo "----------------------------------------------------------------"
echo "Suite completed. Original config.yaml restored."