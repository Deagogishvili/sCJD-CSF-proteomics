#!/usr/bin/env bash

# Directory containing the hierarchical-hotnet folder
BASE_DIR="$PWD"
HOTNET_DIR="$BASE_DIR/hierarchical-hotnet"

# Array of module directories
MODULES=(
    "01_blue"
    "02_brown"
    "03_green"
    "04_grey"
    "05_red"
    "06_turquoise"
    "07_yellow"
)

# Function to run a module
run_module() {
    local module=$1
    echo "================================================"
    echo "Starting analysis for module: $module"
    echo "================================================"
    
    # Change to module directory
    cd "$HOTNET_DIR/$module"
    
    # Check if directory exists
    if [ ! -d "$HOTNET_DIR/$module" ]; then
        echo "Error: Module directory $module not found!"
        return 1
    fi
    
    # Check if run script exists
    if [ ! -f "run_commands_parallel.sh" ]; then
        echo "Error: run_commands_parallel.sh not found in $module!"
        return 1
    fi
    
    # Make sure the script is executable
    chmod +x run_commands_parallel.sh
    
    # Run the module's script
    ./run_commands_parallel.sh
    
    # Check if the script executed successfully
    if [ $? -eq 0 ]; then
        echo "Successfully completed analysis for $module"
    else
        echo "Error occurred while running $module"
        return 1
    fi
    
    echo ""
}

# Check if we're in the correct directory
if [ ! -d "$HOTNET_DIR" ]; then
    echo "Error: hierarchical-hotnet directory not found!"
    echo "Please make sure you're running this script from the 10_hierarchical_hotnet directory"
    exit 1
fi

# Main execution
echo "Starting HotNet analysis for all modules..."
echo "Total modules to process: ${#MODULES[@]}"
echo ""

# Process each module
for module in "${MODULES[@]}"; do
    run_module "$module"
    
    # If module fails, exit the script
    if [ $? -ne 0 ]; then
        echo "Stopping due to error in module $module"
        exit 1
    fi
done

echo "================================================"
echo "All modules processed successfully!"
echo "================================================"