#!/bin/bash

# Noah-MP Spin-up Automation Script with Auto-Resume and CPU Temperature Control
# Usage: ./spinup_noahmp.sh [number_of_spinup_cycles]

# =============================================================================
# Configuration
# =============================================================================
SPINUP_CYCLES=${1:-20}  # Default 20 cycles if not specified
BASE_DIR="/land1/user/gychoi/LIS/run/run_NoahMP5.0"
RUN_DIR="$BASE_DIR/spin"
OUTPUT_BASE="$BASE_DIR/OUTPUT_era5/SURFACEMODEL"
SPIN_OUTPUT="$OUTPUT_BASE/spin"
CONFIG_FILE="$RUN_DIR/lis.config.era5"
RUN_SCRIPT="$RUN_DIR/run_LIS_noahmp36"

# Original restart file name pattern
RESTART_BASE="LIS_RST_NoahMP50_201101010000.d01"
RESTART_EXT=".nc"

# CPU Temperature Control Settings
MAX_CPU_TEMP=80          # Maximum allowed CPU temperature (°C)
TARGET_CPU_TEMP=60       # Target CPU temperature for optimal performance (°C)
MIN_SLEEP_TIME=30        # Minimum sleep time between cycles (seconds)
MAX_SLEEP_TIME=300       # Maximum sleep time for cooling (seconds)
TEMP_CHECK_INTERVAL=10   # Interval for temperature checks during cooling (seconds)

# Global variables for function communication
LAST_CYCLE=0
START_CYCLE=1

# =============================================================================
# CPU Temperature and System Monitoring Functions
# =============================================================================

get_cpu_temperature() {
    local temp=0
    
    # Try multiple methods to get CPU temperature
    if command -v sensors >/dev/null 2>&1; then
        # Method 1: lm-sensors (most reliable)
        temp=$(sensors 2>/dev/null | grep -i "package id 0\|core 0\|cpu" | head -1 | grep -o '+[0-9]\+\.[0-9]\+°C' | head -1 | sed 's/+//g' | sed 's/°C//g')
        
        # If package id 0 not found, try other patterns
        if [[ -z "$temp" ]]; then
            temp=$(sensors 2>/dev/null | grep -i "temp1\|cpu" | head -1 | grep -o '+[0-9]\+\.[0-9]\+°C' | head -1 | sed 's/+//g' | sed 's/°C//g')
        fi
    fi
    
    # Method 2: Try /sys/class/thermal if sensors failed
    if [[ -z "$temp" || "$temp" == "0" ]]; then
        if [[ -f /sys/class/thermal/thermal_zone0/temp ]]; then
            local temp_millic=$(cat /sys/class/thermal/thermal_zone0/temp 2>/dev/null)
            if [[ -n "$temp_millic" && "$temp_millic" -gt 0 ]]; then
                temp=$(echo "scale=1; $temp_millic / 1000" | bc 2>/dev/null)
            fi
        fi
    fi
    
    # Method 3: Try /proc/cpuinfo or other thermal zones
    if [[ -z "$temp" || "$temp" == "0" ]]; then
        # Try other thermal zones
        for zone in /sys/class/thermal/thermal_zone*/temp; do
            if [[ -f "$zone" ]]; then
                local temp_millic=$(cat "$zone" 2>/dev/null)
                if [[ -n "$temp_millic" && "$temp_millic" -gt 10000 && "$temp_millic" -lt 150000 ]]; then
                    temp=$(echo "scale=1; $temp_millic / 1000" | bc 2>/dev/null)
                    break
                fi
            fi
        done
    fi
    
    # Return temperature as integer (remove decimal part)
    if [[ -n "$temp" && "$temp" != "0" ]]; then
        echo "${temp%.*}"
    else
        echo "0"
    fi
}

get_cpu_load() {
    local load_1min=$(uptime | awk -F'load average:' '{print $2}' | awk '{print $1}' | sed 's/,//')
    echo "${load_1min%.*}"  # Return as integer
}

get_memory_usage() {
    local mem_percent=$(free | grep Mem | awk '{printf("%.0f", ($3/$2) * 100.0)}')
    echo "$mem_percent"
}

log_system_stats() {
    local temp=$(get_cpu_temperature)
    local load=$(get_cpu_load)
    local mem=$(get_memory_usage)
    
    log_message "System Status - CPU: ${temp}°C, Load: ${load}, Memory: ${mem}%"
}

calculate_dynamic_sleep() {
    local current_temp=$(get_cpu_temperature)
    local current_load=$(get_cpu_load)
    local sleep_time=$MIN_SLEEP_TIME
    
    # Temperature-based adjustment
    if [[ "$current_temp" -gt "$MAX_CPU_TEMP" ]]; then
        sleep_time=$MAX_SLEEP_TIME
        # stderr로 출력하여 함수 반환값과 분리
        log_message "High temperature detected (${current_temp}°C), extending sleep to ${sleep_time}s" >&2
    elif [[ "$current_temp" -gt "$TARGET_CPU_TEMP" ]]; then
        # Linear interpolation between MIN and MAX sleep times
        local temp_excess=$((current_temp - TARGET_CPU_TEMP))
        local temp_range=$((MAX_CPU_TEMP - TARGET_CPU_TEMP))
        local sleep_increase=$(echo "scale=0; ($temp_excess * ($MAX_SLEEP_TIME - $MIN_SLEEP_TIME)) / $temp_range" | bc)
        sleep_time=$((MIN_SLEEP_TIME + sleep_increase))
        # stderr로 출력하여 함수 반환값과 분리
        log_message "Elevated temperature (${current_temp}°C), adjusting sleep to ${sleep_time}s" >&2
    fi
    
    # Load-based adjustment (additional factor)
    local cpu_cores=$(nproc)
    if [[ "$current_load" -gt $((cpu_cores * 2)) ]]; then
        sleep_time=$((sleep_time + 60))  # Add 1 minute for high load
        # stderr로 출력하여 함수 반환값과 분리
        log_message "High system load detected (${current_load}), adding extra cooling time" >&2
    fi
    
    # Ensure sleep time is within bounds
    if [[ "$sleep_time" -lt "$MIN_SLEEP_TIME" ]]; then
        sleep_time=$MIN_SLEEP_TIME
    elif [[ "$sleep_time" -gt "$MAX_SLEEP_TIME" ]]; then
        sleep_time=$MAX_SLEEP_TIME
    fi
    
    # 순수한 숫자만 반환 (로그 메시지는 stderr로 이미 출력됨)
    echo "$sleep_time"
}

wait_for_cooling() {
    local current_temp=$(get_cpu_temperature)
    
    if [[ "$current_temp" -eq 0 ]]; then
        log_message "Warning: Cannot read CPU temperature, using fixed sleep time"
        sleep $MIN_SLEEP_TIME
        return 0
    fi
    
    log_message "Current CPU temperature: ${current_temp}°C"
    
    # If temperature is acceptable, use dynamic sleep
    if [[ "$current_temp" -le "$MAX_CPU_TEMP" ]]; then
        local sleep_time=$(calculate_dynamic_sleep)
        log_message "Cooling for ${sleep_time} seconds..."
        sleep "$sleep_time"
        return 0
    fi
    
    # If temperature is too high, actively cool down
    log_message "CPU temperature too high (${current_temp}°C > ${MAX_CPU_TEMP}°C), active cooling..."
    
    local cooling_start=$(date +%s)
    while [[ "$current_temp" -gt "$TARGET_CPU_TEMP" ]]; do
        log_message "Cooling... Current temperature: ${current_temp}°C (target: ≤${TARGET_CPU_TEMP}°C)"
        sleep "$TEMP_CHECK_INTERVAL"
        current_temp=$(get_cpu_temperature)
        
        # Safety check - don't cool forever
        local cooling_time=$(( $(date +%s) - cooling_start ))
        if [[ "$cooling_time" -gt 1800 ]]; then  # 30 minutes max
            log_message "Warning: Extended cooling time (${cooling_time}s). Proceeding with current temperature: ${current_temp}°C"
            break
        fi
    done
    
    local final_temp=$(get_cpu_temperature)
    local total_cooling_time=$(( $(date +%s) - cooling_start ))
    log_message "Cooling completed. Final temperature: ${final_temp}°C (cooled for ${total_cooling_time}s)"
}

intelligent_cycle_pause() {
    local cycle_num=$1
    local total_cycles=$2
    
    # Log system status before cooling
    log_system_stats
    
    # Different cooling strategies based on progress
    if [[ "$cycle_num" -eq "$total_cycles" ]]; then
        log_message "Final cycle completed - skipping pause"
        return 0
    elif [[ "$cycle_num" -gt $((total_cycles - 3)) ]]; then
        log_message "Near completion - using conservative cooling"
        wait_for_cooling
    else
        log_message "Standard inter-cycle cooling period"
        wait_for_cooling
    fi
    
    # Log system status after cooling
    log_system_stats
}

# =============================================================================
# Original Functions (with minimal modifications)
# =============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_system_requirements() {
    log_message "Checking system requirements..."
    
    # Check if bc is available for temperature calculations
    if ! command -v bc >/dev/null 2>&1; then
        log_message "Warning: 'bc' calculator not found. Installing or use manual calculations."
        # Try to continue without bc
    fi
    
    # Check temperature monitoring capability
    local test_temp=$(get_cpu_temperature)
    if [[ "$test_temp" -eq 0 ]]; then
        log_message "Warning: CPU temperature monitoring not available"
        log_message "         Will use fixed sleep times instead of dynamic cooling"
    else
        log_message "CPU temperature monitoring available: ${test_temp}°C"
        log_message "Temperature thresholds: Target=${TARGET_CPU_TEMP}°C, Max=${MAX_CPU_TEMP}°C"
    fi
    
    log_message "System requirements check completed"
}

check_directories() {
    log_message "Checking directory structure..."
    
    if [[ ! -d "$RUN_DIR" ]]; then
        log_message "ERROR: Run directory not found: $RUN_DIR"
        exit 1
    fi
    
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_message "ERROR: Config file not found: $CONFIG_FILE"
        exit 1
    fi
    
    if [[ ! -f "$RUN_SCRIPT" ]]; then
        log_message "ERROR: Run script not found: $RUN_SCRIPT"
        exit 1
    fi
    
    # Create output directories if they don't exist
    mkdir -p "$OUTPUT_BASE"
    mkdir -p "$SPIN_OUTPUT"
    
    log_message "Directory structure verified."
}

find_last_completed_cycle() {
    local max_cycle=0
    local found_files=()
    
    # Look for restart files with pattern: RESTART_BASE_N.nc
    for file in "$OUTPUT_BASE"/${RESTART_BASE}_*${RESTART_EXT}; do
        if [[ -f "$file" ]]; then
            # Extract cycle number from filename
            local basename=$(basename "$file")
            local cycle_num=$(echo "$basename" | sed "s/${RESTART_BASE}_\([0-9]\+\)${RESTART_EXT}/\1/")
            
            if [[ "$cycle_num" =~ ^[0-9]+$ ]]; then
                found_files+=("$file (cycle $cycle_num)")
                if [[ $cycle_num -gt $max_cycle ]]; then
                    max_cycle=$cycle_num
                fi
            fi
        fi
    done
    
    # Also check for original restart file (cycle 0)
    if [[ -f "$OUTPUT_BASE/${RESTART_BASE}${RESTART_EXT}" ]]; then
        found_files+=("$OUTPUT_BASE/${RESTART_BASE}${RESTART_EXT} (initial/cycle 0)")
    fi
    
    if [[ ${#found_files[@]} -gt 0 ]]; then
        log_message "Found existing restart files:"
        for file in "${found_files[@]}"; do
            log_message "  - $file"
        done
    else
        log_message "No existing restart files found."
    fi
    
    # Use global variable instead of echo
    LAST_CYCLE=$max_cycle
}

determine_start_cycle() {
    log_message "Searching for existing restart files..."
    find_last_completed_cycle
    
    local start_cycle=$((LAST_CYCLE + 1))
    
    log_message "Last completed cycle: $LAST_CYCLE"
    
    if [[ $LAST_CYCLE -eq 0 ]]; then
        # Check if initial restart file exists
        if [[ ! -f "$OUTPUT_BASE/${RESTART_BASE}${RESTART_EXT}" ]]; then
            log_message "ERROR: No initial restart file found: $OUTPUT_BASE/${RESTART_BASE}${RESTART_EXT}"
            log_message "Please ensure the initial restart file exists before running spin-up."
            exit 1
        fi
        log_message "Starting fresh spin-up from cycle 1"
    else
        if [[ $start_cycle -gt $SPINUP_CYCLES ]]; then
            log_message "INFO: Requested cycles ($SPINUP_CYCLES) already completed."
            log_message "Last completed cycle: $LAST_CYCLE"
            log_message "If you want to run more cycles, specify a higher number."
            exit 0
        fi
        log_message "Will resume spin-up from cycle $start_cycle (after completed cycle $LAST_CYCLE)"
    fi
    
    # Use global variable instead of echo
    START_CYCLE=$start_cycle
}

backup_config() {
    local backup_file="${CONFIG_FILE}.backup_$(date +%Y%m%d_%H%M%S)"
    cp "$CONFIG_FILE" "$backup_file"
    log_message "Config file backed up to: $backup_file"
}

update_restart_file_in_config() {
    local cycle_num=$1
    local restart_filename
    
    if [[ $cycle_num -eq 1 ]]; then
        restart_filename="${RESTART_BASE}${RESTART_EXT}"
    else
        restart_filename="${RESTART_BASE}_$((cycle_num-1))${RESTART_EXT}"
    fi
    
    log_message "Updating config file with restart file: $restart_filename"
    
    # Update the restart file line in config
    sed -i "s|Noah-MP.5.0 restart file:.*|Noah-MP.5.0 restart file:                  ../OUTPUT_era5/SURFACEMODEL/$restart_filename|" "$CONFIG_FILE"
    
    # Verify the update
    local updated_line=$(grep "Noah-MP.5.0 restart file:" "$CONFIG_FILE")
    log_message "Updated line: $updated_line"
}

copy_and_rename_restart() {
    local cycle_num=$1
    local source_dir="$SPIN_OUTPUT/SURFACEMODEL/201101"
    local source_file="$source_dir/${RESTART_BASE}${RESTART_EXT}"
    local target_file="$OUTPUT_BASE/${RESTART_BASE}_${cycle_num}${RESTART_EXT}"
    
    if [[ -f "$source_file" ]]; then
        cp "$source_file" "$target_file"
        log_message "Copied restart file: $target_file"
    else
        log_message "ERROR: Source restart file not found: $source_file"
        return 1
    fi
}

preserve_last_years() {
    local cycle_num=$1
    local total_cycles=$2
    
    # Preserve last 3 cycles for equilibrium analysis
    local cycles_to_preserve=3
    local preserve_threshold=$((total_cycles - cycles_to_preserve + 1))
    
    if [[ $cycle_num -ge $preserve_threshold ]]; then
        local years_from_end=$((total_cycles - cycle_num + 1))
        local preserve_dir="$SPIN_OUTPUT/last_${years_from_end}"
        
        log_message "Preserving cycle $cycle_num results as last_${years_from_end} for equilibrium analysis"
        
        # Create preservation directory
        mkdir -p "$preserve_dir"
        
        # Copy entire SURFACEMODEL directory
        if [[ -d "$SPIN_OUTPUT/SURFACEMODEL" ]]; then
            cp -r "$SPIN_OUTPUT/SURFACEMODEL" "$preserve_dir/"
            log_message "Results preserved in: $preserve_dir"
        fi
        
        return 0  # Indicate that this cycle was preserved
    fi
    
    return 1  # Indicate that this cycle was not preserved
}

clean_spin_output() {
    local cycle_num=$1
    local was_preserved=${2:-false}
    
    if [[ "$was_preserved" == "true" ]]; then
        log_message "Cycle $cycle_num was preserved - cleaning only temporary files"
        # Only clean the original SURFACEMODEL directory, keep preserved copies
        if [[ -d "$SPIN_OUTPUT/SURFACEMODEL" ]]; then
            rm -rf "$SPIN_OUTPUT/SURFACEMODEL"
            log_message "Original spin output cleaned, preserved copies maintained."
        fi
    else
        log_message "Cleaning spin output directory..."
        if [[ -d "$SPIN_OUTPUT/SURFACEMODEL" ]]; then
            rm -rf "$SPIN_OUTPUT/SURFACEMODEL"
            log_message "Spin output cleaned."
        fi
    fi
}

run_model() {
    local cycle_num=$1
    log_message "Starting spin-up cycle $cycle_num/$SPINUP_CYCLES"
    
    # Log system status before model run
    log_system_stats
    
    cd "$RUN_DIR" || exit 1
    
    # Make sure the run script is executable
    chmod +x "$RUN_SCRIPT"
    
    # Run the model
    log_message "Executing: $RUN_SCRIPT"
    local start_time=$(date +%s)
    
    if ! ./"$(basename "$RUN_SCRIPT")"; then
        log_message "ERROR: Model run failed for cycle $cycle_num"
        return 1
    fi
    
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    log_message "Model run completed for cycle $cycle_num (runtime: ${runtime}s)"
    
    # Log system status after model run
    log_system_stats
}

check_initial_restart() {
    # This function is now integrated into determine_start_cycle()
    return 0
}

# =============================================================================
# Main Script
# =============================================================================

main() {
    log_message "Starting Noah-MP Spin-up Automation with Auto-Resume and CPU Temperature Control"
    log_message "Target spin-up cycles: $SPINUP_CYCLES"
    log_message "Base directory: $BASE_DIR"
    
    # Initial checks
    check_system_requirements
    check_directories
    backup_config
    
    # Determine where to start
    determine_start_cycle
    
    log_message "Starting from cycle: $START_CYCLE"
    
    if [[ $START_CYCLE -gt $SPINUP_CYCLES ]]; then
        log_message "All requested cycles already completed. Exiting."
        return 0
    fi
    
    # Main spin-up loop
    for ((cycle=START_CYCLE; cycle<=SPINUP_CYCLES; cycle++)); do
        log_message "=== Spin-up Cycle $cycle/$SPINUP_CYCLES ==="
        
        # Update config file with appropriate restart file
        update_restart_file_in_config $cycle
        
        # Clean previous spin output first (before running new cycle)
        if [[ $cycle -gt $START_CYCLE ]]; then
            # Check if previous cycle was preserved
            local prev_cycle=$((cycle - 1))
            local cycles_to_preserve=3
            local preserve_threshold=$((SPINUP_CYCLES - cycles_to_preserve + 1))
            
            if [[ $prev_cycle -ge $preserve_threshold ]]; then
                clean_spin_output $prev_cycle true
            else
                clean_spin_output $prev_cycle false
            fi
        fi
        
        # Run the model
        if ! run_model $cycle; then
            log_message "ERROR: Spin-up failed at cycle $cycle"
            log_message "You can resume from this point by running the script again."
            exit 1
        fi
        
        # Copy and rename the new restart file
        if ! copy_and_rename_restart $cycle; then
            log_message "ERROR: Failed to copy restart file for cycle $cycle"
            exit 1
        fi
        
        # Preserve results if this is one of the last 3 cycles
        preserve_last_years $cycle $SPINUP_CYCLES
        was_preserved=$?
        
        log_message "Cycle $cycle completed successfully"
        log_message "Restart file saved as: ${RESTART_BASE}_${cycle}${RESTART_EXT}"
        
        if [[ $was_preserved -eq 0 ]]; then
            local years_from_end=$((SPINUP_CYCLES - cycle + 1))
            log_message "Results preserved for equilibrium analysis as: last_${years_from_end}"
        fi
        
        # Intelligent cooling pause between cycles (replaces simple sleep 30)
        intelligent_cycle_pause $cycle $SPINUP_CYCLES
    done
    
    # Final summary with preserved analysis data
    log_message "=== Spin-up Automation Completed ==="
    log_message "Total cycles completed: $SPINUP_CYCLES"
    log_message "Final restart file: ${RESTART_BASE}_${SPINUP_CYCLES}${RESTART_EXT}"
    log_message ""
    log_message "EQUILIBRIUM ANALYSIS DATA:"
    log_message "Last 3 years preserved in: $SPIN_OUTPUT"
    log_message "  └── last_3/SURFACEMODEL/  (3rd from end - Year N-2)"
    log_message "  └── last_2/SURFACEMODEL/  (2nd from end - Year N-1)" 
    log_message "  └── last_1/SURFACEMODEL/  (Final year - Year N)"
    log_message ""
    log_message "TO CHECK EQUILIBRIUM:"
    log_message "Compare the following variables across last_3, last_2, last_1:"
    log_message "  • LAI seasonal patterns"
    log_message "  • Soil moisture in all layers (SoilMoist_tavg)"
    log_message "  • Carbon pools (LEAFC, STEMC, ROOTC, SOILC)"
    log_message "  • Water/energy fluxes (Evap, GPP, NEE)"
    log_message ""
    log_message "IF NOT IN EQUILIBRIUM:"
    log_message "Run additional cycles: $0 [higher_number_of_cycles]"
    log_message "The script will automatically resume from where it left off"
    log_message ""
    log_message "IF IN EQUILIBRIUM:"
    log_message "Use restart file: $OUTPUT_BASE/${RESTART_BASE}_${SPINUP_CYCLES}${RESTART_EXT}"
    log_message ""
    log_message "ANALYSIS TIP:"
    log_message "Use Python/NCO tools to compare annual means between preserved years:"
    log_message "  ncra last_1/SURFACEMODEL/201??/LIS_HIST_*.nc year1_mean.nc"
    log_message "  ncra last_2/SURFACEMODEL/201??/LIS_HIST_*.nc year2_mean.nc"
    log_message "  ncra last_3/SURFACEMODEL/201??/LIS_HIST_*.nc year3_mean.nc"
}

# =============================================================================
# Script execution
# =============================================================================

# Check if help is requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 [number_of_spinup_cycles]"
    echo ""
    echo "Noah-MP Spin-up Automation Script with Auto-Resume and CPU Temperature Control"
    echo ""
    echo "Arguments:"
    echo "  number_of_spinup_cycles  Number of spin-up cycles to run (default: 20)"
    echo ""
    echo "Features:"
    echo "  • Automatically detects completed cycles and resumes from where it left off"
    echo "  • Preserves last 3 cycles for equilibrium analysis"
    echo "  • Safe restart capability - no risk of losing previous work"
    echo "  • Dynamic CPU temperature monitoring and cooling"
    echo "  • Adaptive sleep times based on system load and temperature"
    echo ""
    echo "Temperature Control:"
    echo "  • Target temperature: ${TARGET_CPU_TEMP}°C"
    echo "  • Maximum temperature: ${MAX_CPU_TEMP}°C"
    echo "  • Sleep time range: ${MIN_SLEEP_TIME}-${MAX_SLEEP_TIME} seconds"
    echo ""
    echo "Examples:"
    echo "  $0 30    # Run 30 spin-up cycles (or resume if some already completed)"
    echo "  $0       # Run 20 spin-up cycles (default)"
    echo ""
    echo "Resume behavior:"
    echo "  If you previously ran 15 cycles and now specify 30, it will:"
    echo "  - Detect cycles 1-15 are already completed"
    echo "  - Resume from cycle 16 and run through cycle 30"
    exit 0
fi

# Validate input
if [[ $# -gt 1 ]]; then
    log_message "ERROR: Too many arguments. Use -h for help."
    exit 1
fi

if [[ $# -eq 1 && ! "$1" =~ ^[0-9]+$ ]]; then
    log_message "ERROR: Argument must be a positive integer."
    exit 1
fi

# Run the main function
main
